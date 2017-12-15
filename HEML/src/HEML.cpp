#include <Ciphertext.h>
#include <CZZ.h>
#include <NTL/BasicThreadPool.h>
#include <NTL/tools.h>
#include <NTL/ZZ.h>
#include <Params.h>
#include <Plaintext.h>
#include <Scheme.h>
#include <SecretKey.h>
#include <TimeUtils.h>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <unistd.h>
#include <ios>
#include <fstream>
#include <string>

#include "MemoryUsage.h"
#include "CipherGD.h"
#include "GD.h"

using namespace std;

/*
 * run: ./HEML trainfile(string) isYfirst(bool) iter(long) learnPortion(double) approx(long) isEncrypted(bool) testfile(string)
 *
 * example: ./HEML "../data/data103x1579.txt" 1 7 1 7 1
 * example: ./HEML "../data/1_training_data_csv" 1 7 1 7 1 "../data/1_testing_data_csv"
 *
 * parameters:
 * trainfile - path to train file
 * isYfirst - {0,1} y parameter first OR last
 * iter - number of iterations
 * learnPortion - portion of data used for learning (randomly chosen from sample set)
 * approx - {3,5,7} polynomial degree approximation used
 * isEncrypted - {0,1} encrypted OR unencrypted learning
 * testfile - path to test file (checks if number of arguments <= 7 then testfile = trainfile)
 *
 * current files that in data folder (filename isYfirst):
 * "../data/data5x500.txt" false
 * "../data/data9x1253.txt" false
 * "../data/data15x1500.txt" false
 * "../data/data16x101.txt" false
 * "../data/data27x148.txt" false
 * "../data/data51x653.txt" false
 * "../data/data67x216.txt" false
 * "../data/data103x1579.txt" true
 * "../data/1_training_data.csv" true
 * "../data/1_testing_data.csv" true
 *
 * FYI: approx 3 suggested iter: 4, 9, 18, 36, ...
 * FYI: approx 5 suggested iter: 3, 7, 14, 28, ...
 * FYI: approx 7 suggested iter: 3, 7, 14, 28, ...
 */
int main(int argc, char **argv) {
	string trainfile(argv[1]);
	bool isYfirst = atoi(argv[2]);
	long iter = atol(argv[3]);
	double learnPortion = atof(argv[4]);
	long approx = atol(argv[5]);
	bool isEncrypted = atoi(argv[6]);

	string testfile = argc > 7 ? string(argv[7]) : trainfile;

	TimeUtils timeutils;
	SetNumThreads(8);
	/*
	 * gammaDownCnst > 0 : gamma = gammaUpCnst / gammaDownCnst / learnDim -> constant gamma
	 * gammaDownCnst < 0 : gamma = gammaUpCnst / (i + |gammaDownCnst|) / learnDim -> decreasing gamma
	 */
	double gammaUpCnst = 1;
	double gammaDownCnst = -1;

	long factorDim = 0;
	long sampleDim = 0;
	long** xyData = GD::xyDataFromFile(trainfile, factorDim, sampleDim, isYfirst);

	long factorDimTest = 0;
	long sampleDimTest = 0;
	long** xyDataTest;
	if(argc > 7) {
		xyDataTest = GD::xyDataFromFile(testfile, factorDimTest, sampleDimTest, isYfirst);
	} else {
		xyDataTest = xyData;
		factorDimTest = factorDim;
		sampleDimTest = sampleDim;
	}
	long learnDim = (long)((double)sampleDim * learnPortion);

	cout << "iter: " << iter << endl;
	cout << "isEncrypted: " << isEncrypted << endl;
	cout << "approx: " << approx << endl;
	cout << "gammaUpCnst: " << gammaUpCnst << endl;
	cout << "gammaDownCnst: " << gammaDownCnst << endl;
	cout << "sampleDim: " << sampleDim << endl;
	cout << "factorDim: " << factorDim << endl;
	cout << "learnDim: " << learnDim << endl;

	long** xyDataLearn = GD::RandomxyDataLearn(xyData, learnDim, sampleDim, factorDim);

	double alpha0, alpha1, alpha2;
	double etaprev, eta;
	double gamma;

	alpha0 = 0.01;

	alpha1 = (1. + sqrt(1. + 4.0 * alpha0 * alpha0)) / 2.0;
	alpha2 = (1. + sqrt(1. + 4.0 * alpha1 * alpha1)) / 2.0;

	if(!isEncrypted) {
		double* vData = new double[factorDim];
		double* wData = new double[factorDim];
		for (long i = 0; i < factorDim; ++i) {
			double tmp = 0.0;
			for (long j = 0; j < learnDim; ++j) {
				tmp += xyDataLearn[j][i];
			}
			tmp /= learnDim;
			wData[i] = tmp;
			vData[i] = tmp;
		}

		for (long k = 0; k < iter; ++k) {
			etaprev = (1 - alpha0) / alpha1;
			eta = (1 - alpha1) / alpha2;
			gamma = gammaDownCnst > 0 ? gammaUpCnst / gammaDownCnst / learnDim : gammaUpCnst / (k - gammaDownCnst) / learnDim;

			GD::stepNLGD(xyDataLearn, wData, vData, factorDim, learnDim, gamma, eta);

			timeutils.start("check on train data");
			GD::check(xyData, wData, factorDim, sampleDim);
			double auctrain = GD::calcuateAUC(xyData, wData, factorDim, sampleDim, 100);
			cout << "auc train: " << auctrain << endl;
			timeutils.stop("check on train data");

			timeutils.start("check on test data");
			GD::check(xyDataTest, wData, factorDimTest, sampleDimTest);
			double auctest = GD::calcuateAUC(xyDataTest, wData, factorDimTest, sampleDimTest, 100);
			cout << "auc test: " << auctest << endl;
			timeutils.stop("check on test data");

			alpha0 = alpha1;
			alpha1 = alpha2;
			alpha2 = (1. + sqrt(1. + 4.0 * alpha1 * alpha1)) / 2.0;
		}
	} else {
		long fdimBits = (long)ceil(log2(factorDim));
		long ldimBits = (long)ceil(log2(learnDim));
		long wBits = 37;
		long xyBits = 37;
		long lBits = 5;
		long pBits = 18;
		long aBits = 2;
		long logQ = (approx == 3) ? (ldimBits + xyBits) + iter * (2 * wBits + xyBits + pBits + aBits) + lBits
				: (ldimBits + xyBits) + iter * (3 * wBits + xyBits + pBits + aBits) + lBits;
		long logN = Params::suggestlogN(80, logQ);
		long bBits = min(logN - 1 - ldimBits, fdimBits);
		long batch = 1 << bBits;

		long sBits = ldimBits + bBits;
		long slots =  1 << sBits;
		long cnum = (long)ceil((double)factorDim / batch);

		cout << "fdimBits: " << fdimBits << endl;
		cout << "ldimBits: " << ldimBits << endl;
		cout << "xyBits: " << xyBits << endl;
		cout << "wBits: " << wBits << endl;
		cout << "lBits: " << lBits << endl;
		cout << "pBits: " << pBits << endl;
		cout << "aBits: " << aBits << endl;
		cout << "bBits: " << bBits << endl;
		cout << "slots: " << slots << endl;
		cout << "cnum: " << cnum << endl;

		size_t currentPreSchemeSize = getCurrentRSS( ) / 1048576;
		size_t peakPreSchemeSize = getPeakRSS() / 1048576;
		cout << "Current Memory Usage Before Scheme Generation: " << currentPreSchemeSize << "MB"<< endl;
		cout << "Peak Memory Usage Before Scheme Generation: " << peakPreSchemeSize << "MB"<< endl;

		timeutils.start("Scheme generating...");
		Params params(logN, logQ);
		Context context(params);
		SecretKey secretKey(params);
		Scheme scheme(secretKey, context);
		CipherGD cipherGD(scheme, secretKey);
		scheme.addLeftRotKeys(secretKey);
		scheme.addRightRotKeys(secretKey);
		timeutils.stop("Scheme generation");

		size_t currentAfterSchemeSize = getCurrentRSS( ) / 1048576;
		size_t peakAfterSchemeSize = getPeakRSS() / 1048576;
		cout << "Current Memory Usage After Scheme Generation: " << currentAfterSchemeSize << "MB"<< endl;
		cout << "Peak Memory Usage After Scheme Generation: " << peakAfterSchemeSize << "MB"<< endl;
		cout << "Scheme Size is approximately " << currentAfterSchemeSize - currentPreSchemeSize << "MB" << endl;

		cout << "HEAAN PARAMETER logQ: " << logQ << endl;
		cout << "HEAAN PARAMETER logN: " << logN << endl;
		cout << "HEAAN PARAMETER h: " << params.h << endl;
		cout << "HEAAN PARAMETER sigma: " << params.sigma << endl;

		timeutils.start("Polynomial generating...");
		ZZX poly = cipherGD.generateAuxPoly(slots, batch, pBits);
		timeutils.stop("Polynomial generation");

		Ciphertext* cxyData = new Ciphertext[cnum];
		timeutils.start("Encrypting xyData...");
		cipherGD.encxyData(cxyData, xyDataLearn, slots, factorDim, learnDim, batch, cnum, xyBits);
		timeutils.stop("xyData encryption");

		Ciphertext* cwData = new Ciphertext[cnum];
		Ciphertext* cvData = new Ciphertext[cnum];

		timeutils.start("Encrypting wData and vData...");
		cipherGD.encwData(cwData, cxyData, cnum, sBits, ldimBits, bBits, xyBits, wBits);
		for (long i = 0; i < cnum; ++i) {
			cvData[i] = cwData[i];
		}
		timeutils.stop("wData and vData encryption");

		size_t currentAfterCipherSize = getCurrentRSS( ) / 1048576;
		size_t peakAfterCipherSize = getPeakRSS() / 1048576;
		cout << "Current Memory Usage After Ciphertexts Generation: " << currentAfterCipherSize << "MB"<< endl;
		cout << "Peak Memory Usage After Ciphertexts Generation: " << peakAfterCipherSize << "MB"<< endl;
		cout << "Total Ciphertexts Size is approximately " << currentAfterCipherSize - currentAfterSchemeSize << "MB" << endl;

		double* dwData = new double[factorDim];
		for (long k = 0; k < iter; ++k) {

			cout << " !!! START " << k + 1 << " ITERATION !!! " << endl;

			etaprev = (1 - alpha0) / alpha1;
			eta = (1 - alpha1) / alpha2;
			gamma = gammaDownCnst > 0 ? gammaUpCnst / gammaDownCnst / learnDim : gammaUpCnst / (k - gammaDownCnst) / learnDim;

			cout << "cwData logq: " << cwData[0].logq << endl;
			timeutils.start("Encrypting NLGD step...");
			cipherGD.encNLGDiteration(approx, cxyData, cwData, cvData, poly, cnum, gamma, eta, etaprev, sBits, bBits, xyBits, wBits, pBits, aBits);
			timeutils.stop("NLGD step ");
			cout << "cwData logq: " << cwData[0].logq << endl;

			alpha0 = alpha1;
			alpha1 = alpha2;
			alpha2 = (1. + sqrt(1. + 4.0 * alpha1 * alpha1)) / 2.0;

			cout << " !!! STOP " << k + 1 << " ITERATION !!! " << endl;

			size_t currentAfterIterSize = getCurrentRSS( ) / 1048576;
			size_t peakAfterIterSize = getPeakRSS() / 1048576;
			cout << "Current Memory Usage After Iteration " << (k+1) << ": " << currentAfterIterSize << "MB"<< endl;
			cout << "Peak Memory Usage After Iteration " << (k+1) << ": " << peakAfterIterSize << "MB"<< endl;

			timeutils.start("decrypting");
			cipherGD.decwData(dwData, cwData, factorDim, batch, cnum, wBits);
			timeutils.stop("decrypting");

			timeutils.start("check on train data");
			GD::check(xyData, dwData, factorDim, sampleDim);
			double auctrain = GD::calcuateAUC(xyData, dwData, factorDim, sampleDim, 100);
			cout << "auc train: " << auctrain << endl;
			timeutils.stop("check on train data");

			timeutils.start("check on test data");
			GD::check(xyDataTest, dwData, factorDimTest, sampleDimTest);
			double auctest = GD::calcuateAUC(xyDataTest, dwData, factorDimTest, sampleDimTest, 100);
			cout << "auc test: " << auctest << endl;
			timeutils.stop("check on test data");

			size_t currentAfterDecSize = getCurrentRSS( ) / 1048576;
			size_t peakAfterDecSize = getPeakRSS() / 1048576;
			cout << "Current Memory Usage After Decryption " << (k+1) << ": " << currentAfterDecSize << "MB"<< endl;
			cout << "Peak Memory Usage After Decryption " << (k+1) << ": " << peakAfterDecSize << "MB"<< endl;
		}
		delete[] dwData;
	}
	size_t currentAfterDecSize = getCurrentRSS( ) / 1048576;
	size_t peakAfterDecSize = getPeakRSS() / 1048576;
	cout << "Current Memory Usage Final: " << currentAfterDecSize << "MB"<< endl;
	cout << "Peak Memory Usage Final: " << peakAfterDecSize << "MB"<< endl;

	return 0;
}
