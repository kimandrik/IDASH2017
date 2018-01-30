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
#include <StringUtils.h>

#include "MemoryUsage.h"
#include "CipherGD.h"
#include "GD.h"

using namespace std;

/*
 * run: ./HEML trainfile(string) isYfirst(bool) iter(long) learnPortion(double) approx(long) testfile(string)
 *
 * example: ./HEML "../data/data103x1579.txt" 1 7 1 7
 * example: ./HEML "../data/1_training_data_csv" 1 7 1 7 "../data/1_testing_data_csv"
 *
 * parameters:
 * trainfile - path to train file
 * isYfirst - {0,1} y parameter first OR last
 * iter - number of iterations
 * learnPortion - portion of data used for learning (randomly chosen from sample set)
 * approx - {3,5,7} polynomial degree approximation used
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
	long approxDeg = atol(argv[5]);

	string testfile = argc > 6 ? string(argv[6]) : trainfile;

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
	long factorDimTest = 0;
	long sampleDimTest = 0;

	double** xyData = GD::xyDataFromFile(trainfile, factorDim, sampleDim, isYfirst);
	double** xyDataTest;
	if(argc > 7) {
		xyDataTest = GD::xyDataFromFile(testfile, factorDimTest, sampleDimTest, isYfirst);
	} else {
		xyDataTest = xyData;
		factorDimTest = factorDim;
		sampleDimTest = sampleDim;
	}
	long learnDim = (long)((double)sampleDim * learnPortion);

	cout << "iter: " << iter << endl;
	cout << "approxDeg: " << approxDeg << endl;
	cout << "gammaUpCnst: " << gammaUpCnst << endl;
	cout << "gammaDownCnst: " << gammaDownCnst << endl;
	cout << "sampleDim: " << sampleDim << endl;
	cout << "factorDim: " << factorDim << endl;
	cout << "learnDim: " << learnDim << endl;

	double alpha0, alpha1, eta, gamma;
	double auctrain, auctest;
	alpha0 = 0.01;
	alpha1 = (1. + sqrt(1. + 4.0 * alpha0 * alpha0)) / 2.0;

	long fdimBits = (long)ceil(log2(factorDim));
	long ldimBits = (long)ceil(log2(learnDim));
	long wBits = 30;
	long pBits = 20;
	long lBits = 5;
	long aBits = 3;
	long logQ = (approxDeg == 3) ? (ldimBits + wBits + lBits) + iter * (3 * wBits + 2 * pBits + aBits)
			: (ldimBits + wBits + lBits) + iter * (4 * wBits + 2 * pBits + aBits);
	long logN = Params::suggestlogN(80, logQ);
	long bBits = min(logN - 1 - ldimBits, fdimBits);
	long batch = 1 << bBits;

	long sBits = ldimBits + bBits;
	long slots =  1 << sBits;
	long cnum = (long)ceil((double)factorDim / batch);

	cout << "fdimBits: " << fdimBits << endl;
	cout << "ldimBits: " << ldimBits << endl;
	cout << "wBits: " << wBits << endl;
	cout << "pBits: " << pBits << endl;
	cout << "lBits: " << lBits << endl;
	cout << "aBits: " << aBits << endl;
	cout << "bBits: " << bBits << endl;
	cout << "slots: " << slots << endl;
	cout << "cnum: " << cnum << endl;

	double** xyDataLearn = GD::RandomxyDataLearn(xyData, learnDim, sampleDim, factorDim);

	double* wData = new double[factorDim];
	double* vData = new double[factorDim];

	double* dwData = new double[factorDim];
	double* dvData = new double[factorDim];

	long ldimPow = 1 << ldimBits;
	for (long i = 0; i < factorDim; ++i) {
		double tmp = 0.0;
		for (long j = 0; j < learnDim; ++j) {
			tmp += xyDataLearn[j][i];
		}
		tmp /= ldimPow;
		wData[i] = tmp;
		vData[i] = tmp;
	}

	cout << "HEAAN PARAMETER logQ: " << logQ << endl;
	cout << "HEAAN PARAMETER logN: " << logN << endl;

	timeutils.start("Scheme generating...");
	Params params(logN, logQ);
	Context context(params);
	SecretKey secretKey(params);
	Scheme scheme(secretKey, context);
	CipherGD cipherGD(scheme, secretKey);
	scheme.addLeftRotKeys(secretKey);
	scheme.addRightRotKeys(secretKey);
	timeutils.stop("Scheme generation");

//	size_t currentAfterSchemeSize = getCurrentRSS( ) >> 20;
//	size_t peakAfterSchemeSize = getPeakRSS() >> 20;
//	cout << "Current Memory Usage After Scheme Generation: " << currentAfterSchemeSize << "MB"<< endl;
//	cout << "Peak Memory Usage After Scheme Generation: " << peakAfterSchemeSize << "MB"<< endl;

	Ciphertext* cxyData = new Ciphertext[cnum];
	Ciphertext* cwData = new Ciphertext[cnum];
	Ciphertext* cvData = new Ciphertext[cnum];

	timeutils.start("Encrypting xyData...");
	cipherGD.encxyData(cxyData, xyDataLearn, slots, factorDim, learnDim, batch, cnum, wBits);
	timeutils.stop("xyData encryption");

	timeutils.start("Encrypting wData and vData...");
	cipherGD.encwData(cwData, cxyData, cnum, sBits, bBits);
	for (long i = 0; i < cnum; ++i) {
		cvData[i] = cwData[i];
	}
	timeutils.stop("wData and vData encryption");

	timeutils.start("Polynomial generating...");
	ZZX poly = cipherGD.generateAuxPoly(slots, batch, pBits);
	timeutils.stop("Polynomial generation");

	for (long k = 0; k < iter; ++k) {
		cout << " !!! START " << k + 1 << " ITERATION !!! " << endl;
		eta = (1 - alpha0) / alpha1;
		gamma = gammaDownCnst > 0 ? gammaUpCnst / gammaDownCnst / learnDim : gammaUpCnst / (k - gammaDownCnst) / learnDim;

//		GD::trueNLGDiteration(xyData, wData, vData, factorDim, learnDim, gamma, eta);
		GD::plainNLGDiteration(approxDeg, xyData, wData, vData, factorDim, learnDim, gamma, eta);

		cout << "cwData logq: " << cwData[0].logq << endl;
		timeutils.start("Enc NLGD");
		cipherGD.encNLGDiteration(approxDeg, cxyData, cwData, cvData, poly, cnum, gamma, eta, sBits, bBits, wBits, pBits, aBits);
		timeutils.stop("Enc NLGD");
		cout << "cwData logq: " << cwData[0].logq << endl;

		alpha0 = alpha1;
		alpha1 = (1. + sqrt(1. + 4.0 * alpha0 * alpha0)) / 2.0;

		cout << " !!! STOP " << k + 1 << " ITERATION !!! " << endl;

		cipherGD.decwData(dwData, cwData, factorDim, batch, cnum, wBits);
		timeutils.start("Encrypted check on train data");
		GD::check(xyData, dwData, factorDim, sampleDim);
		auctrain = GD::calcuateAUC(xyData, dwData, factorDim, sampleDim, 100);
		cout << "auc train: " << auctrain << endl;
		timeutils.stop("Encrypted check on train data");

		timeutils.start("Encrypted check on test data");
		GD::check(xyDataTest, dwData, factorDimTest, sampleDimTest);
		auctest = GD::calcuateAUC(xyDataTest, dwData, factorDimTest, sampleDimTest, 100);
		cout << "auc test: " << auctest << endl;
		timeutils.stop("Encrypted check on test data");

		timeutils.start("Plain check on train data");
		GD::check(xyData, wData, factorDim, sampleDim);
		auctrain = GD::calcuateAUC(xyData, wData, factorDim, sampleDim, 100);
		cout << "auc train: " << auctrain << endl;
		timeutils.stop("Encrypted check on train data");

		timeutils.start("Plain check on test data");
		GD::check(xyDataTest, wData, factorDimTest, sampleDimTest);
		auctest = GD::calcuateAUC(xyDataTest, wData, factorDimTest, sampleDimTest, 100);
		cout << "auc test: " << auctest << endl;
		timeutils.stop("Plain check on test data");
	}
	return 0;
}
