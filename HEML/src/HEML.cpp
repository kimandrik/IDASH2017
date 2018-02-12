#include <Ciphertext.h>
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

	//-----------------------------------------
	string trainfile(argv[1]);
	bool isYfirst = atoi(argv[2]);
	long iter = atol(argv[3]);
	double learnPortion = atof(argv[4]);
	long approxDeg = atol(argv[5]);
	string testfile = argc > 6 ? string(argv[6]) : trainfile;
	//-----------------------------------------

	SetNumThreads(8);

	//-----------------------------------------

	long factorDimTrain = 0, sampleDimTrain = 0, factorDimTest = 0, sampleDimTest = 0;
	double **xyDataTrain, **xyDataTest;

	if(argc > 7) {
		xyDataTest = GD::xyDataFromFile(testfile, factorDimTest, sampleDimTest, isYfirst);
		xyDataTrain = GD::xyDataFromFile(trainfile, factorDimTrain, sampleDimTrain, isYfirst);
		if(factorDimTest != factorDimTrain) {
			invalid_argument("factor dimensions of learn and test datasets do not match");
		}
	} else {
		xyDataTest = GD::xyDataFromFile(testfile, factorDimTest, sampleDimTest, isYfirst);
		sampleDimTrain = (long)((double)sampleDimTest * learnPortion);
		factorDimTrain = factorDimTest;
		xyDataTrain = GD::RandomxyDataLearn(xyDataTest, sampleDimTrain, sampleDimTest, factorDimTest);
	}
	GD::normalizexyData2(xyDataTrain, xyDataTest, factorDimTrain, sampleDimTrain, sampleDimTest);

	//-----------------------------------------

	cout << "number of iterations: " << iter << endl;
	cout << "sigmoid approximation polynomial degree: " << approxDeg << endl;
	cout << "factorDimTrain: " << factorDimTrain << endl;
	cout << "sampleDimTrain: " << sampleDimTrain << endl;

	long fdimBits = (long)ceil(log2(factorDimTrain));
	long sdimBits = (long)ceil(log2(sampleDimTrain));

	long wBits = 30;
	long pBits = 20;
	long lBits = 5;
	long aBits = 3;

//	long logQ = (approxDeg == 3) ? (sdimBits + wBits + lBits) + iter * (3 * wBits + 2 * pBits + aBits)
//			: (sdimBits + wBits + lBits) + iter * (4 * wBits + 2 * pBits + aBits);

	long logQ = (approxDeg == 3) ? (wBits + lBits) + iter * (3 * wBits + 2 * pBits + aBits)
			: (wBits + lBits) + iter * (4 * wBits + 2 * pBits + aBits);

	long logN = Params::suggestlogN(80, logQ);
	long bBits = min(logN - 1 - sdimBits, fdimBits);
	long batch = 1 << bBits;
	long sBits = sdimBits + bBits;
	long slots =  1 << sBits;
	long cnum = (long)ceil((double)factorDimTrain / batch);

	cout << "fdimBits: " << fdimBits << endl;
	cout << "sdimBits: " << sdimBits << endl;
	cout << "wBits: " << wBits << endl;
	cout << "pBits: " << pBits << endl;
	cout << "lBits: " << lBits << endl;
	cout << "aBits: " << aBits << endl;
	cout << "bBits: " << bBits << endl;
	cout << "slots: " << slots << endl;
	cout << "cnum: " << cnum << endl;

	double* wData = new double[factorDimTrain];
	double* vData = new double[factorDimTrain];
	double* dwData = new double[factorDimTrain];
	double* dvData = new double[factorDimTrain];

//	long sdimPow = 1 << sdimBits;
//	for (long i = 0; i < factorDimTrain; ++i) {
//		double tmp = 0.0;
//		for (long j = 0; j < sampleDimTrain; ++j) {
//			tmp += xyDataTrain[j][i];
//		}
//		tmp /= sdimPow;
//		wData[i] = tmp;
//		vData[i] = tmp;
//	}

	for (long i = 0; i < factorDimTrain; ++i) {
		wData[i] = 0.0;
		vData[i] = 0.0;
	}

	//	size_t currentAfterSchemeSize = getCurrentRSS( ) >> 20;
	//	size_t peakAfterSchemeSize = getPeakRSS() >> 20;
	//	cout << "Current Memory Usage After Scheme Generation: " << currentAfterSchemeSize << "MB"<< endl;
	//	cout << "Peak Memory Usage After Scheme Generation: " << peakAfterSchemeSize << "MB"<< endl;

	Ciphertext* cxyData = new Ciphertext[cnum];
	Ciphertext* cwData = new Ciphertext[cnum];
	Ciphertext* cvData = new Ciphertext[cnum];

	//-----------------------------------------

	cout << "HEAAN PARAMETER logQ: " << logQ << endl;
	cout << "HEAAN PARAMETER logN: " << logN << endl;

	TimeUtils timeutils;
	timeutils.start("Scheme generating...");
	Params params(logN, logQ);
	Context context(params);
	SecretKey secretKey(params);
	Scheme scheme(secretKey, context);
	scheme.addLeftRotKeys(secretKey);
	scheme.addRightRotKeys(secretKey);
	timeutils.stop("Scheme generation");

	CipherGD cipherGD(scheme, secretKey);
	timeutils.start("Encrypting xyData...");
	cipherGD.encxyData(cxyData, xyDataTrain, slots, factorDimTrain, sampleDimTrain, batch, cnum, wBits);
	timeutils.stop("xyData encryption");

//	timeutils.start("Encrypting wData and vData...");
//	cipherGD.encwData(cwData, cxyData, cnum, sBits, bBits);
//	for (long i = 0; i < cnum; ++i) {
//		cvData[i] = cwData[i];
//	}
//	timeutils.stop("wData and vData encryption");

	timeutils.start("Encrypting wData and vData...");
	cipherGD.encwData0(cwData, cnum, slots, wBits);
	for (long i = 0; i < cnum; ++i) {
		cvData[i] = cwData[i];
	}
	timeutils.stop("wData and vData encryption");

	timeutils.start("Polynomial generating...");
	ZZX poly = cipherGD.generateAuxPoly(slots, batch, pBits);
	timeutils.stop("Polynomial generation");

	//-----------------------------------------

	double gammaUpCnst = 70;
	double gammaDownCnst = -5;
	double alpha0, alpha1, eta, gamma;
	double auctrain, auctest, mse, nmse;

	/*
	 * gammaDownCnst > 0 : gamma = gammaUpCnst / gammaDownCnst / learnDim -> constant gamma
	 * gammaDownCnst < 0 : gamma = gammaUpCnst / (i + |gammaDownCnst|) / learnDim -> decreasing gamma
	 */
	cout << "gammaUpCnst: " << gammaUpCnst << endl;
	cout << "gammaDownCnst: " << gammaDownCnst << endl;

	alpha0 = 0.01;
	alpha1 = (1. + sqrt(1. + 4.0 * alpha0 * alpha0)) / 2.0;

	for (long k = 0; k < iter; ++k) {
		cout << " !!! START " << k + 1 << " ITERATION !!! " << endl;
		//-----------------------------------------
		eta = (1 - alpha0) / alpha1;
		gamma = gammaDownCnst > 0 ? gammaUpCnst / gammaDownCnst / sampleDimTrain : gammaUpCnst / (k - gammaDownCnst) / sampleDimTrain;
		//-----------------------------------------

		cout << "cwData logq before: " << cwData[0].logq << endl;
		timeutils.start("Enc NLGD");
		cipherGD.encNLGDiteration(approxDeg, cxyData, cwData, cvData, poly, cnum, gamma, eta, sBits, bBits, wBits, pBits, aBits);
		timeutils.stop("Enc NLGD");
		cout << "cwData logq after: " << cwData[0].logq << endl;
		cipherGD.decwData(dwData, cwData, factorDimTrain, batch, cnum, wBits);
		timeutils.start("Encrypted check on train data");
		GD::check(xyDataTrain, dwData, factorDimTrain, sampleDimTrain);
		auctrain = GD::calculateAUC(xyDataTrain, dwData, factorDimTrain, sampleDimTrain);
		cout << "auc train: " << auctrain << endl;
		timeutils.stop("Encrypted check on train data");
		timeutils.start("Encrypted check on test data");
		GD::check(xyDataTest, dwData, factorDimTest, sampleDimTest);
		auctest = GD::calculateAUC(xyDataTest, dwData, factorDimTest, sampleDimTest);
		cout << "auc test: " << auctest << endl;
		timeutils.stop("Encrypted check on test data");

		//-----------------------------------------

//		GD::plainNLGDiteration(approxDeg, xyDataTrain, wData, vData, factorDimTrain, sampleDimTrain, gamma, eta);
//		timeutils.start("Plain check on train data");
//		GD::check(xyDataTrain, wData, factorDimTrain, sampleDimTrain);
//		auctrain = GD::calculateAUC(xyDataTrain, wData, factorDimTrain, sampleDimTrain);
//		cout << "auc train: " << auctrain << endl;
//		timeutils.stop("Plain check on train data");
//		timeutils.start("Plain check on test data");
//		GD::check(xyDataTest, wData, factorDimTest, sampleDimTest);
//		auctest = GD::calculateAUC(xyDataTest, wData, factorDimTest, sampleDimTest);
//		cout << "auc test: " << auctest << endl;
//		timeutils.stop("Plain check on test data");

		//-----------------------------------------

		GD::trueNLGDiteration(xyDataTrain, wData, vData, factorDimTrain, sampleDimTrain, gamma, eta);
		timeutils.start("True check on train data");
		GD::check(xyDataTrain, wData, factorDimTrain, sampleDimTrain);
		auctrain = GD::calculateAUC(xyDataTrain, wData, factorDimTrain, sampleDimTrain);
		cout << "auc train: " << auctrain << endl;
		timeutils.stop("True check on train data");
		timeutils.start("True check on test data");
		GD::check(xyDataTest, wData, factorDimTest, sampleDimTest);
		auctest = GD::calculateAUC(xyDataTest, wData, factorDimTest, sampleDimTest);
		cout << "auc test: " << auctest << endl;
		timeutils.stop("True check on test data");

		//-----------------------------------------

		mse = GD::calculateMSE(wData, dwData, factorDimTrain);
		nmse = GD::calculateNMSE(wData, dwData, factorDimTrain);
		cout << "mse: " << mse << endl;
		cout << "nmse: " << nmse << endl;

		//-----------------------------------------
		alpha0 = alpha1;
		alpha1 = (1. + sqrt(1. + 4.0 * alpha0 * alpha0)) / 2.0;
		//-----------------------------------------
		cout << " !!! STOP " << k + 1 << " ITERATION !!! " << endl;
	}
	return 0;
}
