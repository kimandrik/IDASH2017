#include <iostream>
#include <string>
#include <cmath>
#include <chrono>

#include <NTL/BasicThreadPool.h>

#include "Ciphertext.h"
#include "Context.h"
#include "NTL/ZZX.h"
#include "Scheme.h"
#include "SecretKey.h"
#include "TimeUtils.h"

#include "CipherGD.h"
#include "GD.h"
#include "TestGD.h"

using namespace std;
using namespace NTL;

chrono::high_resolution_clock::time_point t11;
chrono::high_resolution_clock::time_point t12;
chrono::high_resolution_clock::time_point t21;
chrono::high_resolution_clock::time_point t22;

void HEML(double** zData, long factorDim, long sampleDim, long numIter, long kdeg, double learningRate);

void ML(double** zData, long factorDim, long sampleDim, long numIter, long kdeg, double learningRate);

// > ./SEER 20 3 1.0 0 
// AUROC ~ 0.72

int main(int argc, char **argv) {
	SetNumThreads(1);
	string trainfile("data/SEER.csv");
	long numIter = atol(argv[1]);
	long kdeg = atol(argv[2]);
	double learningRate = atof(argv[3]);
	bool isEncrypted = atol(argv[4]);
	//-----------------------------------------
	long sampleDim = 0, factorDim = 0;
	double** zData = GD::zDataFromFile(trainfile, factorDim, sampleDim, false);
	GD::shuffleZData(zData, factorDim, sampleDim);
	GD::normalizeZData(zData, factorDim, sampleDim);
	//-----------------------------------------
	if(isEncrypted) {
		cout << "*** Start HEML with data " << trainfile << endl;
		HEML(zData, factorDim, sampleDim, numIter, kdeg, learningRate);
	} else {
		cout << "*** Start ML with data " << trainfile << endl;
		ML(zData, factorDim, sampleDim, numIter, kdeg, learningRate);
	}

	return 0;
}

void HEML(double** zData, long factorDim, long sampleDim, long numIter, long kdeg, double learningRate) {
	//-----------------------------------------
	long numIterPerBoot = 4;
	long numBoot = numIter / numIterPerBoot;
	double gammaUp = learningRate;
	double gammaDown = learningRate;
	bool isInitZero = false;
	//-----------------------------------------
	long sampleDimTest = (long)(0.3 * (double)sampleDim);
	long sampleDimTrain = sampleDim - sampleDimTest;
	long fdimBits = (long)ceil(log2(factorDim));
	long sdimBits = (long)ceil(log2(sampleDim));
	//-----------------------------------------
	long wBits = 30;
	long pBits = 15;
	long lBits = 5;
	long aBits = 3;
	long kBits = (long)ceil(log2(kdeg));
	//-----------------------------------------
	long logQ = isInitZero ? (wBits + lBits) + numIterPerBoot * ((kBits + 1) * wBits + 2 * pBits + aBits) :
			(sdimBits + wBits + lBits) + numIterPerBoot * ((kBits + 1) * wBits + 2 * pBits + aBits);
	long logT = 3;
	long logI = 4;
	long logq = wBits + 5;
	long bitForBoot = logT + 1 + 3 * (logq + logI) + (logI + logT) * (logq + logI) + logq + 2 * logI + 100;
	long logQBoot = logQ + bitForBoot;
	long logN = TestGD::suggestLogN(80, logQBoot);
	long bBits = min(logN - 1 - sdimBits, fdimBits);
	long batch = 1 << bBits;
	long sBits = sdimBits + bBits;
	long slots =  1 << sBits;
	long cnum = (long)ceil((double)factorDim / batch);
	//-----------------------------------------
	cout << "batch = " << batch << ", slots = " << slots << ", cnum = " << cnum << endl;
	cout << "HEAAN PARAMETER logQ: " << logQ << endl;
	cout << "HEAAN PARAMETER logQBoot: " << logQBoot << endl;
	cout << "HEAAN PARAMETER logN: " << logN << endl;
	//-----------------------------------------
	TimeUtils timeutils;
	timeutils.start("Scheme generating...");
	Context context(logN, logQBoot);
	SecretKey secretKey(logN);
	Scheme scheme(secretKey, context);
	scheme.addLeftRotKeys(secretKey);
	scheme.addRightRotKeys(secretKey);
	scheme.addConjKey(secretKey);
	scheme.addBootKey(secretKey, bBits, logq + logI);
	timeutils.stop("Scheme generation");
	CipherGD cipherGD(scheme, secretKey);
	//-----------------------------------------
	timeutils.start("Polynomial generating...");
	ZZX poly = cipherGD.generateAuxPoly(slots, batch, pBits);
	timeutils.stop("Polynomial generation");

	double* pwData = new double[factorDim];
	double* pvData = new double[factorDim];

	double* cwData = new double[factorDim];
	double* cvData = new double[factorDim];

	Ciphertext* encZData = new Ciphertext[cnum];
	Ciphertext* encWData = new Ciphertext[cnum];
	Ciphertext* encVData = new Ciphertext[cnum];

	double **zDataTrain, **zDataTest;

	zDataTrain = new double*[sampleDimTrain];
	zDataTest = new double*[sampleDimTest];

	double avgAUCENC = 0.;
	double avgAUCTRUE = 0.;

	for (long i = 0; i < sampleDimTest; ++i) zDataTest[i] = zData[i];
	for (long i = 0; i < sampleDimTrain; ++i) zDataTrain[i] = zData[sampleDimTest + i];

	timeutils.start("Encrypting zData");
	cipherGD.encZData(encZData, zDataTrain, slots, factorDim, sampleDimTrain, batch, cnum, wBits, logQ);
	timeutils.stop("Encrypting zData");

	timeutils.start("Encrypting wData and vData");
	cipherGD.encWVDataZero(encWData, encVData, cnum, slots, wBits, logQ);
	timeutils.stop("Encrypting wData and vData");

	GD::initialWDataVDataZero(pwData, pvData, factorDim);

	double alpha0, alpha1, eta, gamma, auctrain, auctest;
	
	alpha0 = 0.01;
	alpha1 = (1. + sqrt(1. + 4.0 * alpha0 * alpha0)) / 2.0;

	t11 = chrono::high_resolution_clock::now();
	for(long i = 0; i < numBoot; i++) {
		Ciphertext* encZDataCOPY = new Ciphertext[cnum];
		for(long j = 0; j < cnum; j++) encZDataCOPY[j] = encZData[j];
		for(long j = 0; j < numIterPerBoot; j++) {
			long iter = i * numIterPerBoot + j;
			cout << " !!! START " << iter + 1 << " ITERATION !!! " << endl;
			
			eta = (1 - alpha0) / alpha1;
			gamma = gammaDown > 0 ? gammaUp / gammaDown / sampleDimTrain : gammaUp / (iter - gammaDown) / sampleDimTrain;
				
			cout << "encWData.logq before: " << encWData[0].logq << endl;
			timeutils.start("Enc NLGD");
			cipherGD.encNLGDiteration(kdeg, encZDataCOPY, encWData, encVData, poly, cnum, gamma, eta, sBits, bBits, wBits, pBits, aBits);
			timeutils.stop("Enc NLGD");
			cout << "encWData.logq after: " << encWData[0].logq << endl;
			
			GD::plainNLGDiteration(kdeg, zDataTrain, pwData, pvData, factorDim, sampleDimTrain, gamma, eta);

			cout << "----ENCRYPTED-----" << endl;
			cipherGD.decWData(cwData, encWData, factorDim, batch, cnum, wBits);
			GD::calculateAUC(zDataTest, cwData, factorDim, sampleDimTest);
			cout << "------------------" << endl;

			cout << "-------PLAIN-------" << endl;
			GD::calculateAUC(zDataTest, pwData, factorDim, sampleDimTest);
			cout << "------------------" << endl;

			alpha0 = alpha1;
			alpha1 = (1. + sqrt(1. + 4.0 * alpha0 * alpha0)) / 2.0;
			cout << " !!! STOP " << iter + 1 << " ITERATION !!! " << endl;
		}
		for(long j = 0; j < cnum; j++) {
			t21 = chrono::high_resolution_clock::now();
			scheme.bootstrapAndEqual(encWData[j], logq, logQBoot, logT, logI);
			t22 = chrono::high_resolution_clock::now();
			cout << "Bootstrapping Time: " << chrono::duration_cast<chrono::minutes>(t22 - t21).count() << " min" << endl;
		}
		for(long j = 0; j < cnum; j++) {
			t21 = chrono::high_resolution_clock::now();
			scheme.bootstrapAndEqual(encVData[j], logq, logQBoot, logT, logI);
			t22 = chrono::high_resolution_clock::now();
			cout << "Bootstrapping Time: " << chrono::duration_cast<chrono::minutes>(t22 - t21).count() << " min" << endl;
		}
	}
	t12 = chrono::high_resolution_clock::now();
	cout << "Total HEML Time: " << chrono::duration_cast<chrono::minutes>(t12 - t11).count() << " min" << endl;
}

void ML(double** zData, long factorDim, long sampleDim, long numIter, long kdeg, double learningRate) {
	//-----------------------------------------
	double gammaUp = learningRate;
	double gammaDown = learningRate;
	bool isInitZero = false;
	long fold = 3;
	//-----------------------------------------
	long sampleDimTest = sampleDim / fold;
	long sampleDimTrain = sampleDim - sampleDimTest;
	long fdimBits = (long)ceil(log2(factorDim));
	long sdimBits = (long)ceil(log2(sampleDim));
	//-----------------------------------------
	long wBits = 30;
	long pBits = 15;
	long lBits = 5;
	long aBits = 3;
	long kBits = (long)ceil(log2(kdeg));
	//-----------------------------------------
	double **zDataTrain, **zDataTest;

	zDataTrain = new double*[sampleDimTrain];
	zDataTest = new double*[sampleDimTest];

	double* pwData = new double[factorDim];
	double* pvData = new double[factorDim];

	double* twData = new double[factorDim];
	double* tvData = new double[factorDim];

	double avgAUCPLAIN = 0.;
	double avgAUCTRUE = 0.;

	for (long fnum = 0; fnum < fold; ++fnum) {
		cout << " !!! START " << fnum + 1 << " FOLD !!! " << endl;

		for (long i = 0; i < sampleDimTest; ++i) {
			zDataTest[i] = zData[fnum * sampleDimTest + i];
		}
		for (long j = 0; j < fnum; ++j) {
			for (long i = 0; i < sampleDimTest; ++i) {
				zDataTrain[j * sampleDimTest + i] = zData[j * sampleDimTest + i];
			}
		}
		for (long i = (fnum + 1) * sampleDimTest; i < sampleDim; ++i) {
			zDataTrain[i - sampleDimTest] = zData[i];
		}

		if(isInitZero) {
			GD::initialWDataVDataZero(pwData, pvData, factorDim);
			GD::initialWDataVDataZero(twData, tvData, factorDim);
		} else {
			GD::initialWDataVDataAverage(pwData, pvData, zDataTrain, factorDim, sampleDimTrain);
			GD::initialWDataVDataAverage(twData, tvData, zDataTrain, factorDim, sampleDimTrain);
		}

		//-----------------------------------------

		double alpha0, alpha1, eta, gamma, auctrain, auctest, mse, nmse;

		alpha0 = 0.01;
		alpha1 = (1. + sqrt(1. + 4.0 * alpha0 * alpha0)) / 2.0;

		for (long iter = 0; iter < numIter; ++iter) {
			eta = (1 - alpha0) / alpha1;
			gamma = gammaDown > 0 ? gammaUp / gammaDown / sampleDimTrain : gammaUp / (iter - gammaDown) / sampleDimTrain;
			//-----------------------------------
			GD::plainNLGDiteration(kdeg, zDataTrain, pwData, pvData, factorDim, sampleDimTrain, gamma, eta);
			GD::trueNLGDiteration(zDataTrain, twData, tvData, factorDim, sampleDimTrain, gamma, eta);
			//-----------------------------------
			alpha0 = alpha1;
			alpha1 = (1. + sqrt(1. + 4.0 * alpha0 * alpha0)) / 2.0;
		}

		cout << "-------TRUE-------" << endl;
		avgAUCTRUE += GD::calculateAUC(zDataTest, twData, factorDim, sampleDimTest);
		cout << "------------------" << endl;

		cout << "-------PLAIN-------" << endl;
		avgAUCPLAIN += GD::calculateAUC(zDataTest, pwData, factorDim, sampleDimTest);
		cout << "------------------" << endl;

		cout << " !!! STOP " << fnum + 1 << " FOLD !!! " << endl;
		cout << "------------------" << endl;
	}
	cout << "average AUC (True) : " << avgAUCTRUE / to_double(fold) << endl;
	cout << "average AUC (Plain) : " << avgAUCPLAIN / to_double(fold) << endl;
}