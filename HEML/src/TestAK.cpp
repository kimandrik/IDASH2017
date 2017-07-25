#include "TestAK.h"

#include <NTL/BasicThreadPool.h>
#include <NTL/ZZ.h>
#include <Params.h>
#include <PubKey.h>
#include <Scheme.h>
#include <SchemeAlgo.h>
#include <SchemeAux.h>
#include <SecKey.h>
#include <TimeUtils.h>
#include <NumUtils.h>
#include <cmath>
#include <iostream>

#include "CipherSGD.h"
#include "SGD.h"

using namespace NTL;

//-----------------------------------------

void TestAK::testAK(long logN, long logl, long logp, long L) {
	cout << "!!! START TEST SGD !!!" << endl;
	//-----------------------------------------
	TimeUtils timeutils;
	SetNumThreads(8);
	//-----------------------------------------
	SGD sgd;

	string filename = "data103x1579.txt";
//	string filename = "data15x1500.txt";
//	string filename = "data5x500.txt";
//	string filename = "data9x1253.txt";

	long factorDim = 0;
	long sampleDim = 0;

	long** xyData = sgd.xyDataFromFile(filename, factorDim, sampleDim); //dim = 103, sampledim = 1579

	long sdimBits = (long)ceil(log2(sampleDim)); // log(1579) = 11
	long sampleDimPo2 = (1 << sdimBits); // 1579 -> 2048
	long learnDim = (1 << (sdimBits - 1)); //1024
//	long slots =  learnDim; // N /2
	long slots =  (1 << (logN-1)); // N /2
	long wBatch = slots / learnDim;
	long fdimBits = (long)ceil(log2(factorDim)); // log(103) = 7
	long factorDimPo2 = (1 << fdimBits); //103 -> 128

	cout << "factor dimension: " << factorDim << endl;
	cout << "factor dimension power of 2: " << factorDimPo2 << endl;

	cout << "sample dimension: " << sampleDim << endl;
	cout << "sample dimension power of 2: " << sampleDimPo2 << endl;

	cout << "learn dimension: " << learnDim << endl;

	cout << "slots: " << slots << endl;
	cout << "w batch: " << wBatch << endl;

	double** vData = new double*[wBatch];
	double** wData = new double*[wBatch];
	for (long l = 0; l < wBatch; ++l) {
		wData[l] = new double[factorDim];
		vData[l] = new double[factorDim];
		for (long i = 0; i < factorDim; ++i) {
			double tmp = (1.0 - 2.0 * (double)rand() / RAND_MAX) / 64.0;
//			double tmp = 0.0;
			wData[l][i] = tmp;
			vData[l][i] = tmp;
		}
	}

	long iter = 20;
	long enciter = 3;
	long totaliter = iter + enciter;

	double* alpha = new double[iter + 2];
	alpha[0] = 0.0;
	for (long i = 1; i < iter + 2; ++i) {
		alpha[i] = (1. + sqrt(1. + 4.0 * alpha[i-1] * alpha[i-1])) / 2.0;
	}

//	double lambda = 2.0;

	timeutils.start("sgd");
	for (long k = 0; k < iter; ++k) {

//		double gamma = 5.0 / (k+1);
		double gamma = 1.0 / 5.0;
		double eta = (1. - alpha[k+1]) / alpha[k+2];

		NTL_EXEC_RANGE(wBatch, first, last);
		for (long l = first; l < last; ++l) {
//			sgd.stepQGD(wData[l], xyData, gamma, lambda, factorDim, learnDim);
//			sgd.stepLGD(wData[l], xyData, gamma, lambda, factorDim, learnDim);
//			sgd.stepMLGD(wData[l], vData[l], xyData, gamma, lambda, factorDim, learnDim, eta);
			sgd.stepNLGD(wData[l], vData[l], xyData, gamma, factorDim, learnDim, eta);
		}
		NTL_EXEC_RANGE_END;
	}
	timeutils.stop("sgd");

	double* w = sgd.waverage(wData, wBatch, factorDim);

	sgd.check(w, xyData, factorDim, sampleDim);

	//-----------------------------------------
	Params params(logN, logl, logp, L);
	SecKey secretKey(params);
	PubKey publicKey(params, secretKey);
	SchemeAux schemeaux(params);
	Scheme scheme(params, publicKey, schemeaux);
	SchemeAlgo algo(scheme);
	CipherSGD csgd(scheme, algo, secretKey);
	//-----------------------------------------

	timeutils.start("Enc zdata");
	Cipher* cxyData = csgd.encxyData(xyData, slots, wBatch, factorDim, learnDim);
	timeutils.stop("Enc zdata");

	timeutils.start("Enc wdata");
	Cipher* cwData = csgd.encwData(wData, slots, wBatch, factorDim, learnDim);
	timeutils.stop("Enc wdata");

	//-----------------------------------------
	for (long k = iter; k < totaliter; ++k) {
		ZZ pgamma = ZZ(0);
		double lambda = 2.0;
		cout << k << endl;
		timeutils.start("Enc sgd step");
		csgd.encStepLGD(cxyData, cwData, pgamma, lambda, slots, wBatch, factorDim, learnDim);
		timeutils.stop("Enc sgd step");
	}

	timeutils.start("Enc w out");
	Cipher* cw = csgd.encwaverage(cwData, wBatch, factorDim);
	timeutils.stop("Enc w out");

	timeutils.start("Dec w");
	double* dw = csgd.decw(secretKey, cw, factorDim);
	timeutils.stop("Dec w");

	sgd.check(dw, xyData, factorDim, sampleDim);
	//-----------------------------------------
	cout << "!!! END TEST SGD !!!" << endl;
}
