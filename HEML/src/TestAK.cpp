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

#include "CipherGD.h"
#include "GD.h"

using namespace NTL;

//-----------------------------------------

void TestAK::testAK(long logN, long logl, long logp, long L) {
	cout << "!!! START TEST SGD !!!" << endl;
	//-----------------------------------------
	TimeUtils timeutils;
	SetNumThreads(8);
	//-----------------------------------------
	GD sgd;

	string filename = "data103x1579.txt";
//	string filename = "data15x1500.txt";
//	string filename = "data5x500.txt";
//	string filename = "data9x1253.txt";

	long factorDim = 0;
	long sampleDim = 0;

	long** xyData = sgd.xyDataFromFile(filename, factorDim, sampleDim); // factorDim=103, sampleDim=1579

	long sdimBits = (long)ceil(log2(sampleDim)); // 11
	long sampleDimPo2 = (1 << sdimBits); // 2048

	long fdimBits = (long)ceil(log2(factorDim)); // 7
	long factorDimPo2 = (1 << fdimBits); // 128

	long learnDim = sampleDim; // 1579
	long ldimBits = (long)ceil(log2(learnDim)); //11
	long learnDimPo2 = (1 << ldimBits); // 2048

	long slots =  (1 << (logN-1)); // 2^16
	long wBatch = slots / sampleDimPo2; // 32

	cout << "factorDim: " << factorDim << endl;
	cout << "factorDimPo2: " << factorDimPo2 << endl;

	cout << "sampleDim: " << sampleDim << endl;
	cout << "sampleDimPo2: " << sampleDimPo2 << endl;

	cout << "learnDim: " << learnDim << endl;
	cout << "learnDimPo2: " << learnDimPo2 << endl;

	cout << "slots: " << slots << endl;
	cout << "wBatch: " << wBatch << endl;

	double** vData = new double*[wBatch];
	double** wData = new double*[wBatch];
	for (long l = 0; l < wBatch; ++l) {
		wData[l] = new double[factorDim];
		vData[l] = new double[factorDim];
		for (long i = 0; i < factorDim; ++i) {
//			double tmp = (1.0 - 2.0 * (double)rand() / RAND_MAX) / 64.0; // take small initial values
			double tmp = 0.0;
			wData[l][i] = tmp;
			vData[l][i] = tmp;
		}
	}

	long iter = 10;
	long enciteradded = 3;
	long totaliter = iter + enciteradded;

	double* alpha = new double[totaliter + 2];
	alpha[0] = 0.0;
	for (long i = 1; i < totaliter + 2; ++i) {
		alpha[i] = (1. + sqrt(1. + 4.0 * alpha[i-1] * alpha[i-1])) / 2.0;
	}

	timeutils.start("sgd");
	for (long k = 0; k < iter; ++k) {

		double lambda = 0.0;
		double gamma = 0.001 / (1.0 + k);
		double eta = (1. - alpha[k+1]) / alpha[k+2];

		for (long l = 0; l < wBatch; ++l) {
//			sgd.stepLGD(xyData, wData[l], factorDim, learnDim, lambda, gamma);
			sgd.stepNLGD(xyData, wData[l], vData[l], factorDim, learnDim, lambda, gamma, eta);
		}
	}
	timeutils.stop("sgd");

	double* w = sgd.waverage(wData, factorDim, wBatch);

	sgd.check(xyData, w, factorDim, sampleDim);

	//-----------------------------------------
	Params params(logN, logl, logp, L);
	SecKey secretKey(params);
	PubKey publicKey(params, secretKey);
	SchemeAux schemeaux(params);
	Scheme scheme(params, publicKey, schemeaux);
	SchemeAlgo algo(scheme);
	CipherGD csgd(scheme, algo, secretKey);
	//-----------------------------------------

	timeutils.start("Enc zdata");
	Cipher* cxyData = csgd.encxyData(xyData, slots, factorDim, learnDim, wBatch);
	timeutils.stop("Enc zdata");

	timeutils.start("Enc wdata");
	Cipher* cwData = csgd.encwData(wData, slots, factorDim, learnDim, wBatch);
	timeutils.stop("Enc wdata");

	//-----------------------------------------
	for (long k = iter; k < totaliter; ++k) {

		double lambda = 0.0;
		double gamma = 1.0 / 5.0;
		double eta = (1. - alpha[k+1]) / alpha[k+2];

		timeutils.start("Enc sgd step");
		csgd.encStepLGD(cxyData, cwData, slots, factorDim, learnDim, wBatch, lambda, gamma);
		timeutils.stop("Enc sgd step");
	}

	timeutils.start("Enc w out");
	Cipher* cw = csgd.encwaverage(cwData, factorDim, wBatch);
	timeutils.stop("Enc w out");

	timeutils.start("Dec w");
	double* dw = csgd.decw(secretKey, cw, factorDim);
	timeutils.stop("Dec w");

	sgd.check(xyData, dw, factorDim, sampleDim);
	//-----------------------------------------
	cout << "!!! END TEST SGD !!!" << endl;
}
