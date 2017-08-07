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

void TestAK::testNLGDWB() {
	cout << "!!! START TEST NLGD WB !!!" << endl;
	//-----------------------------------------
	TimeUtils timeutils;
	SetNumThreads(8);
	//-----------------------------------------
	GD sgd;

	string filename = "data/data5x500.txt";     // false   415/500 done
//	string filename = "data/data9x1253.txt";    // false   775/1253 unclear
//	string filename = "data/data15x1500.txt";   // false   1270/1500 done
//	string filename = "data/data16x101.txt";    // false   101/101 done
//	string filename = "data/data27x148.txt";    // false   132/148 done
//	string filename = "data/data43x3247.txt";   // false   3182/3247
//	string filename = "data/data45x296.txt";    // false   257/296
//	string filename = "data/data51x653.txt";    // false   587/653
//	string filename = "data/data67x216.txt";    // false   216/216
//	string filename = "data/data103x1579.txt";  // true    1086/1579

	long factorDim = 0;
	long sampleDim = 0;

	long** xyData = sgd.xyDataFromFile(filename, factorDim, sampleDim, false);

	long sdimBits = (long)ceil(log2(sampleDim));
	long sampleDimPo2 = (1 << sdimBits);
	long fdimBits = (long)ceil(log2(factorDim));
	long factorDimPo2 = (1 << fdimBits);
	long learnDim = sampleDim;
	long ldimBits = (long)ceil(log2(learnDim));
	long learnDimPo2 = (1 << ldimBits);

	cout << "factorDim: " << factorDim << endl;
	cout << "fdimBits: " << fdimBits << endl;
	cout << "factorDimPo2: " << factorDimPo2 << endl;
	cout << "sampleDim: " << sampleDim << endl;
	cout << "sdimBits: " << sdimBits << endl;
	cout << "sampleDimPo2: " << sampleDimPo2 << endl;
	cout << "learnDim: " << learnDim << endl;
	cout << "ldimBits: " << ldimBits << endl;
	cout << "learnDimPo2: " << learnDimPo2 << endl;

	long wBatch = 1;
	long iter = fdimBits;
//	long iter = 5000;
	long logl = 5;
	long logp = 32;
	long L = 5 * iter + 1;
	long logN = Params::suggestlogN(80, logl, logp, L);
//	long logN = max(12, ldimBits);
	bool encrypted = false;
	long slots =  learnDimPo2 * wBatch;

	cout << "logl: " << logl << endl;
	cout << "logp: " << logp << endl;
	cout << "L: " << L << endl;
	cout << "logN: " << logN << endl;
	cout << "slots: " << slots << endl;
	cout << "wBatch: " << wBatch << endl;

	double** vData = new double*[wBatch];
	double** wData = new double*[wBatch];
	for (long l = 0; l < wBatch; ++l) {
		wData[l] = new double[factorDim];
		vData[l] = new double[factorDim];
		for (long i = 0; i < factorDim; ++i) {
			double tmp = (0.5 - 1.0 * (double)rand() / RAND_MAX) / factorDim;
//			double tmp = 0;
			wData[l][i] = tmp;
			vData[l][i] = tmp;
		}
	}

	double* alpha = new double[iter + 2]; // just constansts for Nesterov GD
	alpha[0] = 0.1;
	for (long i = 1; i < iter + 2; ++i) {
		alpha[i] = (1. + sqrt(1. + 4.0 * alpha[i-1] * alpha[i-1])) / 2.0;
	}

	if(!encrypted) {
		for (long k = 0; k < iter; ++k) {

			double gamma = 2.0 / learnDim / (1.0 + k);
			double eta = (1. - alpha[k+1]) / alpha[k+2];

			for (long l = 0; l < wBatch; ++l) {
				sgd.stepNLGD(xyData, wData[l], vData[l], factorDim, learnDim, gamma, eta);
				double* w = sgd.wsum(wData, factorDim, wBatch);
				sgd.check(xyData, w, factorDim, sampleDim);
			}
		}
	} else {
		//-----------------------------------------
		Params params(logN, logl, logp, L);
		SecKey secretKey(params);
		PubKey publicKey(params, secretKey);
		SchemeAux schemeaux(params);
		Scheme scheme(params, publicKey, schemeaux);
		SchemeAlgo algo(scheme);
		CipherGD csgd(scheme, algo, secretKey);
		//-----------------------------------------

		timeutils.start("Enc xyData");
		Cipher* cxyData = csgd.encxyDataWB(xyData, slots, factorDim, learnDim, wBatch);
		timeutils.stop("Enc xyData");

		timeutils.start("Enc wData");
		Cipher* cwData = csgd.encwDataWB(wData, slots, factorDim, learnDim, wBatch);
		timeutils.stop("Enc wData");

		Cipher* cvData = new Cipher[factorDim];
		for (long i = 0; i < factorDim; ++i) {cvData[i] = cwData[i];}

		for (long k = 0; k < iter; ++k) {
			double gamma = 2.0 / learnDim / (1.0 + k);
			double eta = (1. - alpha[k+1]) / alpha[k+2];

			timeutils.start("Enc sgd step");
			csgd.encStepNLGD7WB(cxyData, cwData, cvData, slots, factorDim, learnDim, wBatch, gamma, eta);
			timeutils.stop("Enc sgd step");
			csgd.debugcheck("c wData: ", secretKey, cwData, 5, wBatch);
		}

		timeutils.start("Enc w out");
		csgd.encwsumWB(cwData, factorDim, wBatch);
		timeutils.stop("Enc w out");

		timeutils.start("Dec w");
		double* dw = csgd.decWB(secretKey, cwData, factorDim);
		timeutils.stop("Dec w");

		sgd.check(xyData, dw, factorDim, sampleDim);
	}
	//-----------------------------------------
	cout << "!!! END TEST NLGD WB !!!" << endl;
}

void TestAK::testNLGDXYB() {
	cout << "!!! START TEST NLGD XYB !!!" << endl;
	//-----------------------------------------
	TimeUtils timeutils;
	SetNumThreads(8);
	//-----------------------------------------
	GD gd;

//	string filename = "data/data5x500.txt";     // false   415/500
//	string filename = "data/data9x1253.txt";    // false   775/1253 not good results
//	string filename = "data/data15x1500.txt";   // false   1270/1500
//	string filename = "data/data16x101.txt";    // false   101/101
//	string filename = "data/data27x148.txt";    // false   132/148
//	string filename = "data/data43x3247.txt";   // false   3182/3247
//	string filename = "data/data45x296.txt";    // false   257/296
//	string filename = "data/data51x653.txt";    // false   587/653
//	string filename = "data/data67x216.txt";    // false   216/216 slow convergence
	string filename = "data/data103x1579.txt";  // true    1086/1579 has many troubles with ip

	long factorDim = 0; // 103
	long sampleDim = 0; // 1579

	long** xyData = gd.xyDataFromFile(filename, factorDim, sampleDim, false);

	long sdimBits = (long)ceil(log2(sampleDim)); // 11
	long sampleDimPo2 = (1 << sdimBits); // 2048
	cout << "sampleDim: " << sampleDim << endl;
	cout << "sdimBits: " << sdimBits << endl;
	cout << "sampleDimPo2: " << sampleDimPo2 << endl;

	long fdimBits = (long)ceil(log2(factorDim));  // 7
	long factorDimPo2 = (1 << fdimBits); // 128
	cout << "factorDim: " << factorDim << endl;
	cout << "fdimBits: " << fdimBits << endl;
	cout << "factorDimPo2: " << factorDimPo2 << endl;

	//-----------------------------------------
	bool isEncrypted = true;
	bool isAllsample = true;
	bool is3approx = true; // 3 approx, 7 approx
	bool isFast = true; // always true
	long iter = fdimBits; // 7
	//-----------------------------------------

	long learnDim = isAllsample ? sampleDim : (1 << (sdimBits - 1)); // 1579
	long ldimBits = (long)ceil(log2(learnDim)); // 11
	long learnDimPo2 = (1 << ldimBits); // 2048
	cout << "learnDim: " << learnDim << endl;
	cout << "ldimBits: " << ldimBits << endl;
	cout << "learnDimPo2: " << learnDimPo2 << endl;

	long logl = 10;
	long logp = 32;
	long L = is3approx & isFast ? 4 * iter + 1 : !is3approx & !isFast ? 6 * iter + 1 : 5 * iter + 1;
	long logN = Params::suggestlogN(80, logl, logp, L);
	cout << "logl: " << logl << endl;
	cout << "logp: " << logp << endl;
	cout << "L: " << L << endl;
	cout << "logN: " << logN << endl;

	long xybatchBits = min(logN - 1 - ldimBits, fdimBits);
	long xyBatch = (1 << xybatchBits);
	cout << "xyBatchBits: " << xybatchBits << endl;
	cout << "xyBatch: " << xyBatch << endl;

	long slots =  learnDimPo2 * xyBatch;
	long cnum = factorDimPo2 / xyBatch;
	cout << "slots: " << slots << endl;
	cout << "cnum: " << cnum << endl;

	double* vData = new double[factorDim];
	double* wData = new double[factorDim];
	for (long i = 0; i < factorDim; ++i) {
		double tmp = 0.0; // averages
		wData[i] = tmp;
		vData[i] = tmp;
	}

	double* alpha = new double[iter + 2];
	alpha[0] = 0.01;
	for (long i = 1; i < iter + 2; ++i) {
		alpha[i] = (1. + sqrt(1. + 4.0 * alpha[i-1] * alpha[i-1])) / 2.0;
	}

	double* eta = new double[iter + 1];
	for (long i = 0; i < iter + 1; ++i) {
		eta[i] = (1. - alpha[i]) / alpha[i+1];
	}

	double* gamma = new double[iter];
	for (long i = 0; i < iter; ++i) {
		gamma[i] = 1. / learnDim / 3.;
//		gamma[i] = 1. / learnDim / (1. + k);
//		gamma[i] = 1. / learnDim / (2. + k);
	}

	if(!isEncrypted) {
		for (long k = 0; k < iter; ++k) {
			gd.stepNLGD(xyData, wData, vData, factorDim, learnDim, gamma[k], eta[k+1]);
			gd.check(xyData, wData, factorDim, sampleDim);
		}
	} else {
		timeutils.start("Scheme generating...");
		Params params(logN, logl, logp, L);
		SecKey secretKey(params);
		PubKey publicKey(params, secretKey);
		SchemeAux schemeaux(params);
		Scheme scheme(params, publicKey, schemeaux);
		SchemeAlgo algo(scheme);
		CipherGD cipherGD(scheme, algo, secretKey);
		timeutils.stop("Scheme generated");

		timeutils.start("Polynomial generating...");
		CZZ* pvals = new CZZ[slots];
		for (long j = 0; j < learnDim; ++j) {
			pvals[xyBatch * j] = CZZ(params.p);
		}
		CZZ* pdvals = scheme.groupidx(pvals, slots);
		Message msg = scheme.encode(pdvals, slots);
		timeutils.stop("Polynomial generated");

		timeutils.start("Encrypting xyData XYB...");
		Cipher* cxyData = cipherGD.encxyDataXYB(xyData, slots, factorDim, learnDim, learnDimPo2, xyBatch, cnum);
		timeutils.stop("xyData encrypted");

		timeutils.start("Encrypting wData and vData XYB...");
		Cipher* cwData = cipherGD.encwDataXYB(wData, slots, factorDim, learnDim, learnDimPo2, xyBatch, cnum);
		Cipher* cvData = new Cipher[cnum];
		for (long i = 0; i < cnum; ++i) {
			cvData[i] = cwData[i];
		}
		timeutils.stop("wData and vData encrypted");

		for (long k = 0; k < iter; ++k) {
			if(isFast) {
				if(is3approx) {
					timeutils.start("Encrypting NLGD step with 3 approx, 4 levels...");
					cipherGD.encStepNLGD3XYBfast4(cxyData, cwData, cvData, msg.mx, slots, learnDim, learnDimPo2, xybatchBits, xyBatch, cnum, gamma[k], eta[k+1], eta[k]);
					timeutils.stop("NLGD step with 3 approx, 4 levels finished");
				} else {
					timeutils.start("Encrypting NLGD step with 7 approx, 5 levels...");
					cipherGD.encStepNLGD7XYBfast5(cxyData, cwData, cvData, msg.mx, slots, learnDim, learnDimPo2, xybatchBits, xyBatch, cnum, gamma[k], eta[k+1], eta[k]);
					timeutils.stop("NLGD step with 7 approx, 5 levels finished");
				}
			} else {
				if(is3approx) {
					timeutils.start("Encrypting NLGD step with 3 approx, 5 levels...");
					cipherGD.encStepNLGD3XYB5(cxyData, cwData, cvData, msg.mx, slots, learnDim, learnDimPo2, xybatchBits, xyBatch, cnum, gamma[k], eta[k+1]);
					timeutils.stop("NLGD step with 3 approx, 5 levels finished");
				} else {
					timeutils.start("Encrypting NLGD step with 7 approx, 6 levels...");
					cipherGD.encStepNLGD7XYB6(cxyData, cwData, cvData, msg.mx, slots, learnDim, learnDimPo2, xybatchBits, xyBatch, cnum, gamma[k], eta[k+1]);
					timeutils.stop("NLGD step with 7 approx, 6 levels finished");
				}
			}

			timeutils.start("Decrypting wData");
			double* dw = cipherGD.decXYB(secretKey,cwData, factorDim, xyBatch, cnum);
			timeutils.stop("wData decrypted");

			gd.check(xyData, dw, factorDim, sampleDim);
		}
	}
	//-----------------------------------------
	cout << "!!! END TEST NLGD XYB !!!" << endl;
}
