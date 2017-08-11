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

void TestAK::testNLGD(string filename, long iter, double gammaCnst, double gammaUpCnst, double learnPortion, bool is3approx, bool isEncrypted, bool isYfirst, long xyBits, long wBits, long pBits, long lBits) {
	cout << "!!! START TEST NLGD !!!" << endl;
	//-----------------------------------------
	TimeUtils timeutils;
	SetNumThreads(8);
	//-----------------------------------------
	GD gd;

	cout << "iter: " << iter << endl;
	cout << "isEncrypted: " << isEncrypted << endl;
	cout << "is3approx: " << is3approx << endl;
	cout << "gammaCnst: " << gammaCnst << endl;
	cout << "gammaUpCnst: " << gammaUpCnst << endl;

	long factorDim = 0;
	long sampleDim = 0;

	long** xyData = gd.xyDataFromFile(filename, factorDim, sampleDim, isYfirst);

	long sdimBits = (long)ceil(log2(sampleDim));
	long sampleDimPo2 = (1 << sdimBits);
	cout << "sampleDim: " << sampleDim << endl;
	cout << "sdimBits: " << sdimBits << endl;
	cout << "sampleDimPo2: " << sampleDimPo2 << endl;

	long fdimBits = (long)ceil(log2(factorDim));
	long factorDimPo2 = (1 << fdimBits);
	cout << "factorDim: " << factorDim << endl;
	cout << "fdimBits: " << fdimBits << endl;
	cout << "factorDimPo2: " << factorDimPo2 << endl;

	long learnDim = (long)((double)sampleDim * learnPortion);
	long ldimBits = (long)ceil(log2(learnDim));
	long learnDimPo2 = (1 << ldimBits);
	cout << "learnDim: " << learnDim << endl;
	cout << "ldimBits: " << ldimBits << endl;
	cout << "learnDimPo2: " << learnDimPo2 << endl;

	long** xyDataLearn = gd.pickxyDataLearn(xyData, learnDim, sampleDim, factorDim);

	long logq = is3approx ? iter * (2 * wBits + xyBits + pBits) + lBits + ldimBits + xyBits - wBits
			: iter * (3 * wBits + xyBits + pBits) + lBits + ldimBits + xyBits - wBits;

	long logN = Params::suggestlogN(80, logq);
	cout << "logq: " << logq << endl;
	cout << "logN: " << logN << endl;

	long xybatchBits = min(logN - 1 - ldimBits, fdimBits);
	long xyBatch = 1 << xybatchBits;
	cout << "xyBatchBits: " << xybatchBits << endl;
	cout << "xyBatch: " << xyBatch << endl;

	long slotBits = ldimBits + xybatchBits;
	long slots =  1 << slotBits;
	long cnum = factorDimPo2 / xyBatch;
	cout << "slots: " << slots << endl;
	cout << "cnum: " << cnum << endl;

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
	if(gammaCnst > 0) {
		for (long i = 0; i < iter; ++i) {
			gamma[i] = gammaUpCnst / learnDim / gammaCnst;
		}
	} else {
		for (long i = 0; i < iter; ++i) {
			gamma[i] = gammaUpCnst / learnDim / (i - gammaCnst);
		}
	}

	if(!isEncrypted) {
		for (long k = 0; k < iter; ++k) {
			gd.stepNLGD(xyDataLearn, wData, vData, factorDim, learnDim, gamma[k], eta[k+1]);
			gd.check(xyData, wData, factorDim, sampleDim);
		}
	} else {
		timeutils.start("Scheme generating...");
		Params params(logN, logq);
		SecKey secretKey(params);
		PubKey publicKey(params, secretKey);
		SchemeAux schemeaux(params, wBits);
		Scheme scheme(params, publicKey, schemeaux);
		SchemeAlgo algo(scheme);
		CipherGD cipherGD(scheme, algo, secretKey);
		timeutils.stop("Scheme generated");

		timeutils.start("Polynomial generating...");
		ZZ p = power2_ZZ(pBits);
		CZZ* pvals = new CZZ[slots];
		for (long j = 0; j < learnDim; ++j) {
			pvals[xyBatch * j] = CZZ(p);
		}
		CZZ* pdvals = scheme.groupidx(pvals, slots);
		Message msg = scheme.encode(pdvals, slots);
		timeutils.stop("Polynomial generated");

		timeutils.start("Encrypting xyData...");
		Cipher* cxyData = cipherGD.encxyData(xyDataLearn, slots, factorDim, learnDim, learnDimPo2, xyBatch, cnum, xyBits);
		timeutils.stop("xyData encrypted");

		timeutils.start("Encrypting wData and vData...");
		Cipher* cwData = cipherGD.encwData(cxyData, slotBits, ldimBits, xybatchBits, cnum, xyBits, wBits);
		Cipher* cvData = new Cipher[cnum];
		for (long i = 0; i < cnum; ++i) {
			cvData[i] = cwData[i];
		}
		timeutils.stop("wData and vData encrypted");

		for (long k = 0; k < iter; ++k) {
			if(is3approx) {
				timeutils.start("Encrypting NLGD step with degree 3 approx...");
				cipherGD.encStepNLGD3(cxyData, cwData, cvData, msg.mx, slots, learnDim, learnDimPo2, xybatchBits, xyBatch, cnum, gamma[k], eta[k+1], eta[k], xyBits, wBits, pBits);
				timeutils.stop("NLGD step with degree 3 approx, finished");
			} else {
				timeutils.start("Encrypting NLGD step with degree 5 approx...");
				cipherGD.encStepNLGD5(cxyData, cwData, cvData, msg.mx, slots, learnDim, learnDimPo2, xybatchBits, xyBatch, cnum, gamma[k], eta[k+1], eta[k], xyBits, wBits, pBits);
				timeutils.stop("NLGD step with degree 5 approx finished");
			}

			timeutils.start("Decrypting wData");
			double* dw = cipherGD.decwData(secretKey,cwData, factorDim, xyBatch, cnum, wBits);
			timeutils.stop("wData decrypted");

			gd.check(xyData, dw, factorDim, sampleDim);
		}
	}
	//-----------------------------------------
	cout << "!!! END TEST NLGD !!!" << endl;
}
