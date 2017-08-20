#include <NTL/BasicThreadPool.h>
#include <Cipher.h>
#include <CZZ.h>
#include <Message.h>
#include <Params.h>
#include <PubKey.h>
#include <Scheme.h>
#include <SchemeAlgo.h>
#include <SchemeAux.h>
#include <SecKey.h>
#include <TimeUtils.h>
#include <algorithm>
#include <cmath>
#include <iostream>

#include "CipherGD.h"
#include "GD.h"

using namespace std;

void run(string filename, long iter, double gammaDownCnst, double gammaUpCnst, double learnPortion, bool is3approx, bool isEncrypted, bool isYfirst) {
	cout << "!!! START NLGD !!!" << endl;
	//-----------------------------------------
	TimeUtils timeutils;
	SetNumThreads(8);
	//-----------------------------------------

	cout << "iter: " << iter << endl;
	cout << "isEncrypted: " << isEncrypted << endl;
	cout << "is3approx: " << is3approx << endl;
	cout << "gammaDownCnst: " << gammaDownCnst << endl;
	cout << "gammaUpCnst: " << gammaUpCnst << endl;

	long factorDim = 0;
	long sampleDim = 0;

	long** xyData = GD::xyDataFromFile(filename, factorDim, sampleDim, isYfirst);
	cout << "sampleDim: " << sampleDim << endl;
	cout << "factorDim: " << factorDim << endl;

	long fdimBits = (long)ceil(log2(factorDim));
	cout << "fdimBits: " << fdimBits << endl;

	long learnDim = (long)((double)sampleDim * learnPortion);
	cout << "learnDim: " << learnDim << endl;

	long ldimBits = (long)ceil(log2(learnDim));
	cout << "ldimBits: " << ldimBits << endl;

	long** xyDataLearn = GD::RandomxyDataLearn(xyData, learnDim, sampleDim, factorDim);

	long wBits = max(fdimBits + ldimBits + 18, 28);
	long xyBits = wBits + 2;
	long pBits = 16;
	long lBits = 5;
	cout << "xyBits: " << xyBits << endl;
	cout << "wBits: " << wBits << endl;
	cout << "pBits: " << pBits << endl;
	cout << "lBits: " << lBits << endl;

	long logq = is3approx ? (ldimBits + xyBits) + iter * (2 * wBits + xyBits + pBits) + lBits
			: (ldimBits + xyBits) + iter * (3 * wBits + xyBits + pBits) + lBits;

	long logN = Params::suggestlogN(80, logq);
	cout << "logq: " << logq << endl;
	cout << "logN: " << logN << endl;

	long xybatchBits = min(logN - 1 - ldimBits, fdimBits);
	long xyBatch = 1 << xybatchBits;
	cout << "xyBatchBits: " << xybatchBits << endl;

	long slotBits = ldimBits + xybatchBits;
	long slots =  1 << slotBits;
	long cnum = (long)ceil((double)factorDim / xyBatch);
	cout << "slots: " << slots << endl;
	cout << "cnum: " << cnum << endl;

	double alpha0, alpha1, alpha2;
	double eta0, eta1;
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

			eta1 = (1 - alpha1) / alpha2;
			gamma = gammaDownCnst > 0 ? gammaUpCnst / gammaDownCnst / learnDim : gammaUpCnst / (k - gammaDownCnst) / learnDim;

			GD::stepNLGD(xyDataLearn, wData, vData, factorDim, learnDim, gamma, eta1);
			GD::check(xyData, wData, factorDim, sampleDim);

			alpha0 = alpha1;
			alpha1 = alpha2;
			alpha2 = (1. + sqrt(1. + 4.0 * alpha1 * alpha1)) / 2.0;

		}
	} else {
		timeutils.start("Scheme generating...");
		Params params(logN, logq);
		SecKey secretKey(params);
		PubKey publicKey(params, secretKey);
		SchemeAux schemeaux(logN);
		Scheme scheme(params, publicKey, schemeaux);
		CipherGD cipherGD(scheme, secretKey);
		timeutils.stop("Scheme generated");

		timeutils.start("Polynomial generating...");
		ZZ p = power2_ZZ(pBits);
		CZZ* pvals = new CZZ[slots];
		for (long j = 0; j < learnDim; ++j) {
			pvals[xyBatch * j] = CZZ(p);
		}
		CZZ* pdvals = scheme.groupidx(pvals, slots);
		Message msg = scheme.encode(pdvals, slots);

		delete[] pvals;
		delete[] pdvals;

		timeutils.stop("Polynomial generated");

		Cipher* cxyData = new Cipher[cnum];
		timeutils.start("Encrypting xyData...");
		cipherGD.encxyData(cxyData, xyDataLearn, slots, factorDim, learnDim, xyBatch, cnum, xyBits);
		timeutils.stop("xyData encrypted");

		Cipher* cwData = new Cipher[cnum];
		Cipher* cvData = new Cipher[cnum];
		timeutils.start("Encrypting wData and vData...");
		cipherGD.encwData(cwData, cxyData, slotBits, ldimBits, xybatchBits, cnum, xyBits, wBits);
		for (long i = 0; i < cnum; ++i) {
			cvData[i] = cwData[i];
		}
		timeutils.stop("wData and vData encrypted");

		double* dwData = new double[factorDim];
		for (long k = 0; k < iter; ++k) {

			eta0 = (1 - alpha0) / alpha1;
			eta1 = (1 - alpha1) / alpha2;
			gamma = gammaDownCnst > 0 ? gammaUpCnst / gammaDownCnst / learnDim : gammaUpCnst / (k - gammaDownCnst) / learnDim;

			if(is3approx) {
				timeutils.start("Encrypting NLGD step with degree 3 approx...");
				cipherGD.encStepNLGD3(cxyData, cwData, cvData, msg.mx, slotBits, xybatchBits, cnum, gamma, eta1, eta0, xyBits, wBits, pBits);
				timeutils.stop("NLGD step with degree 3 approx, finished");
			} else {
				timeutils.start("Encrypting NLGD step with degree 5 approx...");
				cipherGD.encStepNLGD5(cxyData, cwData, cvData, msg.mx, slotBits, xybatchBits, cnum, gamma, eta1, eta0, xyBits, wBits, pBits);
				timeutils.stop("NLGD step with degree 5 approx finished");
			}

			alpha0 = alpha1;
			alpha1 = alpha2;
			alpha2 = (1. + sqrt(1. + 4.0 * alpha1 * alpha1)) / 2.0;

			timeutils.start("Decrypting wData");
			cipherGD.decwData(dwData, cwData, factorDim, xyBatch, cnum, wBits);
			timeutils.stop("wData decrypted");

			GD::check(xyData, dwData, factorDim, sampleDim);
		}
	}
	//-----------------------------------------
	cout << "!!! END NLGD !!!" << endl;
}

int main() {

//	string filename = "data/data5x500.txt";    bool isYfirst = false; //  421/500
//	string filename = "data/data9x1253.txt";   bool isYfirst = false; //  1147/1253
//	string filename = "data/data15x1500.txt";  bool isYfirst = false; //  1277/1500
//	string filename = "data/data16x101.txt";   bool isYfirst = false; //  101/101
//	string filename = "data/data27x148.txt";   bool isYfirst = false; //  132/148
//	string filename = "data/data43x3247.txt";  bool isYfirst = false; //  3182/3247
//	string filename = "data/data45x296.txt";   bool isYfirst = false; //  257/296
//	string filename = "data/data51x653.txt";   bool isYfirst = false; //  590/653
//	string filename = "data/data67x216.txt";   bool isYfirst = false; //  216/216
	string filename = "data/data103x1579.txt"; bool isYfirst = true;  //  1086/1579

	long iter = 7;

	/*
	 * gammaDownCnst > 0 : gamma = gammaUpCnst / gammaDownCnst / learnDim - constant gamma
	 * gammaDownCnst < 0 : gamma = gammaUpCnst / (i + |gammaDownCnst|) / learnDim - decreasing gamma
	 */
	double gammaDownCnst = -2.;
	double gammaUpCnst = 2.;

	/*
	 * portion used in learning (randomly chosen from sample set)
	 */
	double learnPortion = 0.9;

	/*
	 * false : use 5 degree polynomial approximation of sigmoid function
	 * true  : use 3 degree polynomial approximation of sigmoid function
	 */
	bool is3approx = false;

	/*
	 * false : unencrypted Nesterov Logistic Gradient Descent
	 * true  : encrypted Nesterov Logistic Gradient Descent
	 */
	bool isEncrypted = false;

	run(filename, iter, gammaDownCnst, gammaUpCnst, learnPortion, is3approx, isEncrypted, isYfirst);

	return 0;
}
