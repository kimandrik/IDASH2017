#include <Cipher.h>
#include <CZZ.h>
#include <Message.h>
#include <NTL/BasicThreadPool.h>
#include <NTL/tools.h>
#include <NTL/ZZ.h>
#include <Params.h>
#include <PubKey.h>
#include <Scheme.h>
#include <SchemeAux.h>
#include <TestScheme.h>
#include <SecKey.h>
#include <TimeUtils.h>
#include <cmath>
#include <iostream>

#include "CipherGD.h"
#include "GD.h"

using namespace std;

/*
 * run: ./HEML filename(string) isYfirst(bool) iter(long) learnPortion(double) approx(long) isEncrypted(bool)
 *
 * example: ./HEML "../data/data103x1579.txt" 1 7 1 7 1
 *
 * parameters:
 * filename - path to file
 * isYfirst - {0,1} y parameter first OR last
 * iter - number of iterations
 * learnPortion - portion of data used for learning (randomly chosen from sample set)
 * approx - {3,5,7} polynomial degree approximation used
 * isEncrypted - {0,1} encrypted OR unencrypted learning
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
 *
 * FYI: approx 3 suggested iter: 4, 9, 18, 36, ...
 * FYI: approx 5 suggested iter: 3, 7, 14, 28, ...
 * FYI: approx 7 suggested iter: 3, 7, 14, 28, ...
 */

int main(int argc, char **argv) {

	string filename(argv[1]);
	bool isYfirst = atoi(argv[2]);
	long iter = atol(argv[3]);
	double learnPortion = atof(argv[4]);
	long approx = atol(argv[5]);
	bool isEncrypted = atoi(argv[6]);

	TimeUtils timeutils;
	SetNumThreads(8);

	cout << "iter: " << iter << endl;
	cout << "isEncrypted: " << isEncrypted << endl;
	cout << "approx: " << approx << endl;

	/*
	 * gammaDownCnst > 0 : gamma = gammaUpCnst / gammaDownCnst / learnDim -> constant gamma
	 * gammaDownCnst < 0 : gamma = gammaUpCnst / (i + |gammaDownCnst|) / learnDim -> decreasing gamma
	 */
	double gammaUpCnst = 1;
	double gammaDownCnst = -1;

	cout << "gammaUpCnst: " << gammaUpCnst << endl;
	cout << "gammaDownCnst: " << gammaDownCnst << endl;

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

	long wBits = 37;
	long xyBits = 37;

	cout << "xyBits: " << xyBits << endl;
	cout << "wBits: " << wBits << endl;

	long lBits = 5;
	long pBits = 18;
	long aBits = 2;
	cout << "lBits: " << lBits << endl;
	cout << "pBits: " << pBits << endl;
	cout << "aBits: " << aBits << endl;

	long logq = (approx == 3) ? (ldimBits + xyBits) + iter * (2 * wBits + xyBits + pBits + aBits) + lBits
			: (ldimBits + xyBits) + iter * (3 * wBits + xyBits + pBits + aBits) + lBits;

	long logN = Params::suggestlogN(80, logq);
	cout << "logq: " << logq << endl;
	cout << "logN: " << logN << endl;

	long bBits = min(logN - 1 - ldimBits, fdimBits);
	long batch = 1 << bBits;
	cout << "bBits: " << bBits << endl;

	long sBits = ldimBits + bBits;
	long slots =  1 << sBits;
	long cnum = (long)ceil((double)factorDim / batch);
	cout << "slots: " << slots << endl;
	cout << "cnum: " << cnum << endl;

	double alpha0, alpha1, alpha2;
	double etaprev, eta;
	double gamma;

	alpha0 = 0.01;

	alpha1 = (1. + sqrt(1. + 4.0 * alpha0 * alpha0)) / 2.0;
	alpha2 = (1. + sqrt(1. + 4.0 * alpha1 * alpha1)) / 2.0;

	cout << "!!! START NLGD !!!" << endl;

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
		timeutils.stop("Scheme generation");

		timeutils.start("Polynomial generating...");
		ZZ p = power2_ZZ(pBits);
		CZZ* pvals = new CZZ[slots];
		for (long j = 0; j < slots; j += batch) {
			pvals[j] = CZZ(p);
		}
		CZZ* pdvals = scheme.groupidx(pvals, slots);
		delete[] pvals;

		Message msg = scheme.encode(pdvals, slots);
		delete[] pdvals;

		timeutils.stop("Polynomial generation");

		Cipher* cxyData = new Cipher[cnum];
		timeutils.start("Encrypting xyData...");
		cipherGD.encxyData(cxyData, xyDataLearn, slots, factorDim, learnDim, batch, cnum, xyBits);
		timeutils.stop("xyData encryption");

		Cipher* cwData = new Cipher[cnum];
		Cipher* cvData = new Cipher[cnum];
		timeutils.start("Encrypting wData and vData...");
		cipherGD.encwData(cwData, cxyData, cnum, sBits, ldimBits, bBits, xyBits, wBits);
		for (long i = 0; i < cnum; ++i) {
			cvData[i] = cwData[i];
		}
		timeutils.stop("wData and vData encryption");

		double* dwData = new double[factorDim];
		for (long k = 0; k < iter; ++k) {

			etaprev = (1 - alpha0) / alpha1;
			eta = (1 - alpha1) / alpha2;
			gamma = gammaDownCnst > 0 ? gammaUpCnst / gammaDownCnst / learnDim : gammaUpCnst / (k - gammaDownCnst) / learnDim;

			if(approx == 3) {
				timeutils.start("Encrypting NLGD step with degree 3 approx...");
				cipherGD.encStepNLGD3(cxyData, cwData, cvData, msg.mx, cnum, gamma, eta, etaprev, sBits, bBits, xyBits, wBits, pBits, aBits);
				timeutils.stop("NLGD step with degree 3 approx");
			} else if(approx == 5) {
				timeutils.start("Encrypting NLGD step with degree 5 approx...");
				cipherGD.encStepNLGD5(cxyData, cwData, cvData, msg.mx, cnum, gamma, eta, etaprev, sBits, bBits, xyBits, wBits, pBits, aBits);
				timeutils.stop("NLGD step with degree 5 approx");
			} else {
				timeutils.start("Encrypting NLGD step with degree 7 approx...");
				cipherGD.encStepNLGD7(cxyData, cwData, cvData, msg.mx, cnum, gamma, eta, etaprev, sBits, bBits, xyBits, wBits, pBits, aBits);
				timeutils.stop("NLGD step with degree 7 approx");
			}

			alpha0 = alpha1;
			alpha1 = alpha2;
			alpha2 = (1. + sqrt(1. + 4.0 * alpha1 * alpha1)) / 2.0;

			timeutils.start("check");
			cipherGD.decwData(dwData, cwData, factorDim, batch, cnum, wBits);
			GD::check(xyData, dwData, factorDim, sampleDim);
			timeutils.stop("check");
		}
	}
	cout << "!!! END NLGD !!!" << endl;

	return 0;
}
