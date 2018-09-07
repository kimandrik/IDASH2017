/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/

#include "TestGD.h"

#include "Ciphertext.h"
#include "NTL/ZZX.h"
#include "Scheme.h"
#include "SecretKey.h"
#include "TimeUtils.h"
#include <cmath>

#include "CipherGD.h"
#include "GD.h"

long TestGD::suggestLogN(long lambda, long logQ) {
	long NBnd = ceil(logQ * (lambda + 110) / 3.6);
	double logNBnd = log2((double)NBnd);
	return (long)ceil(logNBnd);
}

void TestGD::testEncNLGD(double** zDataTrain, double** zDataTest, long factorDim, long sampleDimTrain, long sampleDimTest,
		bool isYfirst, long numIter, long kdeg, double gammaUp, double gammaDown, bool isInitZero) {
	cout << "isYfirst = " << isYfirst << ", number of iterations = " << numIter << ", g_k = " << kdeg << endl;
	cout << "factors = " << factorDim << ", train samples = " << sampleDimTrain << ", test samples = " << sampleDimTest << endl;
	cout << "gammaUp = " << gammaUp << ", gammaDown = " << gammaDown << ", isInitZero = " << isInitZero << endl;

	long fdimBits = (long)ceil(log2(factorDim));
	long sdimBits = (long)ceil(log2(sampleDimTrain));

	long wBits = 30;
	long pBits = 20;
	long lBits = 5;
	long aBits = 3;
	long kBits = (long)ceil(log2(kdeg));

	long logQ = isInitZero ? (wBits + lBits) + numIter * ((kBits + 1) * wBits + 2 * pBits + aBits) :
			(sdimBits + wBits + lBits) + numIter * ((kBits + 1) * wBits + 2 * pBits + aBits);

	long logN = TestGD::suggestLogN(80, logQ);
	long bBits = min(logN - 1 - sdimBits, fdimBits);
	long batch = 1 << bBits;
	long sBits = sdimBits + bBits;
	long slots =  1 << sBits;
	long cnum = (long)ceil((double)factorDim / batch);

	cout << "batch = " << batch << ", slots = " << slots << ", cnum = " << cnum << endl;

	cout << "HEAAN PARAMETER logQ: " << logQ << endl;
	cout << "HEAAN PARAMETER logN: " << logN << endl;

	TimeUtils timeutils;
	timeutils.start("Scheme generating...");
	Ring ring(logN, logQ);
	SecretKey secretKey(ring);
	Scheme scheme(secretKey, ring);
	scheme.addLeftRotKeys(secretKey);
	scheme.addRightRotKeys(secretKey);
	timeutils.stop("Scheme generation");
	CipherGD cipherGD(scheme, secretKey);

	timeutils.start("Polynomial generating...");
	uint64_t* rpoly = cipherGD.generateAuxPoly(slots, batch, pBits);
	timeutils.stop("Polynomial generation");

	double* pwData = new double[factorDim];
	double* pvData = new double[factorDim];

	double* twData = new double[factorDim];
	double* tvData = new double[factorDim];

	double* cwData = new double[factorDim];
	double* cvData = new double[factorDim];

	Ciphertext* encZData = new Ciphertext[cnum];
	Ciphertext* encWData = new Ciphertext[cnum];
	Ciphertext* encVData = new Ciphertext[cnum];

	GD::normalizezData2(zDataTrain, zDataTest, factorDim, sampleDimTrain, sampleDimTest);

	timeutils.start("Encrypting zData...");
	cipherGD.encZData(encZData, zDataTrain, slots, factorDim, sampleDimTrain, batch, cnum, wBits, logQ);
	timeutils.stop("zData encryption");

	timeutils.start("Encrypting wData and vData...");
	if(isInitZero) {
		cipherGD.encWVDataZero(encWData, encVData, cnum, slots, wBits, logQ);
	} else {
		cipherGD.encWVDataAverage(encWData, encVData, encZData, cnum, sBits, bBits);
	}
	timeutils.stop("wData and vData encryption");

	if(isInitZero) {
		GD::initialWDataVDataZero(pwData, pvData, factorDim);
		GD::initialWDataVDataZero(twData, tvData, factorDim);
	} else {
		GD::initialWDataVDataAverage(pwData, pvData, zDataTrain, factorDim, sampleDimTrain);
		GD::initialWDataVDataAverage(twData, tvData, zDataTrain, factorDim, sampleDimTrain);
	}

	double alpha0, alpha1, eta, gamma;
	double enccor, encauc, truecor, trueauc;

	alpha0 = 0.01;
	alpha1 = (1. + sqrt(1. + 4.0 * alpha0 * alpha0)) / 2.0;

	for (long iter = 0; iter < numIter; ++iter) {
		cout << " !!! START " << iter + 1 << " ITERATION !!! " << endl;
		eta = (1 - alpha0) / alpha1;
		gamma = gammaDown > 0 ? gammaUp / gammaDown / sampleDimTrain : gammaUp / (iter - gammaDown) / sampleDimTrain;

		cout << "encWData.logq before: " << encWData[0].logq << endl;
		timeutils.start("Enc NLGD");
		cipherGD.encNLGDiteration(kdeg, encZData, encWData, encVData, rpoly, cnum, gamma, eta, sBits, bBits, wBits, pBits, aBits);
		timeutils.stop("Enc NLGD");
		cout << "encWData.logq after: " << encWData[0].logq << endl;

		cout << "----ENCRYPTED-----" << endl;
		cipherGD.decWData(cwData, encWData, factorDim, batch, cnum, wBits);
		GD::calculateAUC(zDataTest, cwData, factorDim, sampleDimTest, enccor, encauc);
		cout << "------------------" << endl;

		GD::trueNLGDiteration(zDataTrain, twData, tvData, factorDim, sampleDimTrain, gamma, eta);
//		cout << "-------TRUE-------" << endl;
//		GD::calculateAUC(zDataTest, twData, factorDim, sampleDimTest, truecor, trueauc);
//		cout << "------------------" << endl;

		alpha0 = alpha1;
		alpha1 = (1. + sqrt(1. + 4.0 * alpha0 * alpha0)) / 2.0;
		cout << " !!! STOP " << iter + 1 << " ITERATION !!! " << endl;
	}
	cout << "----ENCRYPTED-----" << endl;
	cipherGD.decWData(cwData, encWData, factorDim, batch, cnum, wBits);
	GD::calculateAUC(zDataTest, cwData, factorDim, sampleDimTest, enccor, encauc);
	cout << "------------------" << endl;

	cout << "-------TRUE-------" << endl;
	GD::calculateAUC(zDataTest, twData, factorDim, sampleDimTest, truecor, trueauc);
	cout << "------------------" << endl;

//	GD::calculateMSE(twData, cwData, factorDim);
//	GD::calculateNMSE(twData, cwData, factorDim);

}

void TestGD::testPlainNLGD(double** zDataTrain, double** zDataTest, long factorDim, long sampleDimTrain, long sampleDimTest,
		bool isYfirst, long numIter, long kdeg, double gammaUp, double gammaDown, bool isInitZero) {
	cout << "isYfirst = " << isYfirst << ", number of iterations = " << numIter << ", g_k = " << kdeg << endl;
	cout << "factors = " << factorDim << ", train samples = " << sampleDimTrain << ", test samples = " << sampleDimTest << endl;
	cout << "gammaUp = " << gammaUp << ", gammaDown = " << gammaDown << ", isInitZero = " << isInitZero << endl;

	TimeUtils timeutils;

	double* pwData = new double[factorDim];
	double* pvData = new double[factorDim];

	double* twData = new double[factorDim];
	double* tvData = new double[factorDim];

	GD::normalizezData2(zDataTrain, zDataTest, factorDim, sampleDimTrain, sampleDimTest);

	if(isInitZero) {
		GD::initialWDataVDataZero(pwData, pvData, factorDim);
		GD::initialWDataVDataZero(twData, tvData, factorDim);
	} else {
		GD::initialWDataVDataAverage(pwData, pvData, zDataTrain, factorDim, sampleDimTrain);
		GD::initialWDataVDataAverage(twData, tvData, zDataTrain, factorDim, sampleDimTrain);
	}

	double alpha0, alpha1, eta, gamma;
	double plaincor, plainauc, truecor, trueauc;

	alpha0 = 0.01;
	alpha1 = (1. + sqrt(1. + 4.0 * alpha0 * alpha0)) / 2.0;

	for (long iter = 0; iter < numIter; ++iter) {
		eta = (1 - alpha0) / alpha1;
		gamma = gammaDown > 0 ? gammaUp / gammaDown / sampleDimTrain : gammaUp / (iter - gammaDown) / sampleDimTrain;

		GD::plainNLGDiteration(kdeg, zDataTrain, pwData, pvData, factorDim, sampleDimTrain, gamma, eta);
		GD::trueNLGDiteration(zDataTrain, twData, tvData, factorDim, sampleDimTrain, gamma, eta);

		alpha0 = alpha1;
		alpha1 = (1. + sqrt(1. + 4.0 * alpha0 * alpha0)) / 2.0;
	}

	cout << "------PLAIN-------" << endl;
	GD::calculateAUC(zDataTest, pwData, factorDim, sampleDimTest, plaincor, plainauc);
	cout << "------------------" << endl;

	cout << "-------TRUE-------" << endl;
	GD::calculateAUC(zDataTest, twData, factorDim, sampleDimTest, truecor, trueauc);
	cout << "------------------" << endl;

//	GD::calculateMSE(twData, pwData, factorDim);
//	GD::calculateNMSE(twData, pwData, factorDim);

}

void TestGD::testEncNLGDFOLD(long fold, double** zData, long factorDim, long sampleDim,
			bool isYfirst, long numIter, long kdeg, double gammaUp, double gammaDown, bool isInitZero) {
	cout << "isYfirst = " << isYfirst << ", number of iterations = " << numIter << ", g_k = " << kdeg << endl;
	cout << "factors = " << factorDim << ", samples = " << sampleDim << ", fold = " << fold << endl;
	cout << "gammaUp = " << gammaUp << ", gammaDown = " << gammaDown << ", isInitZero = " << isInitZero << endl;

	long sampleDimTest = sampleDim / fold;
	long sampleDimTrain = sampleDim - sampleDimTest;

	long fdimBits = (long)ceil(log2(factorDim));
	long sdimBits = (long)ceil(log2(sampleDimTrain));

	long wBits = 30;
	long pBits = 20;
	long lBits = 5;
	long aBits = 3;
	long kBits = (long)ceil(log2(kdeg));

	long logQ = isInitZero ? (wBits + lBits) + numIter * ((kBits + 1) * wBits + 2 * pBits + aBits) :
			(sdimBits + wBits + lBits) + numIter * ((kBits + 1) * wBits + 2 * pBits + aBits);

	long logN = TestGD::suggestLogN(80, logQ);
	long bBits = min(logN - 1 - sdimBits, fdimBits);
	long batch = 1 << bBits;
	long sBits = sdimBits + bBits;
	long slots =  1 << sBits;
	long cnum = (long)ceil((double)factorDim / batch);

	cout << "batch = " << batch << ", slots = " << slots << ", cnum = " << cnum << endl;

	cout << "HEAAN PARAMETER logQ: " << logQ << endl;
	cout << "HEAAN PARAMETER logN: " << logN << endl;

	TimeUtils timeutils;
	timeutils.start("Scheme generating...");
	Ring ring(logN, logQ);
	SecretKey secretKey(ring);
	Scheme scheme(secretKey, ring);
	scheme.addLeftRotKeys(secretKey);
	scheme.addRightRotKeys(secretKey);
	timeutils.stop("Scheme generation");
	CipherGD cipherGD(scheme, secretKey);

	timeutils.start("Polynomial generating...");
	uint64_t* rpoly = cipherGD.generateAuxPoly(slots, batch, pBits);
	timeutils.stop("Polynomial generation");

	double* pwData = new double[factorDim];
	double* pvData = new double[factorDim];

	double* twData = new double[factorDim];
	double* tvData = new double[factorDim];

	double* cwData = new double[factorDim];
	double* cvData = new double[factorDim];

	Ciphertext* encZData = new Ciphertext[cnum];
	Ciphertext* encWData = new Ciphertext[cnum];
	Ciphertext* encVData = new Ciphertext[cnum];

	double **zDataTrain, **zDataTest;

	zDataTrain = new double*[sampleDimTrain];
	zDataTest = new double*[sampleDimTest];

	GD::normalizeZData(zData, factorDim, sampleDim);
	GD::shuffleZData(zData, factorDim, sampleDim);

	double enccor, encauc, truecor, trueauc;
	double averenccor = 0, averencauc = 0, avertruecor = 0, avertrueauc = 0;

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

		timeutils.start("Encrypting zData...");
		cipherGD.encZData(encZData, zDataTrain, slots, factorDim, sampleDimTrain, batch, cnum, wBits, logQ);
		timeutils.stop("zData encryption");

		timeutils.start("Encrypting wData and vData...");
		if(isInitZero) {
			cipherGD.encWVDataZero(encWData, encVData, cnum, slots, wBits, logQ);
		} else {
			cipherGD.encWVDataAverage(encWData, encVData, encZData, cnum, sBits, bBits);
		}
		timeutils.stop("wData and vData encryption");

		if(isInitZero) {
			GD::initialWDataVDataZero(pwData, pvData, factorDim);
			GD::initialWDataVDataZero(twData, tvData, factorDim);
		} else {
			GD::initialWDataVDataAverage(pwData, pvData, zDataTrain, factorDim, sampleDimTrain);
			GD::initialWDataVDataAverage(twData, tvData, zDataTrain, factorDim, sampleDimTrain);
		}

		//-----------------------------------------

		double alpha0, alpha1, eta, gamma;

		alpha0 = 0.01;
		alpha1 = (1. + sqrt(1. + 4.0 * alpha0 * alpha0)) / 2.0;

		for (long iter = 0; iter < numIter; ++iter) {
			cout << " !!! START " << iter + 1 << " ITERATION !!! " << endl;
			eta = (1 - alpha0) / alpha1;
			gamma = gammaDown > 0 ? gammaUp / gammaDown / sampleDimTrain : gammaUp / (iter - gammaDown) / sampleDimTrain;

			cout << "encWData.logq before: " << encWData[0].logq << endl;
			timeutils.start("Enc NLGD");
			cipherGD.encNLGDiteration(kdeg, encZData, encWData, encVData, rpoly, cnum, gamma, eta, sBits, bBits, wBits, pBits, aBits);
			timeutils.stop("Enc NLGD");
			cout << "encWData.logq after: " << encWData[0].logq << endl;

//			cout << "----ENCRYPTED-----" << endl;
//			cipherGD.decWData(cwData, encWData, factorDim, batch, cnum, wBits);
//			GD::calculateAUC(zDataTest, cwData, factorDim, sampleDimTest);
//			cout << "------------------" << endl;

			GD::plainNLGDiteration(kdeg, zDataTrain, pwData, pvData, factorDim, sampleDimTrain, gamma, eta);
			GD::trueNLGDiteration(zDataTrain, twData, tvData, factorDim, sampleDimTrain, gamma, eta);

//			cout << "-------TRUE-------" << endl;
//			GD::calculateAUC(zDataTest, twData, factorDim, sampleDimTest, correctness, auc);
//			cout << "------------------" << endl;

			alpha0 = alpha1;
			alpha1 = (1. + sqrt(1. + 4.0 * alpha0 * alpha0)) / 2.0;
			cout << " !!! STOP " << iter + 1 << " ITERATION !!! " << endl;
		}
		cout << "----ENCRYPTED-----" << endl;
		cipherGD.decWData(cwData, encWData, factorDim, batch, cnum, wBits);
		GD::calculateAUC(zDataTest, cwData, factorDim, sampleDimTest, enccor, encauc);
		cout << "------------------" << endl;

		cout << "-------TRUE-------" << endl;
		GD::calculateAUC(zDataTest, twData, factorDim, sampleDimTest, truecor, trueauc);
		cout << "------------------" << endl;

		averenccor += enccor;
		averencauc += encauc;
		avertruecor += truecor;
		avertrueauc += trueauc;

//		GD::calculateMSE(twData, cwData, factorDim);
//		GD::calculateNMSE(twData, cwData, factorDim);

		cout << " !!! STOP " << fnum + 1 << " FOLD !!! " << endl;
		cout << "------------------" << endl;
	}

	cout << "Average Encrypted correctness: " << averenccor << "%" << endl;
	cout << "Average Encrypted AUC: " << averencauc << endl;
	cout << "Average True correctness: " << avertruecor << "%" << endl;
	cout << "Average True AUC: " << avertrueauc << endl;
}

void TestGD::testPlainNLGDFOLD(long fold, double** zData, long factorDim, long sampleDim,
			bool isYfirst, long numIter, long kdeg, double gammaUp, double gammaDown, bool isInitZero) {
	cout << "isYfirst = " << isYfirst << ", number of iterations = " << numIter << ", g_k = " << kdeg << endl;
	cout << "factors = " << factorDim << ", samples = " << sampleDim << ", fold = " << fold << endl;
	cout << "gammaUp = " << gammaUp << ", gammaDown = " << gammaDown << ", isInitZero = " << isInitZero << endl;

	long sampleDimTest = sampleDim / fold;
	long sampleDimTrain = sampleDim - sampleDimTest;

	TimeUtils timeutils;

	double* pwData = new double[factorDim];
	double* pvData = new double[factorDim];

	double* twData = new double[factorDim];
	double* tvData = new double[factorDim];

	double **zDataTrain, **zDataTest;

	zDataTrain = new double*[sampleDimTrain];
	zDataTest = new double*[sampleDimTest];

	GD::normalizeZData(zData, factorDim, sampleDim);
	GD::shuffleZData(zData, factorDim, sampleDim);

	double plaincor, plainauc, truecor, trueauc;
	double averplaincor = 0, averplainauc = 0, avertruecor = 0, avertrueauc = 0;

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

		double alpha0, alpha1, eta, gamma;

		alpha0 = 0.01;
		alpha1 = (1. + sqrt(1. + 4.0 * alpha0 * alpha0)) / 2.0;

		for (long iter = 0; iter < numIter; ++iter) {
			eta = (1 - alpha0) / alpha1;
			cout << eta << endl;
			gamma = gammaDown > 0 ? gammaUp / gammaDown / sampleDimTrain : gammaUp / (iter - gammaDown) / sampleDimTrain;
			GD::plainNLGDiteration(kdeg, zDataTrain, pwData, pvData, factorDim, sampleDimTrain, gamma, eta);
			GD::trueNLGDiteration(zDataTrain, twData, tvData, factorDim, sampleDimTrain, gamma, eta);

			alpha0 = alpha1;
			alpha1 = (1. + sqrt(1. + 4.0 * alpha0 * alpha0)) / 2.0;
		}
		cout << "------PLAIN-------" << endl;
		GD::calculateAUC(zDataTest, pwData, factorDim, sampleDimTest, plaincor, plainauc);
		cout << "------------------" << endl;

		cout << "-------TRUE-------" << endl;
		GD::calculateAUC(zDataTest, twData, factorDim, sampleDimTest, truecor, trueauc);
		cout << "------------------" << endl;

		averplaincor += plaincor;
		averplainauc += plainauc;
		avertruecor += truecor;
		avertrueauc += trueauc;

		cout << " !!! STOP " << fnum + 1 << " FOLD !!! " << endl;
		cout << "------------------" << endl;
	}

	cout << "Average Plain correctness: " << averplaincor << "%" << endl;
	cout << "Average Plain AUC: " << averplainauc << endl;
	cout << "Average True correctness: " << avertruecor << "%" << endl;
	cout << "Average True AUC: " << avertrueauc << endl;
}

