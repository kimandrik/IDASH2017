#include "TestSGD.h"

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

void TestSGD::testSGD(long logN, long logl, long logp, long L) {
	cout << "!!! START TEST SGD !!!" << endl;
	//-----------------------------------------
	TimeUtils timeutils;
	SetNumThreads(8);
	//-----------------------------------------
	SGD sgd;
	string filename = "data.txt";

	long dim = 0;
	long sampledim = 0;

	long** zdata = sgd.zdataFromFile(filename, dim, sampledim); //dim = 103, sampledim = 1579

	long slots =  (1 << (logN-1)); // N /2
	long sampledimbits = (long)ceil(log2(sampledim)); // log(1579) = 11
	long po2sampledim = (1 << sampledimbits); // 1579 -> 2048
//	long learndim = (1 << (sampledimbits - 1)); // 1024
//	long wnum = slots / learndim;
	long learndim = sampledim;
	long wnum = slots / po2sampledim; // N / 2 / 2048
	long dimbits = (long)ceil(log2(dim)); // log(103) = 7
	long po2dim = (1 << dimbits); //103 -> 128

	cout << "dimension: " << dim << endl;
	cout << "power of 2 dimension: " << po2dim << endl;

	cout << "sample dimension: " << sampledim << endl;
	cout << "power of 2 sample dimension: " << po2sampledim << endl;

	cout << "learn dimension: " << learndim << endl;

	cout << "slots: " << slots << endl;
	cout << "wnum: " << wnum << endl;

	double** vdata = new double*[wnum];
	for (long l = 0; l < wnum; ++l) {
		vdata[l] = new double[dim];
	}

	double** wdata = new double*[wnum];
	for (long l = 0; l < wnum; ++l) {
		wdata[l] = new double[dim];
		for (long i = 0; i < dim; ++i) {
			wdata[l][i] = (1.0 - 2.0 * (double)rand() / RAND_MAX) / 2.0; // change to good initial w choice
		}
	}

	long iter = 20;
	long additer = 3;
	long totaliter = iter + additer;

	double* gamma = new double[iter];
	for (long k = 0; k < iter; ++k) {
//		gamma[k] = 0.1 / (k + 1);
//		gamma[k] = 1.0 / (k + 1);
		gamma[k] = (double)learndim / (k + 1);
//		gamma[k] = 10.;
	}

	double eta = 0.6;
	double beta = 0.8;
//	double lambda = 2.0;
	double lambda = 2.0 / learndim;

	timeutils.start("sgd");
	for (long k = 0; k < iter; ++k) {
		NTL_EXEC_RANGE(wnum, first, last);
		for (long l = first; l < last; ++l) {
//			sgd.stepQuadraticRegress(wdata[l], zdata, gamma[k], lambda, dim, learndim);
//			sgd.stepLogRegress(wdata[l], zdata, gamma[k], lambda, dim, learndim);
//			sgd.stepMomentumLogRegress(wdata[l], vdata[l], zdata, gamma[k], lambda, dim, learndim, eta);
			sgd.stepNesterovLogRegress(wdata[l], vdata[l], zdata, gamma[k], lambda, dim, learndim, beta, eta);
		}
		NTL_EXEC_RANGE_END;
	}
	timeutils.stop("sgd");

	double* w = sgd.waverage(wdata, wnum, dim);

	sgd.check(w, zdata, dim, sampledim);

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
	Cipher* czdata = csgd.enczdata(zdata, slots, wnum, dim, learndim, params.p);
	timeutils.stop("Enc zdata");

	timeutils.start("Enc wdata");
	Cipher* cwdata = csgd.encwdata(wdata, slots, wnum, dim, learndim, params.logp);
	timeutils.stop("Enc wdata");

	ZZ* pgamma = csgd.pgammagen(gamma, totaliter, params.logp);

	//-----------------------------------------
	for (long k = iter; k < totaliter; ++k) {
		cout << k << endl;
		timeutils.start("Enc sgd step");
		csgd.encStepQuadraticRegress(czdata, cwdata, pgamma[k], lambda, slots, wnum, dim, learndim);

		timeutils.stop("Enc sgd step");
	}

	timeutils.start("Enc w out");
	Cipher* cw = csgd.encwout(cwdata, wnum, dim);
	timeutils.stop("Enc w out");

	timeutils.start("Dec w");
	double* dw = csgd.decw(secretKey, cw, dim);
	timeutils.stop("Dec w");

	sgd.check(dw, zdata, dim, sampledim);
	//-----------------------------------------
	cout << "!!! END TEST SGD !!!" << endl;
}
