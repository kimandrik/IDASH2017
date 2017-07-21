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
	long sampledimbits = (long)ceil(log2(sampledim));
	long po2sampledim = (1 << sampledimbits); // 1579 -> 2048
	long wnum = slots / po2sampledim; // N / 2 / 2048
	long dimbits = (long)ceil(log2(dim));
	long po2dim = (1 << dimbits); //103 -> 128

	cout << "dimension: " << dim << endl;
	cout << "power of 2 dimension: " << po2dim << endl;

	cout << "sample dimension: " << sampledim << endl;
	cout << "power of 2 sample dimension: " << po2sampledim << endl;

	cout << "slots: " << slots << endl;
	cout << "wnum: " << wnum << endl;

	long iter = 200;
	double** wdata = sgd.wdatagen(wnum, dim);
	double* gamma = sgd.gammagen(iter);

	double lambda = 2.0;
	timeutils.start("sgd");
	for (long k = 0; k < iter; ++k) {
		NTL_EXEC_RANGE(wnum, first, last);
		for (long l = first; l < last; ++l) {
//			sgd.steplogregress(wdata[l], zdata, gamma[k], lambda, dim, sampledim);
			sgd.stepsimpleregress(wdata[l], zdata, gamma[k], lambda, dim, sampledim);
		}
		NTL_EXEC_RANGE_END;
	}
	timeutils.stop("sgd");

	double* w = sgd.wout(wdata, wnum, dim);

	sgd.check(w, zdata, dim, sampledim);

	//-----------------------------------------
	Params params(logN, logl, logp, L);
	SecKey secretKey(params);
	PubKey publicKey(params, secretKey);
	SchemeAux schemeaux(params);
	Scheme scheme(params, publicKey, schemeaux);
	SchemeAlgo algo(scheme);
	CipherSGD csgd(scheme, algo);
	//-----------------------------------------

	timeutils.start("Enc zdata");
	Cipher* czdata = csgd.enczdata(zdata, slots, wnum, dim, sampledim, params.p);
	timeutils.stop("Enc zdata");

	timeutils.start("Enc wdata");
	Cipher* cwdata = csgd.encwdata(wdata, slots, wnum, dim, sampledim, params.logp);
	timeutils.stop("Enc wdata");

	ZZ* pgamma = csgd.pgammagen(gamma, iter, params.logp);

	iter = 100;
	//-----------------------------------------
	for (long k = 0; k < iter; ++k) {
		cout << k << endl;
		timeutils.start("Enc sgd step");
		csgd.encSteplogregress(czdata, cwdata, pgamma[k], lambda, slots, wnum, dim);
		timeutils.stop("Enc sgd step");
	}

	timeutils.start("Enc w out");
	Cipher* cw = csgd.encwout(cwdata, wnum, dim);
	timeutils.start("Enc w out");

	timeutils.start("Dec w");
	double* dw = csgd.decw(secretKey, cw, dim);
	timeutils.start("Dec w");

	sgd.check(dw, zdata, dim, sampledim);
	//-----------------------------------------
	cout << "!!! END TEST SGD !!!" << endl;
}
