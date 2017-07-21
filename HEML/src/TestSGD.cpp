#include "TestSGD.h"

#include <NTL/ZZ.h>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>

#include <NTL/BasicThreadPool.h>

#include "SGD.h"
#include <CZZ.h>
#include <Cipher.h>
#include <Params.h>
#include <PubKey.h>
#include <Scheme.h>
#include <SchemeAlgo.h>
#include <SecKey.h>
#include <StringUtils.h>
#include <TimeUtils.h>
#include <EvaluatorUtils.h>

using namespace NTL;

//-----------------------------------------

void TestSGD::testSGD(long logN, long logl, long logp, long L) {
	cout << "!!! START TEST SGD !!!" << endl;
	//-----------------------------------------
	TimeUtils timeutils;
	Params params(logN, logl, logp, L);
	SecKey secretKey(params);
	PubKey publicKey(params, secretKey);
	SchemeAux schemeaux(params);
	Scheme scheme(params, publicKey, schemeaux);
	SchemeAlgo algo(scheme);
	SGD sgd(scheme, algo);
	//-----------------------------------------
	SetNumThreads(8);
	//-----------------------------------------
	string filename = "data.txt";

	long dim = 0;
	long sampledim = 0;

	long** zdata = sgd.zdataFromFile(filename, dim, sampledim); //dim = 103, sampledim = 1579

	long slots = params.N / 2; // N /2
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

	long iter = 1000;
	double** wdata = sgd.wdatagen(wnum, dim);
	double* gamma = sgd.gammagen(iter);

	double lambda = 2.0;
	for (long k = 0; k < iter; ++k) {
		sgd.step(wdata, zdata, gamma[k], lambda, wnum, dim, sampledim);
	}

	double* w = wgen(wdata, wnum, dim);

	sgd.check(w, zdata, dim, sampledim);

	//-----------------------------------------

	timeutils.start("Enc zdata");
	Cipher* czdata = sgd.enczdata(zdata, slots, wnum, dim, sampledim, params.p);
	timeutils.stop("Enc zdata");

	timeutils.start("Enc wdata");
	Cipher* cwdata = sgd.encwdata(wdata, slots, wnum, dim, sampledim, params.logp);
	timeutils.stop("Enc wdata");

	ZZ* pgamma = sgd.pgammagen(gamma, iter, params.logp)

	iter = 100;
	//-----------------------------------------
	for (long k = 0; k < iter; ++k) {
		timeutils.start("Enc sgd step " + k);
		sgd.encStep(czdata, cwdata, pgamma[k], slots, wnum, dim);
		timeutils.stop("Enc sgd step " + k);
	}

	timeutils.start("Enc w gen");
	Cipher* cw = sgd.encwgen(cwdata, wnum, dim);
	timeutils.start("Enc w gen");

	timeutils.start("Dec w");
	double* dw = sgd.decw(secretKey, cw, dim);
	timeutils.start("Dec w");

	sgd.check(dw, zdata, dim, sampledim);
	//-----------------------------------------
	cout << "!!! END TEST SGD !!!" << endl;
}
