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

	double** wdata = new double*[wnum];
	for (long l = 0; l < wnum; ++l) {
		wdata[l] = new double[dim];
		for (long i = 0; i < dim; ++i) {
			wdata[l][i] = 1.0 - (double)rand() / RAND_MAX; // change to good initial w choice
		}
	}

	long iter = 1000;
	double lambda = 2.0;
	double* alpha = new double[iter];
	for (long k = 0; k < iter; ++k) {
		alpha[k] = 1.0 / (k + 1);
	}

	double* w = sgd.sgd(iter, wnum, wdata, zdata, alpha, lambda, dim, sampledim);

	sgd.check(w, zdata, dim, sampledim);

	//-----------------------------------------

	timeutils.start("Encrypting zdata");
	Cipher* zciphers = sgd.encryptzdata(zdata, slots, wnum, dim, sampledim, params.p);
	timeutils.stop("Encrypting zdata");

	timeutils.start("Encrypting wdata");
	Cipher* wciphers = sgd.encryptwdata(wdata, slots, wnum, dim, sampledim, params.logp);
	timeutils.stop("Encrypting wdata");


	iter = 100;
	//-----------------------------------------
	for (long k = 0; k < iter; ++k) {
		timeutils.start("SGD");
		Cipher* cgrad = sgd.cipherGradient(zciphers, wciphers, dim, slots, wnum);
		timeutils.stop("SGD");
		//-----------------------------------------

		timeutils.start("change grad");
		NTL_EXEC_RANGE(dim, first, last);
		for (long i = first; i < last; ++i) {
		}
		NTL_EXEC_RANGE_END;

		timeutils.stop("change grad");
	}

	//-----------------------------------------
	cout << "!!! END TEST SGD !!!" << endl;
}
