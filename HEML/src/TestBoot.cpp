#include "TestBoot.h"

#include <Cipher.h>
#include <CZZ.h>
#include <EvaluatorUtils.h>
#include <Message.h>
#include <NTL/BasicThreadPool.h>
#include <NTL/RR.h>
#include <NTL/ZZ.h>
#include <NTL/ZZX.h>
#include <Params.h>
#include <PubKey.h>
#include <Scheme.h>
#include <SchemeAlgo.h>
#include <SchemeAux.h>
#include <SecKey.h>
#include <StringUtils.h>
#include <TimeUtils.h>
#include <cstdlib>
#include <iostream>

#include "Boot.h"
#include "BootKey.h"

void TestBoot::testBootstrap() {
	cout << "!!! START TEST BOOTSTRAP ALL !!!" << endl;
	long logq = 1240;
	long nu = 15;
	long msgBits = 39;
	long logq0 = msgBits + nu;
	long logN = 16;
	long logT = 5;
	long logI = 4;
	long logSlots = 0;
	long slots = (1 << logSlots);
	long lkey = logSlots == logN - 1 ? logSlots : logSlots + 1;
	//-----------------------------------------
	TimeUtils timeutils;
	Params params(logN, logq);
	SecKey secretKey(params);
	PubKey publicKey(params, secretKey);
	SchemeAux schemeaux(logN);
	Scheme scheme(params, publicKey, schemeaux);
	SchemeAlgo algo(scheme);

	timeutils.start("Boot Key generating");
	Boot boot(scheme, secretKey);
	boot.generateKey(lkey, logq0 + logI);
	timeutils.stop("Boot Key generated");

	//-----------------------------------------
	SetNumThreads(1);


	CZZ* mvec = EvaluatorUtils::evaluateRandomVals(slots, msgBits);

	//-----------------------------------------
	timeutils.start("Encrypt batch");
	Cipher cipher = scheme.encryptWithBits(mvec, logq0, slots);
	timeutils.stop("Encrypt batch");

	boot.normalizeAndEqual(cipher, params.N);
	cipher.cbits = logq;
	cipher.mod = power2_ZZ(logq);

	timeutils.start("Bootstrap");
	boot.bootstrapAndEqual(cipher, logq0, logT, logI);
	timeutils.stop("Bootstrap");

	timeutils.start("Decrypt batch");
	CZZ* dvec = scheme.decrypt(secretKey, cipher);
	timeutils.stop("Decrypt batch");

	StringUtils::showcompare(mvec, dvec, slots, "m");
	cout << "!!! END TEST BOOTSRTAP ALL !!!" << endl;
}

void TestBoot::testBootstrapOneReal() {
	cout << "!!! START TEST BOOTSTRAP ONE REAL !!!" << endl;
	long logq = 620;
	long logq0 = 29;
	long logN = 15;
	long logT = 2;
	long logI = 4;
	long nu = 6;
	long msgBits = logq0 - nu;

	long logSlots = 0;
	long slots = 1;
	//----------------------------------------

	TimeUtils timeutils;
	timeutils.start("scheme generation");
	Params params(logN, logq);
	SecKey secretKey(params);
	PubKey publicKey(params, secretKey);
	SchemeAux schemeaux(logN);
	Scheme scheme(params, publicKey, schemeaux);
	SchemeAlgo algo(scheme);
	timeutils.stop("scheme generation");

	timeutils.start("Boot Key generating");
	Boot boot(scheme, secretKey);
	timeutils.stop("Boot Key generated");

	//-----------------------------------------
	SetNumThreads(1);

	CZZ* mvec = new CZZ[slots];
	for (long i = 0; i < slots; ++i) {
		ZZ m = RandomBits_ZZ(msgBits);
		mvec[i].r = m;
	}

	//-----------------------------------------
	timeutils.start("Encrypt batch");
	Cipher cipher = scheme.encryptWithBits(mvec, logq0, slots);
	timeutils.stop("Encrypt batch");

	boot.normalizeAndEqual(cipher, params.N);

	cipher.cbits = logq;
	cipher.mod = power2_ZZ(logq);
	timeutils.start("Linear Transform");
	Cipher rotated = cipher;
	for (long i = logSlots; i < scheme.params.logN - 1; ++i) {
		Cipher rot = scheme.leftRotateByPo2(rotated, i);
		scheme.addAndEqual(rotated, rot);
	}
	Cipher cconj = scheme.conjugate(rotated);
	scheme.addAndEqual(rotated, cconj);
	scheme.modSwitchAndEqual(rotated, scheme.params.logN - logSlots);
	timeutils.stop("Linear Transform");

	timeutils.start("remove I part 1");
	Cipher encsin1 = boot.removeIpart(rotated, logq0, logT, logI);
	timeutils.stop("remove I part 1");

	cout << encsin1.cbits << endl;
	timeutils.start("Decrypt batch");
	CZZ* dvec = scheme.decrypt(secretKey, encsin1);
	timeutils.stop("Decrypt batch");

	StringUtils::showcompare(mvec, dvec, slots, "m");
	cout << "!!! END TEST BOOTSRTAP ONE REAL !!!" << endl;
}

void TestBoot::testBoundOfI() {
	cout << "!!! START TEST BOUND OF I !!!" << endl;
	long logq = 200;
	long logq0 = 40;
	long logN = 11;
	long nu = 5;
	long msgBits = logq0 - nu;
	long logSlots = logN - 1;
	long slots = (1 << logSlots);

	//-----------------------------------------
	TimeUtils timeutils;
	Params params(logN, logq);
	SecKey secretKey(params);
	PubKey publicKey(params, secretKey);
	SchemeAux schemeaux(logN);
	Scheme scheme(params, publicKey, schemeaux);
	SchemeAlgo algo(scheme);

	timeutils.start("Boot Key generating");
	Boot boot(scheme, secretKey);
	timeutils.stop("Boot Key generated");

	//-----------------------------------------
	SetNumThreads(1);

	CZZ* mvec = EvaluatorUtils::evaluateRandomVals(slots, msgBits);

	//-----------------------------------------
	timeutils.start("Encrypt batch");
	Cipher cipher = scheme.encryptWithBits(mvec, logq0, slots);
	timeutils.stop("Encrypt batch");

	timeutils.start("Decrypt batch");
	Message msgSmall = scheme.decryptMsg(secretKey, cipher);
	timeutils.stop("Decrypt batch");

	for (long i = 0; i < params.N; ++i) {
		msgSmall.mx.rep[i] = msgSmall.mx.rep[i] % msgSmall.mod;
		if(NumBits(msgSmall.mx.rep[i]) == msgSmall.cbits) msgSmall.mx.rep[i] -= msgSmall.mod;
	}

	boot.normalizeAndEqual(cipher, params.N);

	cipher.cbits = logq;
	cipher.mod = params.q;

	timeutils.start("Decrypt batch");
	Message msgBig = scheme.decryptMsg(secretKey, cipher);
	timeutils.stop("Decrypt batch");

	for (long i = 0; i < params.N; ++i) {
		msgBig.mx.rep[i] = msgBig.mx.rep[i] % msgBig.mod;
		if(NumBits(msgBig.mx.rep[i]) == msgBig.cbits) msgBig.mx.rep[i] -= msgBig.mod;
	}

	ZZ q0 = power2_ZZ(logq0);
	cout << msgBig.mx - msgSmall.mx << endl;
	ZZX pdiv = (msgBig.mx - msgSmall.mx) / q0;
	cout << pdiv << endl;
	ZZ pmax = pdiv.rep[0];
	for (int i = 1; i < deg(pdiv); ++i) {
		if(abs(pmax) < abs(pdiv.rep[i])) {
			pmax = pdiv.rep[i];
		}
	}
	cout << pmax << endl;
	cout << "!!! STOP TEST BOUND OF I !!!" << endl;
}

void TestBoot::test() {
}

