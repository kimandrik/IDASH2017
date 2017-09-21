#include "Boot.h"

#include <CZZ.h>
#include <EvaluatorUtils.h>
#include <NTL/ZZX.h>
#include <NumUtils.h>
#include <Params.h>
#include <Ring2Utils.h>
#include <cmath>
#include <initializer_list>
#include <utility>

Boot::Boot(Scheme& scheme, SecKey& secretKey) : scheme(scheme), secretKey(secretKey) {
	long k = 1 << ((scheme.params.logN - 1) / 2);
	axBabyRot = new ZZX[k];
	bxBabyRot = new ZZX[k];

	ZZX ex;
	for (long i = 0; i < k; ++i) {
		ZZX spow;
		Ring2Utils::inpower(spow, secretKey.sx, scheme.params.rotGroup[i], scheme.params.q, scheme.params.N);
		Ring2Utils::leftShiftAndEqual(spow, scheme.params.logq, scheme.params.qq, scheme.params.N);
		NumUtils::sampleUniform2(axBabyRot[i], scheme.params.N, scheme.params.logqq);
		NumUtils::sampleGauss(ex, scheme.params.N, scheme.params.sigma);
		Ring2Utils::addAndEqual(ex, spow, scheme.params.qq, scheme.params.N);
		Ring2Utils::mult(bxBabyRot[i], secretKey.sx, axBabyRot[i], scheme.params.qq, scheme.params.N);
		Ring2Utils::sub(bxBabyRot[i], ex, bxBabyRot[i], scheme.params.qq, scheme.params.N);
	}
}

Boot::~Boot() {
}

void Boot::generateKey(long l, long pBits) {
	if(bootKeyMap.find(l) == bootKeyMap.end()) {
		bootKeyMap.insert(pair<long, BootKey>(l, BootKey(scheme.params, scheme.aux, secretKey, pBits, l)));
	}
}

void Boot::normalizeAndEqual(Cipher& cipher, long N) {
	for (int i = 0; i < N; ++i) {
		if(NumBits(cipher.ax.rep[i]) == cipher.cbits) cipher.ax.rep[i] -= cipher.mod;
		if(NumBits(cipher.bx.rep[i]) == cipher.cbits) cipher.bx.rep[i] -= cipher.mod;
	}
}

Cipher Boot::leftBabyRotate(Cipher& cipher, long i) {
	ZZ Pmod = cipher.mod << scheme.params.logq;

	ZZX bxrot, bxres, axres;

	Ring2Utils::inpower(bxrot, cipher.bx, scheme.params.rotGroup[i], cipher.mod, scheme.params.N);
	Ring2Utils::inpower(bxres, cipher.ax, scheme.params.rotGroup[i], cipher.mod, scheme.params.N);

	Ring2Utils::mult(axres, bxres, axBabyRot[i], Pmod, scheme.params.N);
	Ring2Utils::multAndEqual(bxres, bxBabyRot[i], Pmod, scheme.params.N);

	Ring2Utils::rightShiftAndEqual(axres, scheme.params.logq, scheme.params.N);
	Ring2Utils::rightShiftAndEqual(bxres, scheme.params.logq, scheme.params.N);

	Ring2Utils::addAndEqual(bxres, bxrot, cipher.mod, scheme.params.N);
	return Cipher(axres, bxres, cipher.mod, cipher.cbits, cipher.slots);
}

void Boot::leftGiantRotateAndEqual(Cipher& cipher, long l, long k, long i) {
	ZZ Pmod = cipher.mod << scheme.params.logq;

	ZZX bxrot, bxres, axres;

	Ring2Utils::inpower(bxrot, cipher.bx, scheme.params.rotGroup[k * i], cipher.mod, scheme.params.N);
	Ring2Utils::inpower(bxres, cipher.ax, scheme.params.rotGroup[k * i], cipher.mod, scheme.params.N);

	Ring2Utils::mult(axres, bxres, bootKeyMap.at(l).axGiantRot[i], Pmod, scheme.params.N);
	Ring2Utils::multAndEqual(bxres, bootKeyMap.at(l).bxGiantRot[i], Pmod, scheme.params.N);

	Ring2Utils::rightShiftAndEqual(axres, scheme.params.logq, scheme.params.N);
	Ring2Utils::rightShiftAndEqual(bxres, scheme.params.logq, scheme.params.N);

	Ring2Utils::addAndEqual(bxres, bxrot, cipher.mod, scheme.params.N);

	cipher.ax = axres;
	cipher.bx = bxres;
}

void Boot::linearTransform(Cipher& res, Cipher& cipher, long size) {
	long logSize = log2(size);
	long k = 1 << (logSize / 2);
	long m = size / k;

	Cipher* cipherRotVec = new Cipher[k];
	cipherRotVec[0] = cipher;

	for (long i = 1; i < k; ++i) {
		cipherRotVec[i] = leftBabyRotate(cipherRotVec[0], i);
	}

	res = scheme.multByPoly(cipherRotVec[0], bootKeyMap.at(logSize).pvec[0]);
	for (long j = 1; j < k; ++j) {
		Cipher cij = scheme.multByPoly(cipherRotVec[j], bootKeyMap.at(logSize).pvec[j]);
		scheme.addAndEqual(res, cij);
	}

	for (long i = 1; i < m; ++i) {
		Cipher ci0 = scheme.multByPoly(cipherRotVec[0], bootKeyMap.at(logSize).pvec[k * i]);
		for (long j = 1; j < k; ++j) {
			Cipher cij = scheme.multByPoly(cipherRotVec[j], bootKeyMap.at(logSize).pvec[j + k * i]);
			scheme.addAndEqual(ci0, cij);
		}
		leftGiantRotateAndEqual(ci0, logSize, k, i);
		scheme.addAndEqual(res, ci0);
	}
	delete[] cipherRotVec;
}

void Boot::linearTransformAndEqual(Cipher& cipher, long size) {
	//TODO
}

void Boot::linearTransformInv(Cipher& res, Cipher& cipher, long size) {
	long logSize = log2(size);
	long k = 1 << (logSize / 2);
	long m = size / k;

	Cipher* cipherRotVec = new Cipher[k];
	cipherRotVec[0] = cipher;
	for (long i = 1; i < k; ++i) {
		cipherRotVec[i] = leftBabyRotate(cipherRotVec[0], i);
	}
	res = scheme.multByPoly(cipherRotVec[0], bootKeyMap.at(logSize).pvecInv[0]);

	for (long j = 1; j < k; ++j) {
		Cipher c0j = scheme.multByPoly(cipherRotVec[j], bootKeyMap.at(logSize).pvecInv[j]);
		scheme.addAndEqual(res, c0j);
	}
	for (long i = 1; i < m; ++i) {
		Cipher ci0 = scheme.multByPoly(cipherRotVec[0], bootKeyMap.at(logSize).pvecInv[k * i]);
		for (long j = 1; j < k; ++j) {
			Cipher cij = scheme.multByPoly(cipherRotVec[j], bootKeyMap.at(logSize).pvecInv[j + k * i]);
			scheme.addAndEqual(ci0, cij);
		}
		leftGiantRotateAndEqual(ci0, logSize, k, i);
		scheme.addAndEqual(res, ci0);
	}
	delete[] cipherRotVec;
}

void Boot::linearTransformInvAndEqual(Cipher& cipher, long size) {
	//TODO
}

/**
 * Evaluate sin(2*pi*x) using degree 7 polynomial approximate for sin
 * sin(x) ~ x - x^3/6 + x^5/120 - x^7/5040
 * This process need 3 level down
 */
Cipher Boot::evaluateEncSin2pix7(Cipher& cipher, long precisionBits) {

	Cipher cipher2 = scheme.square(cipher); //depth 1
	scheme.modSwitchAndEqual(cipher2, precisionBits);
	Cipher cipher4 = scheme.square(cipher2); //depth 2
	scheme.modSwitchAndEqual(cipher4, precisionBits);
	RR c = -4*Pi*Pi*Pi/3;
	ZZ pc = EvaluatorUtils::evaluateVal(c, precisionBits);
	Cipher tmp = scheme.multByConst(cipher, pc);
	scheme.modSwitchAndEqual(tmp, precisionBits); // depth 1

	c = -3/(2*Pi*Pi);
	pc = EvaluatorUtils::evaluateVal(c, precisionBits);
	Cipher cipher13 = scheme.addConst(cipher2, pc);
	scheme.multAndEqual(cipher13, tmp);
	scheme.modSwitchAndEqual(cipher13, precisionBits); // depth 2

	c = -8*Pi*Pi*Pi*Pi*Pi*Pi*Pi/315;
	pc = EvaluatorUtils::evaluateVal(c, precisionBits);
	tmp = scheme.multByConst(cipher, pc);
	scheme.modSwitchAndEqual(tmp, precisionBits); // depth 1

	c = -21/(2*Pi*Pi);
	pc = EvaluatorUtils::evaluateVal(c, precisionBits);
	Cipher cipher57 = scheme.addConst(cipher2, pc);
	scheme.multAndEqual(cipher57, tmp);
	scheme.modSwitchAndEqual(cipher57, precisionBits); // depth 2
	scheme.multAndEqual(cipher57, cipher4);
	scheme.modSwitchAndEqual(cipher57, precisionBits); // depth 3

	scheme.modEmbedAndEqual(cipher13, precisionBits); // depth 3
	scheme.addAndEqual(cipher57, cipher13); // depth 3
	return cipher57;
}

/**
 * Evaluate sin(2*pi*x) using degree 3 polynomial approximate for sin
 * sin(x) ~ x - x^3/6
 * This process need 2 level down
 */
Cipher Boot::evaluateEncSin2pix3(Cipher& cipher, long precisionBits) {

	Cipher cipher2 = scheme.square(cipher); //depth 1
	scheme.modSwitchAndEqual(cipher2, precisionBits);

	RR c = -4*Pi*Pi*Pi/3;
	ZZ pc = EvaluatorUtils::evaluateVal(c, precisionBits);
	Cipher tmp = scheme.multByConst(cipher, pc);
	scheme.modSwitchAndEqual(tmp, precisionBits); // depth 1

	c = -3/(2*Pi*Pi);
	pc = EvaluatorUtils::evaluateVal(c, precisionBits);
	Cipher res = scheme.addConst(cipher2, pc);
	scheme.multAndEqual(res, tmp);
	scheme.modSwitchAndEqual(res, precisionBits); // depth 2

	return res;
}

/**
 * Evaluate cos(2*pi*x) using degree 6 polynomial approximate for cos
 * cos(x) ~ 1 - x^2/2 + x^4/24 - x^6/720
 * */
Cipher Boot::evaluateEncCos2pix6(Cipher& cipher, long precisionBits) {

	Cipher cipher2 = scheme.square(cipher); //depth 1
	scheme.modSwitchAndEqual(cipher2, precisionBits);

	Cipher cipher4 = scheme.square(cipher2); //depth 2
	scheme.modSwitchAndEqual(cipher4, precisionBits);

	RR c = -1/(2*Pi*Pi);
	ZZ pc = EvaluatorUtils::evaluateVal(c, precisionBits);
	Cipher cipher02 = scheme.addConst(cipher2, pc);

	c = -2*Pi*Pi;
	pc = EvaluatorUtils::evaluateVal(c, precisionBits);
	scheme.multByConstAndEqual(cipher02, pc);
	scheme.modSwitchAndEqual(cipher02, precisionBits); // depth 2

	c = -15/(2*Pi*Pi);
	pc = EvaluatorUtils::evaluateVal(c, precisionBits);
	Cipher cipher46 = scheme.addConst(cipher2, pc);

	c = -4*Pi*Pi*Pi*Pi*Pi*Pi/45;
	pc = EvaluatorUtils::evaluateVal(c, precisionBits);
	scheme.multByConstAndEqual(cipher46, pc);
	scheme.modSwitchAndEqual(cipher46, precisionBits); // depth 2

	scheme.multAndEqual(cipher46, cipher4);
	scheme.modSwitchAndEqual(cipher46, precisionBits); // depth 3

	scheme.modEmbedAndEqual(cipher02, precisionBits); // depth 3
	scheme.addAndEqual(cipher46, cipher02); // depth 3
	return cipher46;
}

Cipher Boot::evaluateEncCos2pix2(Cipher& cipher, long precisionBits) {
	Cipher cipher2 = scheme.square(cipher); //depth 1
	scheme.modSwitchAndEqual(cipher2, precisionBits);
	RR c = -1/(2*Pi*Pi);
	ZZ pc = EvaluatorUtils::evaluateVal(c, precisionBits);
	Cipher res = scheme.addConst(cipher2, pc);

	c = -2*Pi*Pi;
	pc = EvaluatorUtils::evaluateVal(c, precisionBits);
	scheme.multByConstAndEqual(res, pc);
	scheme.modSwitchAndEqual(res, precisionBits); // depth 2

	return res;
}

Cipher Boot::evaluateEncSin2x(Cipher& cipherSin, Cipher& cipherCos, long precisionBits) {
	Cipher res = scheme.mult(cipherSin, cipherCos);
	scheme.doubleAndEqual(res);
	scheme.modSwitchAndEqual(res, precisionBits);
	return res;
}

Cipher Boot::evaluateEncCos2x(Cipher& cipherSin, Cipher& cipherCos, long precisionBits) {
	Cipher encSub = scheme.sub(cipherCos, cipherSin);
	Cipher encAdd = scheme.add(cipherCos, cipherSin);
	scheme.multAndEqual(encAdd, encSub);
	scheme.modSwitchAndEqual(encAdd, precisionBits);
	return encAdd;
}

Cipher Boot::removeIpart(Cipher& cipher, long logq0, long logT, long logI) {
	Cipher cms = scheme.modSwitch(cipher, logT);

	Cipher cipherSin = evaluateEncSin2pix7(cms, logq0 + logI);
	Cipher cipherCos = evaluateEncCos2pix6(cms, logq0 + logI);

	Cipher cipherSin2, cipherCos2;
	for (long i = 0; i < logI + logT; ++i) {
		cipherSin2 = evaluateEncSin2x(cipherSin, cipherCos, logq0 + logI);
		cipherCos2 = evaluateEncCos2x(cipherSin, cipherCos, logq0 + logI);
		cipherSin = cipherSin2;
		cipherCos = cipherCos2;
	}
	ZZ temp = EvaluatorUtils::evaluateVal(1/(2*Pi), logq0 + logI);
	scheme.multByConstAndEqual(cipherSin, temp);
	scheme.modSwitchAndEqual(cipherSin, logq0 + 2 * logI);
	return cipherSin;
}

void Boot::removeIpartAndEqual(Cipher& cipher, long logq0, long logT, long logI) {
	//TODO
}

Cipher Boot::bootstrap(Cipher& cipher, long logq0, long logT, long logI) {
	long logSlots = log2(cipher.slots);
	if(logSlots == scheme.params.logN - 1) {
		Cipher clinOdd, clinEven, cresEven, cresOdd;
		linearTransform(clinEven, cipher, scheme.params.N / 2);
		Cipher cshift1 = scheme.multByMonomial(cipher, 2 * scheme.params.N - 1);
		linearTransform(clinOdd, cshift1, scheme.params.N / 2);

		Cipher clinEvenConj = scheme.conjugate(clinEven);
		scheme.addAndEqual(clinEven, clinEvenConj);
		scheme.modSwitchAndEqual(clinEven, logq0 + logI + logSlots);

		Cipher clinOddConj = scheme.conjugate(clinOdd);
		scheme.addAndEqual(clinOdd, clinOddConj);
		scheme.modSwitchAndEqual(clinOdd, logq0 + logI + logSlots);

		Cipher csinEven = removeIpart(clinEven, logq0, logT, logI);
		Cipher csinOdd = removeIpart(clinOdd, logq0, logT, logI);

		linearTransformInv(cresEven, csinEven, scheme.params.N / 2);
		linearTransformInv(cresOdd, csinOdd, scheme.params.N / 2);

		scheme.multByMonomialAndEqual(cresOdd, 1);
		scheme.addAndEqual(cresEven, cresOdd);
		scheme.modSwitchAndEqual(cresEven, logq0 + logI);

		return cresEven;
	} else {
		Cipher clinEven, cresEven;

		Cipher rotated = cipher;
		for (long i = logSlots; i < scheme.params.logN - 1; ++i) {
			Cipher rot = scheme.leftRotateByPo2(rotated, i);
			scheme.addAndEqual(rotated, rot);
		}
		scheme.modSwitchAndEqual(rotated, scheme.params.logN - 1 - logSlots);
		linearTransform(clinEven, rotated, cipher.slots * 2);

		Cipher clinEvenConj = scheme.conjugate(clinEven);
		scheme.addAndEqual(clinEven, clinEvenConj);
		scheme.modSwitchAndEqual(clinEven, logq0 + logI + logSlots + 2);

		Cipher csinEven = removeIpart(clinEven, logq0, logT, logI);
		linearTransformInv(cresEven, csinEven, cipher.slots * 2);
		scheme.modSwitchAndEqual(cresEven, logq0 + logI);
		return cresEven;
	}
}

void Boot::bootstrapAndEqual(Cipher& cipher, long logq0, long logT, long logI) {
	//TODO
}

