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

void Boot::normalize(Cipher& cipher, long N) {
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

void Boot::linearTransform(Cipher& enclin, Cipher& encx, long size) {
	long logSize = log2(size);
	long k = 1 << (logSize / 2);
	long m = size / k;

	Cipher* encxrotvec = new Cipher[k];
	encxrotvec[0] = encx;

	for (long i = 1; i < k; ++i) {
		encxrotvec[i] = leftBabyRotate(encxrotvec[0], i);
	}

	enclin = scheme.multByPoly(encxrotvec[0], bootKeyMap.at(logSize).pvec[0]);
	for (long j = 1; j < k; ++j) {
		Cipher cij = scheme.multByPoly(encxrotvec[j], bootKeyMap.at(logSize).pvec[j]);
		scheme.addAndEqual(enclin, cij);
	}

	for (long i = 1; i < m; ++i) {
		Cipher ci0 = scheme.multByPoly(encxrotvec[0], bootKeyMap.at(logSize).pvec[k * i]);
		for (long j = 1; j < k; ++j) {
			Cipher cij = scheme.multByPoly(encxrotvec[j], bootKeyMap.at(logSize).pvec[j + k * i]);
			scheme.addAndEqual(ci0, cij);
		}
		leftGiantRotateAndEqual(ci0, logSize, k, i);
		scheme.addAndEqual(enclin, ci0);
	}
	delete[] encxrotvec;
}

void Boot::linearTransformInv(Cipher& res, Cipher& enclin, long size) {
	long logSize = log2(size);
	long k = 1 << (logSize / 2);
	long m = size / k;

	Cipher* enclinrotvec = new Cipher[k];
	enclinrotvec[0] = enclin;
	for (long i = 1; i < k; ++i) {
		enclinrotvec[i] = leftBabyRotate(enclinrotvec[0], i);
	}
	res = scheme.multByPoly(enclinrotvec[0], bootKeyMap.at(logSize).pvecInv[0]);

	for (long j = 1; j < k; ++j) {
		Cipher c0j = scheme.multByPoly(enclinrotvec[j], bootKeyMap.at(logSize).pvecInv[j]);
		scheme.addAndEqual(res, c0j);
	}
	for (long i = 1; i < m; ++i) {
		Cipher ci0 = scheme.multByPoly(enclinrotvec[0], bootKeyMap.at(logSize).pvecInv[k * i]);
		for (long j = 1; j < k; ++j) {
			Cipher cij = scheme.multByPoly(enclinrotvec[j], bootKeyMap.at(logSize).pvecInv[j + k * i]);
			scheme.addAndEqual(ci0, cij);
		}
		leftGiantRotateAndEqual(ci0, logSize, k, i);
		scheme.addAndEqual(res, ci0);
	}
	delete[] enclinrotvec;
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
	Cipher encx13 = scheme.addConst(cipher2, pc);
	scheme.multAndEqual(encx13, tmp);
	scheme.modSwitchAndEqual(encx13, precisionBits); // depth 2

	c = -8*Pi*Pi*Pi*Pi*Pi*Pi*Pi/315;
	pc = EvaluatorUtils::evaluateVal(c, precisionBits);
	tmp = scheme.multByConst(cipher, pc);
	scheme.modSwitchAndEqual(tmp, precisionBits); // depth 1

	c = -21/(2*Pi*Pi);
	pc = EvaluatorUtils::evaluateVal(c, precisionBits);
	Cipher encx57 = scheme.addConst(cipher2, pc);
	scheme.multAndEqual(encx57, tmp);
	scheme.modSwitchAndEqual(encx57, precisionBits); // depth 2
	scheme.multAndEqual(encx57, cipher4);
	scheme.modSwitchAndEqual(encx57, precisionBits); // depth 3

	scheme.modEmbedAndEqual(encx13, precisionBits); // depth 3
	scheme.addAndEqual(encx57, encx13); // depth 3
	return encx57;
}

/**
 * Evaluate sin(2*pi*x) using degree 3 polynomial approximate for sin
 * sin(x) ~ x - x^3/6
 * This process need 2 level down
 */
Cipher Boot::evaluateEncSin2pix3(Cipher& enclinx, long precisionBits) {

	Cipher encx2 = scheme.square(enclinx); //depth 1
	scheme.modSwitchAndEqual(encx2, precisionBits);

	RR c = -4*Pi*Pi*Pi/3;
	ZZ pc = EvaluatorUtils::evaluateVal(c, precisionBits);
	Cipher tmp = scheme.multByConst(enclinx, pc);
	scheme.modSwitchAndEqual(tmp, precisionBits); // depth 1

	c = -3/(2*Pi*Pi);
	pc = EvaluatorUtils::evaluateVal(c, precisionBits);
	Cipher res = scheme.addConst(encx2, pc);
	scheme.multAndEqual(res, tmp);
	scheme.modSwitchAndEqual(res, precisionBits); // depth 2

	return res;
}

/**
 * Evaluate cos(2*pi*x) using degree 6 polynomial approximate for cos
 * cos(x) ~ 1 - x^2/2 + x^4/24 - x^6/720
 * */
Cipher Boot::evaluateEncCos2pix6(Cipher& enclinx, long precisionBits) {

	Cipher encx2 = scheme.square(enclinx); //depth 1
	scheme.modSwitchAndEqual(encx2, precisionBits);

	Cipher encx4 = scheme.square(encx2); //depth 2
	scheme.modSwitchAndEqual(encx4, precisionBits);

	RR c = -1/(2*Pi*Pi);
	ZZ pc = EvaluatorUtils::evaluateVal(c, precisionBits);
	Cipher encx02 = scheme.addConst(encx2, pc);

	c = -2*Pi*Pi;
	pc = EvaluatorUtils::evaluateVal(c, precisionBits);
	scheme.multByConstAndEqual(encx02, pc);
	scheme.modSwitchAndEqual(encx02, precisionBits); // depth 2

	c = -15/(2*Pi*Pi);
	pc = EvaluatorUtils::evaluateVal(c, precisionBits);
	Cipher encx46 = scheme.addConst(encx2, pc);

	c = -4*Pi*Pi*Pi*Pi*Pi*Pi/45;
	pc = EvaluatorUtils::evaluateVal(c, precisionBits);
	scheme.multByConstAndEqual(encx46, pc);
	scheme.modSwitchAndEqual(encx46, precisionBits); // depth 2

	scheme.multAndEqual(encx46, encx4);
	scheme.modSwitchAndEqual(encx46, precisionBits); // depth 3

	scheme.modEmbedAndEqual(encx02, precisionBits); // depth 3
	scheme.addAndEqual(encx46, encx02); // depth 3
	return encx46;
}

Cipher Boot::evaluateEncCos2pix2(Cipher& enclinx, long precisionBits) {

	Cipher encx2 = scheme.square(enclinx); //depth 1
	scheme.modSwitchAndEqual(encx2, precisionBits);

	RR c = -1/(2*Pi*Pi);
	ZZ pc = EvaluatorUtils::evaluateVal(c, precisionBits);
	Cipher res = scheme.addConst(encx2, pc);

	c = -2*Pi*Pi;
	pc = EvaluatorUtils::evaluateVal(c, precisionBits);
	scheme.multByConstAndEqual(res, pc);
	scheme.modSwitchAndEqual(res, precisionBits); // depth 2

	return res;
}

Cipher Boot::evaluateEncSin2x(Cipher& encSinx, Cipher& encCosx, long precisionBits) {
	Cipher res = scheme.mult(encSinx, encCosx);
	scheme.doubleAndEqual(res);
	scheme.modSwitchAndEqual(res, precisionBits);
	return res;
}

Cipher Boot::evaluateEncCos2x(Cipher& encSinx, Cipher& encCosx, long precisionBits) {
	Cipher encSub = scheme.sub(encCosx, encSinx);
	Cipher encAdd = scheme.add(encCosx, encSinx);
	scheme.multAndEqual(encAdd, encSub);
	scheme.modSwitchAndEqual(encAdd, precisionBits);
	return encAdd;
}

Cipher Boot::removeIpart(Cipher& cipher, long logq0, long logT, long logI) {
	Cipher cms = scheme.modSwitch(cipher, logT);

	Cipher encsinx = evaluateEncSin2pix7(cms, logq0 + logI);
	Cipher enccosx = evaluateEncCos2pix6(cms, logq0 + logI);

	Cipher encsin2x, enccos2x;
	for (long i = 0; i < logI + logT; ++i) {
		encsin2x = evaluateEncSin2x(encsinx, enccosx, logq0 + logI);
		enccos2x = evaluateEncCos2x(encsinx, enccosx, logq0 + logI);
		encsinx = encsin2x;
		enccosx = enccos2x;
	}
	ZZ temp = EvaluatorUtils::evaluateVal(1/(2*Pi), logq0 + logI);
	scheme.multByConstAndEqual(encsinx, temp);
	scheme.modSwitchAndEqual(encsinx, logq0 + 2 * logI);
	return encsinx;
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
	}
}

