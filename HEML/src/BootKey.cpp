#include "BootKey.h"

#include <NTL/ZZ.h>
#include <NTL/ZZX.h>
#include <Params.h>
#include <Ring2Utils.h>
#include <NumUtils.h>
#include <EvaluatorUtils.h>

BootKey::BootKey(Params& params, SchemeAux& aux, SecKey& secretKey, long pBits, long l) : pBits(pBits) {
	ZZ pmod = power2_ZZ(pBits);


	ZZX ex;
	long lpow = 1 << l;
	long lpow2 = lpow << 1;
	long Mover2l = params.N >> l;

	long lk = (1 << (l/2));
	long lm = lpow / lk;

	axGiantRot = new ZZX[lm];
	bxGiantRot = new ZZX[lm];
	for (long i = 0; i < lm; ++i) {
		ZZX spow;
		Ring2Utils::inpower(spow, secretKey.sx, params.rotGroup[i * lk], params.q, params.N);
		Ring2Utils::leftShiftAndEqual(spow, params.logq, params.qq, params.N);
		NumUtils::sampleUniform2(axGiantRot[i], params.N, params.logqq);
		NumUtils::sampleGauss(ex, params.N, params.sigma);
		Ring2Utils::addAndEqual(ex, spow, params.qq, params.N);
		Ring2Utils::mult(bxGiantRot[i], secretKey.sx, axGiantRot[i], params.qq, params.N);
		Ring2Utils::sub(bxGiantRot[i], ex, bxGiantRot[i], params.qq, params.N);
	}

	pvec = new ZZX[lpow];
	pvecInv = new ZZX[lpow];

	CZZ* pdvals = new CZZ[lpow2];
	for (long j = 0; j < lpow; ++j) {
		for (long i = j; i < lpow; ++i) {
			long deg =((2 * params.N - params.rotGroup[i]) * (i - j) * Mover2l) % (2 * params.N);
			CZZ tmp = EvaluatorUtils::evaluateVal(aux.ksiPowsr[deg], aux.ksiPowsi[deg], pBits);
			long idx = (params.rotGroup[i-j] % (2 * lpow2)  - 1) / 2;
					pdvals[idx] = tmp;
					pdvals[lpow2 - idx - 1] = tmp.conjugate();
		}
		for (long i = 0; i < j; ++i) {
			long deg =((2 * params.N - params.rotGroup[i]) * (lpow + i - j) * Mover2l) % (2 * params.N);
			CZZ tmp = EvaluatorUtils::evaluateVal(aux.ksiPowsr[deg], aux.ksiPowsi[deg], pBits);
			long idx = (params.rotGroup[lpow + i - j] % (2 * lpow2) - 1) / 2;
			pdvals[idx] = tmp;
			pdvals[lpow2 - idx - 1] = tmp.conjugate();
		}

		NumUtils::fftSpecialInv(pdvals, lpow2, aux);

		pvec[j].SetLength(params.N);
		long idx = 0;
		long gap = params.N / lpow2;
		for (long i = 0; i < lpow2; ++i) {
			pvec[j].rep[idx] = pdvals[i].r;
			idx += gap;
		}
	}

	for (long j = 0; j < lpow; ++j) {
		for (long i = j; i < lpow; ++i) {
			long deg =(params.rotGroup[i-j] * i * Mover2l) % (2 * params.N);
			CZZ tmp = EvaluatorUtils::evaluateVal(aux.ksiPowsr[deg], aux.ksiPowsi[deg], pBits);
			long idx = (params.rotGroup[i-j] % (2 * lpow2) - 1) / 2;
			pdvals[idx] = tmp;
			pdvals[lpow2 - idx - 1] = tmp.conjugate();
		}
		for (long i = 0; i < j; ++i) {
			long deg = (params.rotGroup[lpow + i - j] * i * Mover2l) % (2 * params.N);
			CZZ tmp = EvaluatorUtils::evaluateVal(aux.ksiPowsr[deg], aux.ksiPowsi[deg], pBits);
			long idx = (params.rotGroup[lpow + i - j] % (2 * lpow2) - 1) / 2;
			pdvals[idx] = tmp;
			pdvals[lpow2 - idx - 1] = tmp.conjugate();
		}
		NumUtils::fftSpecialInv(pdvals, lpow2, aux);

		pvecInv[j].SetLength(params.N);
		long idx = 0;
		long gap = params.N / lpow2;
		for (long i = 0; i < lpow2; ++i) {
			pvecInv[j].rep[idx] = pdvals[i].r;
			idx += gap;
		}
	}
	delete[] pdvals;

	for (long i = 1; i < lm; ++i) {
		for (long j = 0; j < lk; ++j) {
			pvec[j + lk * i] = Ring2Utils::inpower(pvec[j + lk * i], params.rotGroup[lpow - lk * i], pmod, params.N);
			pvecInv[j + lk * i] = Ring2Utils::inpower(pvecInv[j + lk * i], params.rotGroup[lpow - lk * i], pmod, params.N);
		}
	}

	for (long l = 0; l < params.logN; ++l) {
	}
}

