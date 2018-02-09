#include "CipherGD.h"
#include "GD.h"

#include <Ciphertext.h>
#include <EvaluatorUtils.h>
#include <NTL/BasicThreadPool.h>
#include <NTL/RR.h>
#include <NTL/ZZ.h>

void CipherGD::encxyData(Ciphertext* cxyData, double** xyData, long slots, long factorDim, long sampleDim, long batch, long cnum, long wBits) {
	complex<double>* pxyData = new complex<double>[slots];
	for (long i = 0; i < cnum - 1; ++i) {
		for (long j = 0; j < sampleDim; ++j) {
			for (long l = 0; l < batch; ++l) {
				pxyData[batch * j + l].real(xyData[j][batch * i + l]);
			}
		}
		cxyData[i] = scheme.encrypt(pxyData, slots, wBits, scheme.context.logQ);
	}

	long rest = factorDim - batch * (cnum - 1);
	for (long j = 0; j < sampleDim; ++j) {
		for (long l = 0; l < rest; ++l) {
			pxyData[batch * j + l].real(xyData[j][batch * (cnum - 1) + l]);
		}
		for (long l = rest; l < batch; ++l) {
			pxyData[batch * j + l] = 0;
		}
	}
	cxyData[cnum - 1] = scheme.encrypt(pxyData, slots, wBits, scheme.context.logQ);

	delete[] pxyData;
}

void CipherGD::encwData(Ciphertext* cwData, Ciphertext* cxyData, long cnum, long sBits, long bBits) {
	NTL_EXEC_RANGE(cnum, first, last);
	for (long i = first; i < last; ++i) {
		cwData[i] = cxyData[i];
		for (long l = bBits; l < sBits; ++l) {
			Ciphertext rot = scheme.leftRotateByPo2(cwData[i], l);
			scheme.addAndEqual(cwData[i], rot);
		}
		scheme.divByPo2AndEqual(cwData[i], sBits - bBits);
	}
	NTL_EXEC_RANGE_END;
}

ZZX CipherGD::generateAuxPoly(long slots, long batch, long pBits) {
	complex<double>* pvals = new complex<double>[slots];
	for (long j = 0; j < slots; j += batch) {
		pvals[j].real(1.0);
	}
	ZZX msg = scheme.context.encode(pvals, slots, pBits);
	delete[] pvals;
	return msg;
}

Ciphertext CipherGD::encIP(Ciphertext* cxyData, Ciphertext* cwData, ZZX& poly, long cnum, long bBits, long wBits, long pBits) {
	Ciphertext cip = scheme.modDownTo(cxyData[0], cwData[0].logq);
	scheme.multAndEqual(cip, cwData[0]); // xy * w
	for (long l = 0; l < bBits; ++l) {
		Ciphertext rot = scheme.leftRotateByPo2(cip, l);
		scheme.addAndEqual(cip, rot);
	}

	NTL_EXEC_RANGE(cnum, first, last);
	for (long i = first + 1; i < last; ++i) {
		Ciphertext tmp = scheme.modDownTo(cxyData[i], cwData[i].logq);
		scheme.multAndEqual(tmp, cwData[i]); // xy * w
		for (long l = 0; l < bBits; ++l) {
			Ciphertext rot = scheme.leftRotateByPo2(tmp, l);
			scheme.addAndEqual(tmp, rot);
		}
		scheme.addAndEqual(cip, tmp);
	}
	NTL_EXEC_RANGE_END

	scheme.reScaleByAndEqual(cip, wBits);
	scheme.multByPolyAndEqual(cip, poly, pBits);
	for (long l = 0; l < bBits; ++l) {
		Ciphertext tmp = scheme.rightRotateByPo2(cip, l);
		scheme.addAndEqual(cip, tmp);
	}
	return cip;
}

void CipherGD::encSigmoid(long approxDeg, Ciphertext* cxyData, Ciphertext* cgrad, Ciphertext& cip, long cnum, double gamma, long sBits, long bBits, long wBits, long aBits) {
	Ciphertext cip2 = scheme.square(cip);
	scheme.reScaleByAndEqual(cip2, wBits);

	if(approxDeg == 3) {
		scheme.addConstAndEqual(cip2, degree3[1] / degree3[2], wBits - 2 * aBits);
		NTL_EXEC_RANGE(cnum, first, last);
		for (long i = first; i < last; ++i) {
			cgrad[i] = scheme.multByConst(cxyData[i], gamma  * degree3[2], wBits + 3 * aBits);
			scheme.reScaleByAndEqual(cgrad[i], wBits);
			scheme.modDownToAndEqual(cgrad[i], cip.logq);
			scheme.multAndEqual(cgrad[i], cip);
			scheme.reScaleByAndEqual(cgrad[i], wBits);
			scheme.multAndEqual(cgrad[i], cip2);
			scheme.reScaleByAndEqual(cgrad[i], wBits);

			Ciphertext tmp = scheme.multByConst(cxyData[i], gamma * degree3[0], wBits);
			scheme.reScaleByAndEqual(tmp, wBits);
			scheme.modDownToAndEqual(tmp, cgrad[i].logq);
			scheme.addAndEqual(cgrad[i], tmp);
		}
		NTL_EXEC_RANGE_END;
	} else if (approxDeg == 5) {
		Ciphertext cip4 = scheme.square(cip2);
		scheme.multByConstAndEqual(cip2, degree5[2] / degree5[3], wBits - 2 * aBits);
		scheme.reScaleByAndEqual(cip2, wBits);
		scheme.reScaleByAndEqual(cip4, wBits);
		scheme.addAndEqual(cip4, cip2);
		scheme.addConstAndEqual(cip4, degree5[1] / degree5[3], wBits - 4 * aBits);

		NTL_EXEC_RANGE(cnum, first, last);
		for (long i = first; i < last; ++i) {
			cgrad[i] = scheme.multByConst(cxyData[i], gamma * degree5[3], wBits + 5 * aBits);
			scheme.reScaleByAndEqual(cgrad[i], wBits);
			scheme.modDownToAndEqual(cgrad[i], cip.logq);
			scheme.multAndEqual(cgrad[i], cip);
			scheme.reScaleByAndEqual(cgrad[i], wBits);
			scheme.modDownToAndEqual(cgrad[i], cip4.logq);
			scheme.multAndEqual(cgrad[i], cip4);
			scheme.reScaleByAndEqual(cgrad[i], wBits);

			Ciphertext tmp = scheme.multByConst(cxyData[i], gamma * degree5[0], wBits);
			scheme.reScaleByAndEqual(tmp, wBits);
			scheme.modDownToAndEqual(tmp, cgrad[i].logq);
			scheme.addAndEqual(cgrad[i], tmp);
		}
		NTL_EXEC_RANGE_END;
	} else {
		Ciphertext cip4 = scheme.square(cip2);
		scheme.reScaleByAndEqual(cip4, wBits);
		Ciphertext cip2c = scheme.multByConst(cip2, degree7[3] / degree7[4], wBits - 2 * aBits);
		scheme.reScaleByAndEqual(cip2c, wBits);
		scheme.addAndEqual(cip4, cip2c);
		scheme.addConstAndEqual(cip4, degree7[2] / degree7[4], wBits - 4 * aBits);

		NTL_EXEC_RANGE(cnum, first, last);
		for (long i = first; i < last; ++i) {
			Ciphertext tmp = scheme.multByConst(cxyData[i], gamma * degree7[1], wBits + aBits);
			scheme.reScaleByAndEqual(tmp, wBits);
			scheme.modDownToAndEqual(tmp, cip.logq);
			scheme.multAndEqual(tmp, cip);
			scheme.reScaleByAndEqual(tmp, wBits);

			cgrad[i] = scheme.multByConst(cxyData[i], gamma * degree7[0], wBits);
			scheme.reScaleByAndEqual(cgrad[i], wBits);
			scheme.modDownToAndEqual(cgrad[i], tmp.logq);
			scheme.addAndEqual(tmp, cgrad[i]);

			cgrad[i] = scheme.multByConst(cxyData[i], gamma * degree7[4], wBits + 7 * aBits);
			scheme.reScaleByAndEqual(cgrad[i], wBits);
			scheme.modDownToAndEqual(cgrad[i], cip.logq);
			scheme.multAndEqual(cgrad[i], cip);
			scheme.reScaleByAndEqual(cgrad[i], wBits);
			scheme.multAndEqual(cgrad[i], cip2);
			scheme.reScaleByAndEqual(cgrad[i], wBits);
			scheme.multAndEqual(cgrad[i], cip4);
			scheme.reScaleByAndEqual(cgrad[i], wBits);

			scheme.modDownToAndEqual(tmp, cgrad[i].logq);
			scheme.addAndEqual(cgrad[i], tmp);
		}
		NTL_EXEC_RANGE_END;

	}

	NTL_EXEC_RANGE(cnum, first, last);
	for (long i = first; i < last; ++i) {
		for (long l = bBits; l < sBits; ++l) {
			Ciphertext tmp = scheme.leftRotateByPo2(cgrad[i], l);
			scheme.addAndEqual(cgrad[i], tmp);
		}
	}
	NTL_EXEC_RANGE_END;
}

void CipherGD::encLGDstep(Ciphertext* cwData, Ciphertext* cgrad, long cnum) {
	NTL_EXEC_RANGE(cnum, first, last);
	for (long i = first; i < last; ++i) {
		scheme.modDownToAndEqual(cwData[i], cgrad[i].logq);
		scheme.subAndEqual(cwData[i], cgrad[i]);
	}
	NTL_EXEC_RANGE_END;
}

void CipherGD::encMLGDstep(Ciphertext* cwData, Ciphertext* cvData, Ciphertext* cgrad, double eta, long cnum, long pBits) {
	NTL_EXEC_RANGE(cnum, first, last);
	for (long i = first; i < last; ++i) {
		scheme.multByConstAndEqual(cvData[i], eta, pBits);
		scheme.reScaleByAndEqual(cvData[i], pBits);
		scheme.modDownToAndEqual(cvData[i], cgrad[i].logq);
		scheme.addAndEqual(cvData[i], cgrad[i]);
		scheme.modDownToAndEqual(cwData[i], cvData[i].logq);
		scheme.subAndEqual(cwData[i], cvData[i]);
	}
	NTL_EXEC_RANGE_END;
}

void CipherGD::encNLGDstep(Ciphertext* cwData, Ciphertext* cvData, Ciphertext* cgrad, double eta, long cnum, long pBits) {
	NTL_EXEC_RANGE(cnum, first, last);
	for (long i = first; i < last; ++i) {
		scheme.modDownToAndEqual(cvData[i], cgrad[i].logq);
		Ciphertext ctmpw = scheme.sub(cvData[i], cgrad[i]);
		cvData[i] = scheme.multByConst(ctmpw, 1. - eta, pBits);
		scheme.reScaleByAndEqual(cvData[i], pBits);
		scheme.multByConstAndEqual(cwData[i], eta, pBits);
		scheme.reScaleByAndEqual(cwData[i], pBits);
		scheme.modDownToAndEqual(cwData[i], cvData[i].logq);
		scheme.addAndEqual(cvData[i], cwData[i]);
		cwData[i] = ctmpw;
	}
	NTL_EXEC_RANGE_END;
}

void CipherGD::encLGDiteration(long approxDeg, Ciphertext* cxyData, Ciphertext* cwData, ZZX& poly, long cnum, double gamma, long sBits, long bBits, long wBits, long pBits, long aBits) {
 	Ciphertext* cgrad = new Ciphertext[cnum];
	Ciphertext cip = encIP(cxyData, cwData, poly, cnum, bBits, wBits, pBits);
	encSigmoid(approxDeg, cxyData, cgrad, cip, cnum, gamma, sBits, bBits, wBits, aBits);
	encLGDstep(cwData, cgrad, cnum);
	delete[] cgrad;
}

void CipherGD::encMLGDiteration(long approxDeg, Ciphertext* cxyData, Ciphertext* cwData, Ciphertext* cvData, ZZX& poly, long cnum, double gamma, double eta, long sBits, long bBits, long wBits, long pBits, long aBits) {
 	Ciphertext* cgrad = new Ciphertext[cnum];
	Ciphertext cip = encIP(cxyData, cwData, poly, cnum, bBits, wBits, pBits);
	encSigmoid(approxDeg, cxyData, cgrad, cip, cnum, gamma, sBits, bBits, wBits, aBits);
	encMLGDstep(cwData, cvData, cgrad, eta, cnum, wBits);
	delete[] cgrad;
}

void CipherGD::encNLGDiteration(long approxDeg, Ciphertext* cxyData, Ciphertext* cwData, Ciphertext* cvData, ZZX& poly, long cnum, double gamma, double eta, long sBits, long bBits, long wBits, long pBits, long aBits) {
 	Ciphertext* cgrad = new Ciphertext[cnum];
	Ciphertext cip = encIP(cxyData, cvData, poly, cnum, bBits, wBits, pBits);
	scheme.reScaleByAndEqual(cip, pBits + aBits);
	encSigmoid(approxDeg, cxyData, cgrad, cip, cnum, gamma, sBits, bBits, wBits, aBits);
	encNLGDstep(cwData, cvData, cgrad, eta, cnum, pBits);
	delete[] cgrad;
}

void CipherGD::decwData(double* wData, Ciphertext* cwData, long factorDim, long batch, long cnum, long wBits) {
	for (long i = 0; i < (cnum - 1); ++i) {
		complex<double>* dcw = scheme.decrypt(secretKey, cwData[i]);
		for (long j = 0; j < batch; ++j) {
			wData[batch * i + j] = dcw[j].real();
		}
		delete[] dcw;
	}
	complex<double>* dcw = scheme.decrypt(secretKey, cwData[cnum-1]);
	long rest = factorDim - batch * (cnum - 1);
	for (long j = 0; j < rest; ++j) {
		wData[batch * (cnum - 1) + j] = dcw[j].real();
	}
	delete[] dcw;
}
