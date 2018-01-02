#include "CipherGD.h"

#include <Ciphertext.h>
#include <CZZ.h>
#include <EvaluatorUtils.h>
#include <NTL/BasicThreadPool.h>
#include <NTL/RR.h>
#include <NTL/ZZ.h>

void CipherGD::encxyData(Ciphertext* cxyData, long** xyData, long& slots, long& factorDim, long& learnDim, long& batch, long& cnum, long& xyBits) {
	ZZ precision = power2_ZZ(xyBits);
	CZZ* pxyData = new CZZ[slots];
	for (long i = 0; i < cnum - 1; ++i) {
		for (long j = 0; j < learnDim; ++j) {
			for (long l = 0; l < batch; ++l) {
				pxyData[batch * j + l] = xyData[j][batch * i + l] == -1 ? CZZ(-precision) :
						xyData[j][batch * i + l] == 1 ? CZZ(precision) : CZZ();
			}
		}
		cxyData[i] = scheme.encrypt(pxyData, slots, scheme.context.logQ);
	}

	long rest = factorDim - batch * (cnum - 1);
	for (long j = 0; j < learnDim; ++j) {
		for (long l = 0; l < rest; ++l) {
			pxyData[batch * j + l] = xyData[j][batch * (cnum - 1) + l] == -1 ? CZZ(-precision) :
					xyData[j][batch * (cnum - 1) + l] == 1 ? CZZ(precision) : CZZ();
		}
		for (long l = rest; l < batch; ++l) {
			pxyData[batch * j + l] = CZZ();
		}
	}
	cxyData[cnum - 1] = scheme.encrypt(pxyData, slots, scheme.context.logQ);

	delete[] pxyData;
}

void CipherGD::encwData(Ciphertext* cwData, Ciphertext* cxyData, long& cnum, long& sBits, long& ldimBits, long& bBits, long& xyBits, long& wBits) {
	long downBits = ldimBits + xyBits - wBits;
	NTL_EXEC_RANGE(cnum, first, last);
	for (long i = first; i < last; ++i) {
		cwData[i] = cxyData[i];
		for (long l = bBits; l < sBits; ++l) {
			Ciphertext rot = scheme.leftRotateByPo2(cwData[i], l);
			scheme.addAndEqual(cwData[i], rot);
		}
		scheme.reScaleByAndEqual(cwData[i], downBits);
	}
	NTL_EXEC_RANGE_END;
}

ZZX CipherGD::generateAuxPoly(long& slots, long& batch, long& pBits) {
	ZZ p = power2_ZZ(pBits);
	CZZ* pvals = new CZZ[slots];
	for (long j = 0; j < slots; j += batch) {
		pvals[j] = CZZ(p);
	}
	Plaintext msg = scheme.encode(pvals, slots, scheme.context.logQ);
	delete[] pvals;
	msg.mx >>= scheme.context.logQ;
	return msg.mx;
}

Ciphertext CipherGD::encIP(Ciphertext* cxyData, Ciphertext* cwData, Ciphertext* cgrad, Ciphertext* cprod, ZZX& poly, long& cnum, long& bBits, long& xyBits, long& pBits, long& aBits) {
	long bitsDown = cxyData[0].logq - cwData[0].logq;

	NTL_EXEC_RANGE(cnum, first, last);
	for (long i = first; i < last; ++i) {
		scheme.modDownByAndEqual(cxyData[i], bitsDown);
		cprod[i] = scheme.mult(cxyData[i], cwData[i]);
		for (long l = 0; l < bBits; ++l) {
			cgrad[i] = scheme.leftRotateByPo2(cprod[i], l);
			scheme.addAndEqual(cprod[i], cgrad[i]);
		}
	}
	NTL_EXEC_RANGE_END

	Ciphertext cip = cprod[0];
	for (long i = 1; i < cnum; ++i) {
		scheme.addAndEqual(cip, cprod[i]);
	}
	scheme.reScaleByAndEqual(cip, xyBits);
	scheme.multByPolyAndEqual(cip, poly);
	for (long l = 0; l < bBits; ++l) {
		cprod[0] = scheme.rightRotateByPo2(cip, l);
		scheme.addAndEqual(cip, cprod[0]);
	}
	scheme.reScaleByAndEqual(cip, pBits + aBits); // (xyBits + pBits + aBits, wBits - aBits)
	return cip;
}

void CipherGD::encSigmoid(long& approxDeg, Ciphertext* cxyData, Ciphertext* cgrad, Ciphertext* cprod, Ciphertext& cip, long& cnum, double& gamma, long& sBits, long & bBits, long& xyBits, long& wBits, long& pBits, long& aBits) {
	Ciphertext cip2 = scheme.square(cip);
	scheme.reScaleByAndEqual(cip2, wBits); // cip2 (xyBits + pBits + aBits + wBits, wBits - 2aBits)

	if(approxDeg == 3) {
		RR tmp = to_RR(degree3[1]) / degree3[2];
		ZZ tmpzz = EvaluatorUtils::evalZZ(tmp, wBits - 2 * aBits);
		scheme.addConstAndEqual(cip2, tmpzz); // cip2 (xyBits + pBits + aBits + wBits, wBits - 2aBits)
		tmp = to_RR(gamma) * degree3[0];
		tmpzz = EvaluatorUtils::evalZZ(tmp, wBits);
		tmp = to_RR(gamma) * degree3[2];
		ZZ tmpzz1 = EvaluatorUtils::evalZZ(tmp, wBits + 3 * aBits);

		NTL_EXEC_RANGE(cnum, first, last);
		for (long i = first; i < last; ++i) {
			cprod[i] = scheme.multByConst(cxyData[i], tmpzz); // cprod (0, xyBits + wBits)
			scheme.reScaleByAndEqual(cprod[i], xyBits);// cprod (xyData, wBits)
			scheme.modDownByAndEqual(cprod[i], pBits + aBits + 2 * wBits);// cprod (xyData + pBits + aBits + 2wBits, wBits)

			cgrad[i] = scheme.multByConst(cxyData[i], tmpzz1);// cgrad (0, xyData + wBits + 3aBits)
			scheme.reScaleByAndEqual(cgrad[i], xyBits);// cgrad (xyData, wBits + 3aBits)
			scheme.modDownByAndEqual(cgrad[i], pBits + aBits);// cgrad (xyData + pBits + aBits, wBits + 3aBits)
			scheme.multAndEqual(cgrad[i], cip);// cgrad (xyData + pBits + aBits, 2wBits + 2aBits)
			scheme.reScaleByAndEqual(cgrad[i], wBits);// cgrad (xyData + pBits + aBits + wBits, wBits + 2aBits)
			scheme.multAndEqual(cgrad[i], cip2);// cgrad (xyData + pBits + aBits + 2wBits, 2wBits)
			scheme.reScaleByAndEqual(cgrad[i], wBits);// cgrad (xyData + pBits + aBits + 2wBits, wBits)
			scheme.addAndEqual(cgrad[i], cprod[i]);// cgrad (xyData + pBits + aBits + 2wBits, wBits)
		}
		NTL_EXEC_RANGE_END;
	} else if (approxDeg == 5) {
		Ciphertext cip4 = scheme.square(cip2);
		RR tmp = to_RR(degree5[2]) / degree5[3];
		ZZ tmpzz = EvaluatorUtils::evalZZ(tmp, wBits - 2 * aBits);
		scheme.multByConstAndEqual(cip2, tmpzz); // cip2 (xyBits + pBits + aBits + wBits, 2wBits - 4aBits)
		scheme.reScaleByAndEqual(cip2, wBits); // cip2 (xyBits + pBits + aBits + 2wBits, wBits - 4aBits)
		scheme.reScaleByAndEqual(cip4, wBits); // cip4 (xyBits + pBits + aBits + 2wBits, wBits - 4aBits)
		scheme.addAndEqual(cip4, cip2);
		tmp = to_RR(degree5[1]) / degree5[3];
		tmpzz = EvaluatorUtils::evalZZ(tmp, wBits - 4 * aBits);
		scheme.addConstAndEqual(cip4, tmpzz); // cip4 (xyBits + pBits + aBits + 2wBits, wBits - 4aBits)
		tmp = to_RR(gamma) * degree5[0];
		tmpzz = EvaluatorUtils::evalZZ(tmp, wBits);
		tmp = to_RR(gamma) * degree5[3];
		ZZ tmpzz1 = EvaluatorUtils::evalZZ(tmp, wBits + 5 * aBits);

		NTL_EXEC_RANGE(cnum, first, last);
		for (long i = first; i < last; ++i) {
			cprod[i] = scheme.multByConst(cxyData[i], tmpzz); // cprod (0, xyBits + wBits)
			scheme.reScaleByAndEqual(cprod[i], xyBits);// cprod (xyData, wBits)
			scheme.modDownByAndEqual(cprod[i], pBits + aBits + 3 * wBits);// cprod (xyData + pBits + aBits + 3wBits, wBits)

			cgrad[i] = scheme.multByConst(cxyData[i], tmpzz1);// cgrad (0, xyData + wBits + 5aBits)
			scheme.reScaleByAndEqual(cgrad[i], xyBits);// cgrad (xyData, wBits + 5aBits)
			scheme.modDownByAndEqual(cgrad[i], pBits + aBits);// cgrad (xyData + pBits + aBits, wBits + 5aBits)
			scheme.multAndEqual(cgrad[i], cip);// cgrad (xyData + pBits + aBits, 2wBits + 4aBits)
			scheme.reScaleByAndEqual(cgrad[i], wBits);// cgrad (xyData + pBits + aBits + wBits, wBits + 4aBits)
			scheme.modDownByAndEqual(cgrad[i], wBits);// cgrad (xyData + pBits + aBits + 2wBits, wBits + 4aBits)
			scheme.multAndEqual(cgrad[i], cip4);// cgrad (xyData + pBits + aBits + 2wBits, 2wBits)
			scheme.reScaleByAndEqual(cgrad[i], wBits);// cgrad (xyData + pBits + aBits + 3wBits, wBits)
			scheme.addAndEqual(cgrad[i], cprod[i]);// cgrad (xyData + pBits + aBits + 3wBits, wBits)
		}
		NTL_EXEC_RANGE_END;
	} else {
		Ciphertext cip4 = scheme.square(cip2);
		scheme.reScaleByAndEqual(cip4, wBits); // cip4 (xyBits + pBits + aBits + 2wBits, wBits - 4aBits)
		RR tmp = to_RR(degree7[3]) / degree7[4];
		ZZ tmpzz = EvaluatorUtils::evalZZ(tmp, wBits - 2 * aBits);
		Ciphertext cip2c = scheme.multByConst(cip2, tmpzz); // cip2c (xyBits + pBits + aBits + wBits, 2wBits - 4aBits)
		scheme.reScaleByAndEqual(cip2c, wBits); // cip2c (xyBits + pBits + aBits + 2wBits, wBits - 4aBits)
		scheme.addAndEqual(cip4, cip2c);
		tmp = to_RR(degree7[2]) / degree7[4];
		tmpzz = EvaluatorUtils::evalZZ(tmp, wBits - 4 * aBits);
		scheme.addConstAndEqual(cip4, tmpzz); // cip4 (xyBits + pBits + aBits + 2wBits, wBits - 4aBits)
		tmp = to_RR(gamma) * degree7[0];
		tmpzz = EvaluatorUtils::evalZZ(tmp, wBits);
		tmp = to_RR(gamma) * degree7[1];
		ZZ tmpzz1 = EvaluatorUtils::evalZZ(tmp, wBits + aBits);
		tmp = to_RR(gamma) * degree7[4];
		ZZ tmpzz2 = EvaluatorUtils::evalZZ(tmp, wBits + 7 * aBits);
		NTL_EXEC_RANGE(cnum, first, last);
		for (long i = first; i < last; ++i) {
			cgrad[i] = scheme.multByConst(cxyData[i], tmpzz); // cgrad (0, xyBits + wBits)
			scheme.reScaleByAndEqual(cgrad[i], xyBits);// cgrad (xybits, wBits)
			scheme.modDownByAndEqual(cgrad[i], pBits + aBits + 3 * wBits);// cgrad (xyBits + pBits + aBits + 3wBits, wBits)

			cprod[i] = scheme.multByConst(cxyData[i], tmpzz1);// cprod (0, xyBits + wBits + aBits);
			scheme.reScaleByAndEqual(cprod[i], xyBits);// cprod (xyBits, wBits + aBits);
			scheme.modDownByAndEqual(cprod[i], pBits + aBits);// cprod (xyBits + pBits + aBits, wBits + aBits);
			scheme.multAndEqual(cprod[i], cip);// cprod (xyBits + pBits + aBits, 2wBits);
			scheme.reScaleByAndEqual(cprod[i], wBits);// cprod (xyBits + pBits + aBits + wBits, wBits);
			scheme.modDownByAndEqual(cprod[i], 2 * wBits);// cprod (xyBits + pBits + aBits + 3wBits, wBits);
			scheme.addAndEqual(cprod[i], cgrad[i]);// cprod (xyBits + pBits + aBits + 3wBits, wBits);

			cgrad[i] = scheme.multByConst(cxyData[i], tmpzz2);// cgrad (0, xyData + wBits + 7aBits)
			scheme.reScaleByAndEqual(cgrad[i], xyBits);// cgrad (xyBits, wBits + 7aBits)
			scheme.modDownByAndEqual(cgrad[i], pBits + aBits);// cgrad (xyBits + pBits + aBits, wBits + 7aBits)
			scheme.multAndEqual(cgrad[i], cip);// cgrad (xyData + pBits + aBits, 2wBits + 6aBits)
			scheme.reScaleByAndEqual(cgrad[i], wBits);// cgrad (xyBits + pBits + aBits + wBits, wBits + 6aBits)
			scheme.multAndEqual(cgrad[i], cip2);// cgrad (xyBits + pBits + aBits + wBits, 2wBits + 4aBits)
			scheme.reScaleByAndEqual(cgrad[i], wBits);// cgrad (xyData + pBits + aBits + 2wBits, wBits + 4aBits)
			scheme.multAndEqual(cgrad[i], cip4);// cgrad (xyData + pBits + aBits + 2wBits, 2wBits)
			scheme.reScaleByAndEqual(cgrad[i], wBits);// cgrad (xyData + pBits + aBits + 3wBits, wBits)
			scheme.addAndEqual(cgrad[i], cprod[i]);// cgrad (xyData + pBits + aBits + 3wBits, wBits)
		}
		NTL_EXEC_RANGE_END;

	}

	NTL_EXEC_RANGE(cnum, first, last);
	for (long i = first; i < last; ++i) {
		for (long l = bBits; l < sBits; ++l) {
			cprod[i] = scheme.leftRotateByPo2(cgrad[i], l);
			scheme.addAndEqual(cgrad[i], cprod[i]);
		}
	}
	NTL_EXEC_RANGE_END;
}

void CipherGD::encLGDstep(Ciphertext* cwData, Ciphertext* cgrad, long& cnum, long& wBits, long& bitsDown) {
	NTL_EXEC_RANGE(cnum, first, last);
	for (long i = first; i < last; ++i) {
		scheme.modDownByAndEqual(cwData[i], bitsDown + wBits);
		scheme.subAndEqual(cwData[i], cgrad[i]);
	}
	NTL_EXEC_RANGE_END;
}

void CipherGD::encMLGDstep(Ciphertext* cwData, Ciphertext* cvData, Ciphertext* cgrad, double& eta, long& cnum, long& wBits, long& bitsDown) {
	RR tmp = to_RR(eta);
	ZZ tmpzz = EvaluatorUtils::evalZZ(tmp, wBits);
	NTL_EXEC_RANGE(cnum, first, last);
	for (long i = first; i < last; ++i) {
		scheme.multByConstAndEqual(cvData[i], tmpzz);
		scheme.reScaleByAndEqual(cvData[i], wBits);
		scheme.modDownByAndEqual(cvData[i], bitsDown);
		scheme.addAndEqual(cvData[i], cgrad[i]);
		scheme.modDownByAndEqual(cwData[i], bitsDown + wBits);
		scheme.subAndEqual(cwData[i],cvData[i]);
	}
	NTL_EXEC_RANGE_END;
}

void CipherGD::encNLGDstep(Ciphertext* cwData, Ciphertext* cvData, Ciphertext* cgrad, double& eta, double& etaprev, long& cnum, long& wBits, long& bitsDown) {
	RR tmp = to_RR(1. - eta);
	ZZ tmpzz = EvaluatorUtils::evalZZ(tmp, wBits);
	tmp = to_RR(eta / (1 - etaprev));
	ZZ tmpzz1 = EvaluatorUtils::evalZZ(tmp, wBits);

	NTL_EXEC_RANGE(cnum, first, last);
	for (long i = first; i < last; ++i) {
		scheme.multByConstAndEqual(cvData[i], tmpzz);
		scheme.reScaleByAndEqual(cvData[i], wBits);
		scheme.modDownByAndEqual(cvData[i], bitsDown);
		cgrad[i] = scheme.sub(cvData[i], cgrad[i]);

		scheme.multByConstAndEqual(cwData[i], tmpzz1);
		scheme.reScaleByAndEqual(cwData[i], wBits);
		scheme.modDownByAndEqual(cwData[i], bitsDown);

		cvData[i] = scheme.add(cwData[i], cgrad[i]);
		cwData[i] = cgrad[i];
	}
	NTL_EXEC_RANGE_END;
}

void CipherGD::encLGDiteration(long& approxDeg, Ciphertext* cxyData, Ciphertext* cwData, ZZX& poly, long& cnum, double& gamma, long& sBits, long& bBits, long& xyBits, long& wBits, long& pBits, long& aBits) {
	Ciphertext* cprod = new Ciphertext[cnum];
 	Ciphertext* cgrad = new Ciphertext[cnum];

	Ciphertext cip = encIP(cxyData, cwData, cgrad, cprod, poly, cnum, bBits, xyBits, pBits, aBits);
	encSigmoid(approxDeg, cxyData, cgrad, cprod, cip, cnum, gamma, sBits, bBits, xyBits, wBits, pBits, aBits);
	delete[] cprod;

	long bitsDown = approxDeg == 3 ? xyBits + pBits + aBits + wBits : xyBits + pBits + aBits + 2 * wBits;

	encLGDstep(cwData, cgrad, cnum, wBits, bitsDown);

	delete[] cgrad;
}

void CipherGD::encMLGDiteration(long& approxDeg, Ciphertext* cxyData, Ciphertext* cwData, Ciphertext* cvData, ZZX& poly, long& cnum, double& gamma, double& eta, long& sBits, long& bBits, long& xyBits, long& wBits, long& pBits, long& aBits) {
	Ciphertext* cprod = new Ciphertext[cnum];
 	Ciphertext* cgrad = new Ciphertext[cnum];

	Ciphertext cip = encIP(cxyData, cwData, cgrad, cprod, poly, cnum, bBits, xyBits, pBits, aBits);
	encSigmoid(approxDeg, cxyData, cgrad, cprod, cip, cnum, gamma, sBits, bBits, xyBits, wBits, pBits, aBits);

	delete[] cprod;

	long bitsDown = approxDeg == 3 ? xyBits + pBits + aBits + wBits : xyBits + pBits + aBits + 2 * wBits;

	encMLGDstep(cwData, cvData, cgrad, eta, cnum, wBits, bitsDown);

	delete[] cgrad;
}

void CipherGD::encNLGDiteration(long& approxDeg, Ciphertext* cxyData, Ciphertext* cwData, Ciphertext* cvData, ZZX& poly, long& cnum, double& gamma, double& eta, double& etaprev, long& sBits, long& bBits, long& xyBits, long& wBits, long& pBits, long& aBits) {
	Ciphertext* cprod = new Ciphertext[cnum];
 	Ciphertext* cgrad = new Ciphertext[cnum];

	Ciphertext cip = encIP(cxyData, cvData, cgrad, cprod, poly, cnum, bBits, xyBits, pBits, aBits);

	encSigmoid(approxDeg, cxyData, cgrad, cprod, cip, cnum, gamma, sBits, bBits, xyBits, wBits, pBits, aBits);

	delete[] cprod;

	long bitsDown = approxDeg == 3 ? xyBits + pBits + aBits + wBits : xyBits + pBits + aBits + 2 * wBits;

	encNLGDstep(cwData, cvData, cgrad, eta, etaprev, cnum, wBits, bitsDown);

	delete[] cgrad;
}

void CipherGD::decwData(double* wData, Ciphertext* cwData, long& factorDim, long& batch, long& cnum, long& wBits) {
	for (long i = 0; i < (cnum - 1); ++i) {
		CZZ* dcw = scheme.decrypt(secretKey, cwData[i]);
//		EvaluatorUtils::evalRealArray(wData + batch * i, dcw, batch, wBits);
		for (long l = 0; l < batch; ++l) {
			RR wi = to_RR(dcw[l].r);
			wi.e -= wBits;
			wData[batch * i + l] = to_double(wi);
		}
		delete[] dcw;
	}
	CZZ* dcw = scheme.decrypt(secretKey, cwData[cnum-1]);
	long rest = factorDim - batch * (cnum - 1);
//	EvaluatorUtils::evalRealArray(wData + batch * (cnum - 1), dcw, rest, wBits);
	for (long l = 0; l < rest; ++l) {
		RR wi = to_RR(dcw[l].r);
		wi.e -= wBits;
		wData[batch * (cnum - 1) + l] = to_double(wi);
	}
	delete[] dcw;
}
