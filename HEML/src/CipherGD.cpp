#include "CipherGD.h"

#include <Cipher.h>
#include <CZZ.h>
#include <EvaluatorUtils.h>
#include <NTL/BasicThreadPool.h>
#include <NTL/RR.h>
#include <NTL/ZZ.h>
#include <SchemeAux.h>
#include <cmath>

void CipherGD::encxyData(Cipher*& cxyData, long**& xyData, long& slots, long& factorDim, long& learnDim, long& xyBatch, long& cnum, long& xyBits) {
	ZZ precision = power2_ZZ(xyBits);
	NTL_EXEC_RANGE(cnum, first, last);
	for (long i = first; i < last; ++i) {
		CZZ* pxyData = new CZZ[slots];
		if(i != (cnum - 1)) {
			for (long l = 0; l < xyBatch; ++l) {
				for (long j = 0; j < learnDim; ++j) {
					if(xyData[j][xyBatch * i + l] == -1) {
						pxyData[xyBatch * j + l] = CZZ(-precision);
					} else if(xyData[j][xyBatch * i + l] == 1) {
						pxyData[xyBatch * j + l] = CZZ(precision);
					}
				}
			}
		} else {
			long rest = factorDim - xyBatch * i;
			for (long l = 0; l < rest; ++l) {
				for (long j = 0; j < learnDim; ++j) {
					if(xyData[j][xyBatch * i + l] == -1) {
						pxyData[xyBatch * j + l] = CZZ(-precision);
					} else if(xyData[j][xyBatch * i + l] == 1) {
						pxyData[xyBatch * j + l] = CZZ(precision);
					}
				}
			}
		}
		cxyData[i] = scheme.encrypt(pxyData, slots);
		delete[] pxyData;
	}
	NTL_EXEC_INDEX_END;
}

void CipherGD::encwData(Cipher*& cwData, Cipher*& cxyData, long& slotBits, long& ldimBits, long& xyBatchBits, long& cnum, long& xyBits, long& wBits) {
	long downBits = ldimBits + xyBits - wBits;
	NTL_EXEC_RANGE(cnum, first, last);
	for (long i = first; i < last; ++i) {
		cwData[i] = cxyData[i];
		for (long l = xyBatchBits; l < slotBits; ++l) {
			Cipher rot = scheme.leftRotateByPo2(cwData[i], l);
			scheme.addAndEqual(cwData[i], rot);
		}
		scheme.modSwitchAndEqual(cwData[i], downBits);
	}
	NTL_EXEC_RANGE_END;
}

void CipherGD::encStepNLGD5(Cipher*& cxyData, Cipher*& cwData, Cipher*& cvData, ZZX& poly, long& slots, long& learnDim, long& xybatchBits, long& xyBatch, long& cnum, double& gamma, double& eta, double& etaprev, long& xyBits, long& wBits, long& pBits) {
	Cipher* cprod = new Cipher[cnum];

	long bitsDown = cxyData[0].cbits - cwData[0].cbits;
	NTL_EXEC_RANGE(cnum, first, last);
	for (long i = first; i < last; ++i) {
		scheme.modEmbedAndEqual(cxyData[i], bitsDown);
		cprod[i] = scheme.mult(cxyData[i], cwData[i]);
		for (long l = 0; l < xybatchBits; ++l) {
			Cipher rot = scheme.leftRotateByPo2(cprod[i], l);
			scheme.addAndEqual(cprod[i], rot);
		}
	}
	NTL_EXEC_RANGE_END;

	Cipher cip = cprod[0];
	for (long i = 1; i < cnum; ++i) {
		scheme.addAndEqual(cip, cprod[i]);
	}
	scheme.modSwitchAndEqual(cip, xyBits);

	scheme.multByPolyAndEqual(cip, poly);
	for (long l = 0; l < xybatchBits; ++l) {
		Cipher rot = scheme.rightRotateByPo2(cip, l);
		scheme.addAndEqual(cip, rot);
	}
	scheme.modSwitchAndEqual(cip, pBits);

 	double* coeffs = scheme.aux.taylorCoeffsMap.at(SIGMOIDPRIMEGOOD5);

 	Cipher cip2 = scheme.square(cip);
 	scheme.modSwitchAndEqual(cip2, wBits);

 	ZZ c0zz = EvaluatorUtils::evaluateVal(coeffs[0], wBits);
 	ZZ c1zz = EvaluatorUtils::evaluateVal(coeffs[1], wBits);
 	double c3c5 = coeffs[3] / coeffs[5];
 	ZZ c3c5zz = EvaluatorUtils::evaluateVal(c3c5, wBits);
 	ZZ c5zz = EvaluatorUtils::evaluateVal(coeffs[5], wBits);

 	Cipher cipc0c1 = scheme.multByConst(cip, c1zz);
 	scheme.modSwitchAndEqual(cipc0c1, wBits);
 	scheme.addConstAndEqual(cipc0c1, c0zz);

 	Cipher cip2c3c5 = scheme.addConst(cip2, c3c5zz);
 	Cipher cipc5 = scheme.multByConst(cip, c5zz);
 	scheme.modSwitchAndEqual(cipc5, wBits);
 	scheme.multAndEqual(cip2c3c5, cipc5);
 	scheme.modSwitchAndEqual(cip2c3c5, wBits);

 	Cipher* cgrad = new Cipher[cnum];

	RR tmprr = to_RR(gamma) * (1. - eta);
	ZZ tmpzz = EvaluatorUtils::evaluateVal(tmprr, wBits);

	NTL_EXEC_RANGE(cnum, first, last);
	for (long i = first; i < last; ++i) {
		cprod[i] = scheme.multByConst(cxyData[i], tmpzz);
		scheme.modSwitchAndEqual(cprod[i], xyBits);
		scheme.modEmbed(cprod[i], pBits  + wBits);
		cgrad[i] = scheme.mult(cprod[i], cip2);
		scheme.modSwitchAndEqual(cgrad[i], wBits);
		scheme.multAndEqual(cgrad[i], cip2c3c5);
		scheme.modSwitchAndEqual(cgrad[i], wBits);
		scheme.multAndEqual(cprod[i], cipc0c1);
		scheme.modSwitchAndEqual(cprod[i], wBits);
		scheme.modEmbedAndEqual(cprod[i], wBits);
		scheme.addAndEqual(cgrad[i], cprod[i]);
	}
	NTL_EXEC_RANGE_END;

	long logslots = log2(slots);

	NTL_EXEC_RANGE(cnum, first, last);
	for (long i = first; i < last; ++i) {
		for (long l = xybatchBits; l < logslots; ++l) {
			Cipher rot = scheme.leftRotateByPo2(cgrad[i], l);
			scheme.addAndEqual(cgrad[i], rot);
		}
	}
	NTL_EXEC_RANGE_END;

	RR cnst = to_RR(1. - eta);
	ZZ pcnst = EvaluatorUtils::evaluateVal(cnst, wBits);
	bitsDown = xyBits + pBits + 2 * wBits;
	if(etaprev < 0) {
		RR cnst2 = to_RR(eta) / (1 - etaprev);
		ZZ pcnst2 = EvaluatorUtils::evaluateVal(cnst2, wBits);
		NTL_EXEC_RANGE(cnum, first, last);
		for (long i = first; i < last; ++i) {
			scheme.multByConstAndEqual(cwData[i], pcnst);
			scheme.modSwitchAndEqual(cwData[i], wBits);
			scheme.modEmbedAndEqual(cwData[i], bitsDown);
			cgrad[i] = scheme.sub(cwData[i], cgrad[i]);

			scheme.multByConstAndEqual(cvData[i], pcnst2);
			scheme.modSwitchAndEqual(cvData[i], wBits);
			scheme.modEmbedAndEqual(cvData[i], bitsDown);

			cwData[i] = scheme.add(cvData[i], cgrad[i]);
			cvData[i] = cgrad[i];
		}
		NTL_EXEC_RANGE_END;
	} else {
		NTL_EXEC_RANGE(cnum, first, last);
		for (long i = first; i < last; ++i) {
			cvData[i] = scheme.multByConst(cwData[i], pcnst);
			scheme.modSwitchAndEqual(cvData[i], wBits);
			scheme.modEmbedAndEqual(cvData[i], bitsDown);
			scheme.subAndEqual(cvData[i], cgrad[i]);
			scheme.modEmbedAndEqual(cwData[i], bitsDown + wBits);
			scheme.subAndEqual(cwData[i], cgrad[i]);
		}
		NTL_EXEC_RANGE_END;
	}
	delete[] cprod;
	delete[] cgrad;
}

void CipherGD::encStepNLGD3(Cipher*& cxyData, Cipher*& cwData, Cipher*& cvData, ZZX& poly, long& slots, long& learnDim, long& xybatchBits, long& xyBatch, long& cnum, double& gamma, double& eta, double& etaprev, long& xyBits, long& wBits, long& pBits) {

	Cipher* cprod = new Cipher[cnum];
	long bitsDown = cxyData[0].cbits - cwData[0].cbits;

	NTL_EXEC_RANGE(cnum, first, last);
	for (long i = first; i < last; ++i) {
		scheme.modEmbedAndEqual(cxyData[i], bitsDown);
		cprod[i] = scheme.mult(cxyData[i], cwData[i]);
		for (long l = 0; l < xybatchBits; ++l) {
			Cipher rot = scheme.leftRotateByPo2(cprod[i], l);
			scheme.addAndEqual(cprod[i], rot);
		}
	}
	NTL_EXEC_RANGE_END;
	Cipher cip = cprod[0];
	for (long i = 1; i < cnum; ++i) {
		scheme.addAndEqual(cip, cprod[i]);
	}
	scheme.modSwitchAndEqual(cip, xyBits);

	scheme.multByPolyAndEqual(cip, poly);
	for (long l = 0; l < xybatchBits; ++l) {
		Cipher rot = scheme.rightRotateByPo2(cip, l);
		scheme.addAndEqual(cip, rot);
	}
	scheme.modSwitchAndEqual(cip, pBits);

 	double* coeffs = scheme.aux.taylorCoeffsMap.at(SIGMOIDPRIMEGOOD5);

 	Cipher cip2 = scheme.square(cip);
 	scheme.modSwitchAndEqual(cip2, wBits);

 	ZZ c0zz = EvaluatorUtils::evaluateVal(coeffs[0], wBits);
 	ZZ c1zz = EvaluatorUtils::evaluateVal(coeffs[1], wBits);

 	Cipher cipc0c1 = scheme.multByConst(cip, c1zz);
 	scheme.modSwitchAndEqual(cipc0c1, wBits);
 	scheme.addConstAndEqual(cipc0c1, c0zz);

 	Cipher* cgrad = new Cipher[cnum];

	RR tmprr = to_RR(gamma) * (1. - eta);
	ZZ tmpzz = EvaluatorUtils::evaluateVal(tmprr, wBits);

	RR tmpc3rr = tmprr * coeffs[3];
	ZZ tmpc3zz = EvaluatorUtils::evaluateVal(tmpc3rr, wBits);

	NTL_EXEC_RANGE(cnum, first, last);
	for (long i = first; i < last; ++i) {
		cprod[i] = scheme.multByConst(cxyData[i], tmpc3zz);
		scheme.modSwitchAndEqual(cprod[i], xyBits);
		scheme.modEmbed(cprod[i], pBits);
		scheme.multAndEqual(cprod[i], cip);
		scheme.modSwitchAndEqual(cprod[i], wBits);
		scheme.multAndEqual(cprod[i], cip2);

		cgrad[i] = scheme.multByConst(cxyData[i], tmpzz);
		scheme.modSwitchAndEqual(cgrad[i], xyBits);
		scheme.modEmbed(cgrad[i], pBits + wBits);
		scheme.multAndEqual(cgrad[i], cipc0c1);
		scheme.addAndEqual(cgrad[i], cprod[i]);
		scheme.modSwitchAndEqual(cgrad[i], wBits);
	}
	NTL_EXEC_RANGE_END;

	long logslots = log2(slots);

	NTL_EXEC_RANGE(cnum, first, last);
	for (long i = first; i < last; ++i) {
		for (long l = xybatchBits; l < logslots; ++l) {
			Cipher rot = scheme.leftRotateByPo2(cgrad[i], l);
			scheme.addAndEqual(cgrad[i], rot);
		}
	}
	NTL_EXEC_RANGE_END;

	RR cnst = to_RR(1. - eta);
	ZZ pcnst = EvaluatorUtils::evaluateVal(cnst, wBits);
	bitsDown = xyBits + pBits + 1 * wBits;
	if(etaprev < 0) {
		RR cnst2 = to_RR(eta) / (1 - etaprev);
		ZZ pcnst2 = EvaluatorUtils::evaluateVal(cnst2, wBits);

		NTL_EXEC_RANGE(cnum, first, last);
		for (long i = first; i < last; ++i) {
			scheme.multByConstAndEqual(cwData[i], pcnst);
			scheme.modSwitchAndEqual(cwData[i], wBits);
			scheme.modEmbedAndEqual(cwData[i], bitsDown);
			cgrad[i] = scheme.sub(cwData[i], cgrad[i]);

			scheme.multByConstAndEqual(cvData[i], pcnst2);
			scheme.modSwitchAndEqual(cvData[i], wBits);
			scheme.modEmbedAndEqual(cvData[i], bitsDown);

			cwData[i] = scheme.add(cvData[i], cgrad[i]);
			cvData[i] = cgrad[i];
		}
		NTL_EXEC_RANGE_END;
	} else {
		NTL_EXEC_RANGE(cnum, first, last);
		for (long i = first; i < last; ++i) {
			cvData[i] = scheme.multByConst(cwData[i], pcnst);
			scheme.modSwitchAndEqual(cvData[i], wBits);
			scheme.modEmbedAndEqual(cvData[i], bitsDown);
			scheme.subAndEqual(cvData[i], cgrad[i]);
			scheme.modEmbedAndEqual(cwData[i], bitsDown);
			scheme.subAndEqual(cwData[i], cgrad[i]);
		}
		NTL_EXEC_RANGE_END;
	}
	delete[] cprod;
	delete[] cgrad;
}

void CipherGD::decwData(double*& w, Cipher*& cw, long& factorDim, long& xyBatch, long& cnum, long& wBits) {
	for (long i = 0; i < (cnum - 1); ++i) {
		CZZ* dcw = scheme.decrypt(secretKey, cw[i]);
		for (long l = 0; l < xyBatch; ++l) {
			RR wi = to_RR(dcw[l].r);
			wi.e -= wBits;
			w[xyBatch * i + l] = to_double(wi);
		}
		delete[] dcw;
	}
	CZZ* dcw = scheme.decrypt(secretKey, cw[cnum-1]);
	long rest = factorDim - xyBatch * (cnum - 1);
	for (long l = 0; l < rest; ++l) {
		RR wi = to_RR(dcw[l].r);
		wi.e -= wBits;
		w[xyBatch * (cnum - 1) + l] = to_double(wi);
	}
	delete[] dcw;
}
