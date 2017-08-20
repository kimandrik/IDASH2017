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

void CipherGD::encStepNLGD5(Cipher*& cxyData, Cipher*& cwData, Cipher*& cvData, ZZX& poly, long& slotBits, long& xybatchBits, long& cnum, double& gamma, double& eta, double& etaprev, long& xyBits, long& wBits, long& pBits) {
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

 	ZZ tmpzz = EvaluatorUtils::evaluateVal(coeffs[1], wBits);
 	Cipher cipc0c1 = scheme.multByConst(cip, tmpzz);
 	scheme.modSwitchAndEqual(cipc0c1, wBits);

 	tmpzz = EvaluatorUtils::evaluateVal(coeffs[0], wBits);
 	scheme.addConstAndEqual(cipc0c1, tmpzz);

 	double c3c5 = coeffs[3] / coeffs[5];
 	tmpzz = EvaluatorUtils::evaluateVal(c3c5, wBits);
 	Cipher cip2c3c5 = scheme.addConst(cip2, tmpzz);

 	tmpzz = EvaluatorUtils::evaluateVal(coeffs[5], wBits);
 	Cipher cipc5 = scheme.multByConst(cip, tmpzz);
 	scheme.modSwitchAndEqual(cipc5, wBits);
 	scheme.multAndEqual(cip2c3c5, cipc5);
 	scheme.modSwitchAndEqual(cip2c3c5, wBits);

 	Cipher* cgrad = new Cipher[cnum];

	RR tmprr = to_RR(gamma) * (1. - eta);
	tmpzz = EvaluatorUtils::evaluateVal(tmprr, wBits);

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

	NTL_EXEC_RANGE(cnum, first, last);
	for (long i = first; i < last; ++i) {
		for (long l = xybatchBits; l < slotBits; ++l) {
			Cipher rot = scheme.leftRotateByPo2(cgrad[i], l);
			scheme.addAndEqual(cgrad[i], rot);
		}
	}
	NTL_EXEC_RANGE_END;

	tmprr = to_RR(1. - eta);
	tmpzz = EvaluatorUtils::evaluateVal(tmprr, wBits);
	bitsDown = xyBits + pBits + 2 * wBits;
	if(etaprev < 0) {
		tmprr = to_RR(eta) / (1 - etaprev);
		ZZ tmpzz2 = EvaluatorUtils::evaluateVal(tmprr, wBits);
		NTL_EXEC_RANGE(cnum, first, last);
		for (long i = first; i < last; ++i) {
			scheme.multByConstAndEqual(cwData[i], tmpzz);
			scheme.modSwitchAndEqual(cwData[i], wBits);
			scheme.modEmbedAndEqual(cwData[i], bitsDown);
			cgrad[i] = scheme.sub(cwData[i], cgrad[i]);

			scheme.multByConstAndEqual(cvData[i], tmpzz2);
			scheme.modSwitchAndEqual(cvData[i], wBits);
			scheme.modEmbedAndEqual(cvData[i], bitsDown);

			cwData[i] = scheme.add(cvData[i], cgrad[i]);
			cvData[i] = cgrad[i];
		}
		NTL_EXEC_RANGE_END;
	} else {
		NTL_EXEC_RANGE(cnum, first, last);
		for (long i = first; i < last; ++i) {
			cvData[i] = scheme.multByConst(cwData[i], tmpzz);
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

void CipherGD::encStepNLGD3(Cipher*& cxyData, Cipher*& cwData, Cipher*& cvData, ZZX& poly, long& slotBits, long& xybatchBits, long& cnum, double& gamma, double& eta, double& etaprev, long& xyBits, long& wBits, long& pBits) {

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

 	double* coeffs = scheme.aux.taylorCoeffsMap.at(SIGMOIDPRIMEGOOD3);

 	Cipher cip2 = scheme.square(cip);
 	scheme.modSwitchAndEqual(cip2, wBits);
 	ZZ tmpzz = EvaluatorUtils::evaluateVal(coeffs[1], wBits);
 	Cipher cipc0c1 = scheme.multByConst(cip, tmpzz);

 	tmpzz = EvaluatorUtils::evaluateVal(coeffs[0], wBits);
 	scheme.modSwitchAndEqual(cipc0c1, wBits);
 	scheme.addConstAndEqual(cipc0c1, tmpzz);

 	Cipher* cgrad = new Cipher[cnum];

	RR tmprr = to_RR(gamma) * (1. - eta);
	tmpzz = EvaluatorUtils::evaluateVal(tmprr, wBits);

	tmprr = tmprr * coeffs[3];
	ZZ tmpzz2 = EvaluatorUtils::evaluateVal(tmprr, wBits);

	NTL_EXEC_RANGE(cnum, first, last);
	for (long i = first; i < last; ++i) {
		cprod[i] = scheme.multByConst(cxyData[i], tmpzz2);
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

	NTL_EXEC_RANGE(cnum, first, last);
	for (long i = first; i < last; ++i) {
		for (long l = xybatchBits; l < slotBits; ++l) {
			Cipher rot = scheme.leftRotateByPo2(cgrad[i], l);
			scheme.addAndEqual(cgrad[i], rot);
		}
	}
	NTL_EXEC_RANGE_END;

	tmprr = to_RR(1. - eta);
	tmpzz = EvaluatorUtils::evaluateVal(tmprr, wBits);
	bitsDown = xyBits + pBits + wBits;
	if(etaprev < 0) {
		tmprr = to_RR(eta) / (1 - etaprev);
		tmpzz2 = EvaluatorUtils::evaluateVal(tmprr, wBits);

		NTL_EXEC_RANGE(cnum, first, last);
		for (long i = first; i < last; ++i) {
			scheme.multByConstAndEqual(cwData[i], tmpzz);
			scheme.modSwitchAndEqual(cwData[i], wBits);
			scheme.modEmbedAndEqual(cwData[i], bitsDown);
			cgrad[i] = scheme.sub(cwData[i], cgrad[i]);

			scheme.multByConstAndEqual(cvData[i], tmpzz2);
			scheme.modSwitchAndEqual(cvData[i], wBits);
			scheme.modEmbedAndEqual(cvData[i], bitsDown);

			cwData[i] = scheme.add(cvData[i], cgrad[i]);
			cvData[i] = cgrad[i];
		}
		NTL_EXEC_RANGE_END;
	} else {
		NTL_EXEC_RANGE(cnum, first, last);
		for (long i = first; i < last; ++i) {
			cvData[i] = scheme.multByConst(cwData[i], tmpzz);
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

void CipherGD::decwData(double*& wData, Cipher*& cwData, long& factorDim, long& xyBatch, long& cnum, long& wBits) {
	for (long i = 0; i < (cnum - 1); ++i) {
		CZZ* dcw = scheme.decrypt(secretKey, cwData[i]);
		for (long l = 0; l < xyBatch; ++l) {
			RR wi = to_RR(dcw[l].r);
			wi.e -= wBits;
			wData[xyBatch * i + l] = to_double(wi);
		}
		delete[] dcw;
	}
	CZZ* dcw = scheme.decrypt(secretKey, cwData[cnum-1]);
	long rest = factorDim - xyBatch * (cnum - 1);
	for (long l = 0; l < rest; ++l) {
		RR wi = to_RR(dcw[l].r);
		wi.e -= wBits;
		wData[xyBatch * (cnum - 1) + l] = to_double(wi);
	}
	delete[] dcw;
}
