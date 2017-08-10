#include "CipherGD.h"

#include <NTL/ZZ.h>
#include <NTL/BasicThreadPool.h>

#include <Message.h>
#include <Cipher.h>
#include <CZZ.h>
#include <EvaluatorUtils.h>
#include <NTL/RR.h>
#include <NTL/ZZ.h>
#include <Params.h>
#include <SchemeAux.h>
#include <cmath>
#include <map>

Cipher* CipherGD::encxyData(long**& xyData, long& slots, long& factorDim, long& learnDim, long& learnDimPo2, long& xyBatch, long& cnum, long& xyBits) {
	Cipher* cxyData = new Cipher[cnum];
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
	}
	NTL_EXEC_INDEX_END;

	return cxyData;
}

Cipher* CipherGD::encwData(double*& wData, long& slots, long& factorDim, long& learnDim, long& learnDimPo2, long& xyBatch, long& cnum, long& wBits) {
	Cipher* cwData = new Cipher[cnum];
	NTL_EXEC_RANGE(cnum, first, last);
	for (long i = first; i < last; ++i) {
		CZZ* pwData = new CZZ[slots];
		if(i != (cnum - 1)) {
			for (long l = 0; l < xyBatch; ++l) {
				ZZ tmp = EvaluatorUtils::evaluateVal(wData[xyBatch * i + l], wBits);
				for (long j = 0; j < learnDim; ++j) {
					pwData[xyBatch * j + l] = CZZ(tmp);
				}
			}
		} else {
			long rest = factorDim - xyBatch * i;
			for (long l = 0; l < rest; ++l) {
				ZZ tmp = EvaluatorUtils::evaluateVal(wData[xyBatch * i + l], wBits);
				for (long j = 0; j < learnDim; ++j) {
					pwData[xyBatch * j + l] = CZZ(tmp);
				}
			}
		}
		cwData[i] = scheme.encrypt(pwData, slots);
	}
	NTL_EXEC_RANGE_END;
	return cwData;
}

void CipherGD::encStepNLGD5(Cipher*& cxyData, Cipher*& cwData, Cipher*& cvData, ZZX& poly, long& slots, long& learnDim, long learnDimPo2, long& xybatchBits, long& xyBatch, long& cnum, double& gamma, double& eta, double& etaprev, long& xyBits, long& wBits, long& pBits) {

	Cipher* cprod = new Cipher[cnum];
	long bitsDown = cxyData[0].cbits - cwData[0].cbits;

	NTL_EXEC_RANGE(cnum, first, last);
	for (long i = first; i < last; ++i) {
		cprod[i] = scheme.modEmbed(cxyData[i], bitsDown);
		scheme.multAndEqual(cprod[i], cwData[i]);
		for (long l = 0; l < xybatchBits; ++l) {
			Cipher rot = scheme.leftRotateByPo2(cprod[i], l);
			scheme.addAndEqual(cprod[i], rot);
		}
	}
	NTL_EXEC_RANGE_END;
	Cipher cip = algo.sum(cprod, cnum);
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
		scheme.modEmbed(cprod[i], bitsDown + pBits  + wBits);
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
	if(etaprev < 0) {
		RR cnst2 = to_RR(eta) / (1 - etaprev);
		ZZ pcnst2 = EvaluatorUtils::evaluateVal(cnst2, wBits);

		NTL_EXEC_RANGE(cnum, first, last);
		for (long i = first; i < last; ++i) {
			scheme.multByConstAndEqual(cwData[i], pcnst);
			scheme.modSwitchAndEqual(cwData[i], wBits);
			scheme.modEmbedAndEqual(cwData[i], xyBits + pBits + 2 * wBits);
			cgrad[i] = scheme.sub(cwData[i], cgrad[i]);

			scheme.multByConstAndEqual(cvData[i], pcnst2);
			scheme.modSwitchAndEqual(cvData[i], wBits);
			scheme.modEmbedAndEqual(cvData[i], xyBits + pBits + 2 * wBits);

			cwData[i] = scheme.add(cvData[i], cx	grad[i]);
			cvData[i] = cgrad[i];
		}
		NTL_EXEC_RANGE_END;
	} else {
		NTL_EXEC_RANGE(cnum, first, last);
		for (long i = first; i < last; ++i) {
			cvData[i] = scheme.multByConst(cwData[i], pcnst);
			scheme.modSwitchAndEqual(cvData[i], wBits);
			scheme.modEmbedAndEqual(cvData[i], xyBits + pBits + 2 * wBits);
			scheme.subAndEqual(cvData[i], cgrad[i]);
			scheme.modEmbedAndEqual(cwData[i], xyBits + pBits + 3 * wBits);
			scheme.subAndEqual(cwData[i], cgrad[i]);
		}
		NTL_EXEC_RANGE_END;
	}
}

void CipherGD::encStepNLGD3(Cipher*& cxyData, Cipher*& cwData, Cipher*& cvData, ZZX& poly, long& slots, long& learnDim, long learnDimPo2, long& xybatchBits, long& xyBatch, long& cnum, double& gamma, double& eta, double& etaprev, long& xyBits, long& wBits, long& pBits) {

}

double* CipherGD::decwData(SecKey& secretKey, Cipher*& cw, long& factorDim, long& xyBatch, long& cnum, long& wBits) {
	double* w = new double[factorDim];
	for (long i = 0; i < (cnum - 1); ++i) {
		CZZ* dcw = scheme.decrypt(secretKey, cw[i]);
		for (long l = 0; l < xyBatch; ++l) {
			RR wi = to_RR(dcw[l].r);
			wi.e -= wBits;
			w[xyBatch * i + l] = to_double(wi);
		}
	}
	CZZ* dcw = scheme.decrypt(secretKey, cw[cnum-1]);
	long rest = factorDim - xyBatch * (cnum - 1);
	for (long l = 0; l < rest; ++l) {
		RR wi = to_RR(dcw[l].r);
		wi.e -= wBits;
		w[xyBatch * (cnum - 1) + l] = to_double(wi);
	}
	return w;
}
