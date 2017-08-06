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

Cipher* CipherGD::encxyDataWB(long**& xyData, long& slots, long& factorDim, long& learnDim, long& wBatch) {
	Cipher* cxyData = new Cipher[factorDim];
	NTL_EXEC_RANGE(factorDim, first, last);
	for (long i = first; i < last; ++i) {
		CZZ* pxyData = new CZZ[slots];
		for (long j = 0; j < learnDim; ++j) {
			if(xyData[j][i] == -1) {
				for (long l = 0; l < wBatch; ++l) {
					pxyData[wBatch * j + l] = CZZ(-scheme.params.p);
				}
			} else if(xyData[j][i] == 1) {
				for (long l = 0; l < wBatch; ++l) {
					pxyData[wBatch * j + l] = CZZ(scheme.params.p);
				}
			}
		}
		cxyData[i] = scheme.encrypt(pxyData, slots);
	}
	NTL_EXEC_INDEX_END;
	return cxyData;
}

Cipher* CipherGD::encwDataWB(double**& wData, long& slots, long& factorDim, long& learnDim, long& wBatch) {
	Cipher* cwData = new Cipher[factorDim];
	NTL_EXEC_RANGE(factorDim, first, last);
	for (long i = first; i < last; ++i) {
		CZZ* pwData = new CZZ[slots];
		for (long l = 0; l < wBatch; ++l) {
			CZZ tmp = EvaluatorUtils::evaluateVal(wData[l][i], 0.0, scheme.params.logp);
			for (long j = 0; j < learnDim; ++j) {
				pwData[wBatch * j + l] = tmp;
			}
		}
		cwData[i] = scheme.encrypt(pwData, slots);
	}
	NTL_EXEC_RANGE_END;
	return cwData;
}

Cipher* CipherGD::encxyDataXYB(long**& xyData, long& slots, long& factorDim, long& learnDim, long& learnDimPo2, long& xyBatch, long& cnum) {
	Cipher* cxyData = new Cipher[cnum];
	NTL_EXEC_RANGE(cnum, first, last);
	for (long i = first; i < last; ++i) {
		CZZ* pxyData = new CZZ[slots];
		if(i != (cnum - 1)) {
			for (long l = 0; l < xyBatch; ++l) {
				for (long j = 0; j < learnDim; ++j) {
					if(xyData[j][xyBatch * i + l] == -1) {
						pxyData[xyBatch * j + l] = CZZ(-scheme.params.p);
					} else if(xyData[j][xyBatch * i + l] == 1) {
						pxyData[xyBatch * j + l] = CZZ(scheme.params.p);
					}
				}
			}
		} else {
			long rest = factorDim - xyBatch * i;
			for (long l = 0; l < rest; ++l) {
				for (long j = 0; j < learnDim; ++j) {
					if(xyData[j][xyBatch * i + l] == -1) {
						pxyData[xyBatch * j + l] = CZZ(-scheme.params.p);
					} else if(xyData[j][xyBatch * i + l] == 1) {
						pxyData[xyBatch * j + l] = CZZ(scheme.params.p);
					}
				}
			}
		}
		cxyData[i] = scheme.encrypt(pxyData, slots);
	}
	NTL_EXEC_INDEX_END;

	return cxyData;
}

Cipher* CipherGD::encwDataXYB(double*& wData, long& slots, long& factorDim, long& learnDim, long& learnDimPo2, long& xyBatch, long& cnum) {
	Cipher* cwData = new Cipher[cnum];
	NTL_EXEC_RANGE(cnum, first, last);
	for (long i = first; i < last; ++i) {
		CZZ* pwData = new CZZ[slots];
		if(i != (cnum - 1)) {
			for (long l = 0; l < xyBatch; ++l) {
				CZZ tmp = EvaluatorUtils::evaluateVal(wData[xyBatch * i + l], 0.0, scheme.params.logp);
				for (long j = 0; j < learnDim; ++j) {
					pwData[xyBatch * j + l] = tmp;
				}
			}
		} else {
			long rest = factorDim - xyBatch * i;
			for (long l = 0; l < rest; ++l) {
				CZZ tmp = EvaluatorUtils::evaluateVal(wData[xyBatch * i + l], 0.0, scheme.params.logp);
				for (long j = 0; j < learnDim; ++j) {
					pwData[xyBatch * j + l] = tmp;
				}
			}
		}
		cwData[i] = scheme.encrypt(pwData, slots);
	}
	NTL_EXEC_RANGE_END;
	return cwData;
}

void CipherGD::encStepLGDWB(Cipher*& cxyData, Cipher*& cwData, long& slots, long& factorDim, long& learnDim, long& wBatch, double& lambda, double& gamma) {

	Cipher* cprod = new Cipher[factorDim];

	NTL_EXEC_RANGE(factorDim, first, last);
	for (long i = first; i < last; ++i) {
		cprod[i] = scheme.modEmbed(cxyData[i], cwData[i].level);
		scheme.multAndEqual(cprod[i], cwData[i]);
	}
	NTL_EXEC_RANGE_END;

	Cipher cip = algo.sum(cprod, factorDim);

	scheme.modSwitchOneAndEqual(cip); // cip (-1)

	RR cnst = 1.0 - to_RR(gamma) * lambda; // 1 - gamma * lambda
	ZZ pcnst = pmult(cnst); // p * (1 - gamma * lambda)

	for (long i = 0; i < factorDim; ++i) {
		scheme.multByConstAndEqual(cwData[i], pcnst); // (1 - gamma * lambda) * w
		scheme.modSwitchOneAndEqual(cwData[i]); // (1 - gamma * lambda) * w  (-1)
	}

	double* coeffs = scheme.aux.taylorCoeffsMap.at(SIGMOIDPRIMEGOOD7);
	Cipher* cpows = algo.powerOf2Extended(cip, 2); // ip (-1), ip^2 (-2), ip^4 (-3)

	cnst = to_RR(gamma) * coeffs[0]; // gamma * coeff_0
	pcnst = pmult(cnst); // p * gamma * coeff_0

	Cipher* cgrad = new Cipher[factorDim];
	NTL_EXEC_RANGE(factorDim, first, last);
	for (long i = first; i < last; ++i) {
		cgrad[i] = scheme.multByConst(cxyData[i], pcnst); // p * p * gamma * (z * coeff_0)
		scheme.modSwitchOneAndEqual(cgrad[i]); // p * gamma * (z * coeff_0) (-1)
	}
	NTL_EXEC_RANGE_END;

	for (long t = 1; t < 8; t=t+2) {
		cnst = to_RR(gamma) * coeffs[t]; // gamma * coeff_0
		pcnst = pmult(cnst); // p * gamma * coeff_0

		NTL_EXEC_RANGE(factorDim, first, last);
		for (long i = first; i < last; ++i) {
			Cipher cgradit = scheme.multByConst(cxyData[i], pcnst); // p * p * gamma * (z * coeff_t)
			scheme.modSwitchOneAndEqual(cgradit); // p * gamma * (z * coeff_t) (-1)
			if(bit(t, 0)) {
				scheme.modEmbedAndEqual(cgradit, cpows[0].level);
				scheme.multModSwitchOneAndEqual(cgradit, cpows[0]); // p * gamma * (z * coeff_t * ip^t)
			}
			if(bit(t, 1)) {
				scheme.modEmbedAndEqual(cgradit, cpows[1].level);
				scheme.multModSwitchOneAndEqual(cgradit, cpows[1]); // p * gamma * (z * coeff_t * ip^t)
			}
			if(bit(t, 2)) {
				scheme.modEmbedAndEqual(cgradit, cpows[2].level);
				scheme.multModSwitchOneAndEqual(cgradit, cpows[2]); // p * gamma * (z * coeff_t * ip^t)
			}
			scheme.modEmbedAndEqual(cgrad[i], cgradit.level);
			scheme.addAndEqual(cgrad[i], cgradit); // p * gamma * (z * sigmoid(ip)) (-4)
		}
		NTL_EXEC_RANGE_END;
	}

	long logslots = log2(slots);
	long logwnum = log2(wBatch);

	NTL_EXEC_RANGE(factorDim, first, last);
	for (long i = first; i < last; ++i) {
		for (long l = logwnum; l < logslots; ++l) {
			Cipher rot = scheme.leftRotateByPo2(cgrad[i], l);
			scheme.addAndEqual(cgrad[i], rot); // p * gamma * sum(z * sigmoid(ip)) (-4)
		}
		scheme.modEmbedAndEqual(cwData[i], cgrad[i].level);
		scheme.subAndEqual(cwData[i], cgrad[i]); // w - gamma * grad(w) - gamma * lambda * w (-4)
	}
	NTL_EXEC_RANGE_END;

}

void CipherGD::encStepNLGD7WB(Cipher*& cxyData, Cipher*& cwData, Cipher*& cvData, long& slots, long& factorDim, long& learnDim, long& wBatch, double& gamma, double& eta) {

	Cipher* cprod = new Cipher[factorDim];

	NTL_EXEC_RANGE(factorDim, first, last);
	for (long i = first; i < last; ++i) {
		cprod[i] = scheme.modEmbed(cxyData[i], cwData[i].level);
		scheme.multAndEqual(cprod[i], cwData[i]);
	}
	NTL_EXEC_RANGE_END;

	Cipher cip = algo.sum(cprod, factorDim);

	scheme.modSwitchOneAndEqual(cip); // cip (-1)

	double* coeffs = scheme.aux.taylorCoeffsMap.at(SIGMOIDPRIMEGOOD7);
	Cipher* cpows = algo.powerOf2Extended(cip, 2); // ip (-1), ip^2 (-2), ip^4 (-3)

	RR cnst = to_RR(gamma) * coeffs[0]; // gamma * coeff_0
	ZZ pcnst = pmult(cnst); // p * gamma * coeff_0

	Cipher* cgrad = new Cipher[factorDim];
	NTL_EXEC_RANGE(factorDim, first, last);
	for (long i = first; i < last; ++i) {
		cgrad[i] = scheme.multByConst(cxyData[i], pcnst); // p * p * gamma * (xy * coeff_0)
		scheme.modSwitchOneAndEqual(cgrad[i]); // p * gamma * (xy * coeff_0) (-1)
	}
	NTL_EXEC_RANGE_END;

	for (long t = 1; t < 8; t=t+2) {
		cnst = to_RR(gamma) * coeffs[t]; // gamma * coeff_t
		pcnst = pmult(cnst); // p * gamma * coeff_t

		NTL_EXEC_RANGE(factorDim, first, last);
		for (long i = first; i < last; ++i) {
			Cipher cgradit = scheme.multByConst(cxyData[i], pcnst); // p * p * gamma * (xy * coeff_t)
			scheme.modSwitchOneAndEqual(cgradit); // p * gamma * (xy * coeff_t) (-1)
			if(bit(t, 0)) {
				scheme.modEmbedAndEqual(cgradit, cpows[0].level);
				scheme.multModSwitchOneAndEqual(cgradit, cpows[0]); // p * gamma * (xy * coeff_t * ip^t)
			}
			if(bit(t, 1)) {
				scheme.modEmbedAndEqual(cgradit, cpows[1].level);
				scheme.multModSwitchOneAndEqual(cgradit, cpows[1]); // p * gamma * (xy * coeff_t * ip^t)
			}
			if(bit(t, 2)) {
				scheme.modEmbedAndEqual(cgradit, cpows[2].level);
				scheme.multModSwitchOneAndEqual(cgradit, cpows[2]); // p * gamma * (xy * coeff_t * ip^t)
			}
			scheme.modEmbedAndEqual(cgrad[i], cgradit.level);
			scheme.addAndEqual(cgrad[i], cgradit); // p * gamma * (xy * sigmoid(ip)) (-4)
		}
		NTL_EXEC_RANGE_END;
	}

	long logslots = log2(slots);
	long logwnum = log2(wBatch);

	NTL_EXEC_RANGE(factorDim, first, last);
	for (long i = first; i < last; ++i) {
		for (long l = logwnum; l < logslots; ++l) {
			Cipher rot = scheme.leftRotateByPo2(cgrad[i], l);
			scheme.addAndEqual(cgrad[i], rot); // p * gamma * sum(xy * sigmoid(ip)) (-4)
		}
	}
	NTL_EXEC_RANGE_END;

	cnst = to_RR(eta);
	pcnst = pmult(cnst);

	NTL_EXEC_RANGE(factorDim, first, last);
	for (long i = first; i < last; ++i) {
		scheme.modEmbedAndEqual(cwData[i], cgrad[i].level);
		cgrad[i] = scheme.sub(cwData[i], cgrad[i]); // tmp = w - gamma * grad (-4)
		scheme.modEmbedAndEqual(cvData[i], cgrad[i].level);
		cwData[i] = scheme.sub(cvData[i], cgrad[i]); // w = v - tmp (-4)
		scheme.multByConstAndEqual(cwData[i], pcnst); // w = p * eta * (v - tmp) (-4)
		scheme.modSwitchOneAndEqual(cwData[i]); // w = eta * (v - tmp) (-5)
		cvData[i] = scheme.modEmbed(cgrad[i], cwData[i].level); // v = tmp (-5)
		scheme.addAndEqual(cwData[i], cvData[i]); //  w = tmp + eta * (v - tmp) (-5)
	}
	NTL_EXEC_RANGE_END;
}

void CipherGD::encStepNLGD7XYB6(Cipher*& cxyData, Cipher*& cwData, Cipher*& cvData, ZZX& poly, long& slots, long& learnDim, long learnDimPo2, long& xybatchBits, long& xyBatch, long& cnum, double& gamma, double& eta) {

	Cipher* cprod = new Cipher[cnum];

	NTL_EXEC_RANGE(cnum, first, last);
	for (long i = first; i < last; ++i) {
		cprod[i] = scheme.modEmbed(cxyData[i], cwData[i].level);
		scheme.multAndEqual(cprod[i], cwData[i]);
		for (long l = 0; l < xybatchBits; ++l) {
			Cipher rot = scheme.leftRotateByPo2(cprod[i], l);
			scheme.addAndEqual(cprod[i], rot);
		}
	}
	NTL_EXEC_RANGE_END;
	Cipher cip = algo.sum(cprod, cnum);
	scheme.modSwitchOneAndEqual(cip);

	scheme.multByPolyAndEqual(cip, poly);
	for (long l = 0; l < xybatchBits; ++l) {
		Cipher rot = scheme.rightRotateByPo2(cip, l);
		scheme.addAndEqual(cip, rot);
	}
	scheme.modSwitchOneAndEqual(cip);

 	double* coeffs = scheme.aux.taylorCoeffsMap.at(SIGMOIDPRIMEGOOD7);
	Cipher* cpows = algo.powerOf2Extended(cip, 2); // ip (-1), ip^2 (-2), ip^4 (-3)

	RR cnst = to_RR(gamma) * coeffs[0]; // gamma * coeff_0
	ZZ pcnst = pmult(cnst); // p * gamma * coeff_0

	Cipher* cgrad = new Cipher[cnum];
	NTL_EXEC_RANGE(cnum, first, last);
	for (long i = first; i < last; ++i) {
		cgrad[i] = scheme.multByConst(cxyData[i], pcnst); // p * p * gamma * (xy * coeff_0)
		scheme.modSwitchOneAndEqual(cgrad[i]); // p * gamma * (xy * coeff_0) (-1)
	}
	NTL_EXEC_RANGE_END;

	for (long t = 1; t < 8; t=t+2) {

		cnst = to_RR(gamma) * coeffs[t]; // gamma * coeff_t
		pcnst = pmult(cnst); // p * gamma * coeff_t

		NTL_EXEC_RANGE(cnum, first, last);
		for (long i = first; i < last; ++i) {
			Cipher cgradit = scheme.multByConst(cxyData[i], pcnst); // p * p * gamma * (xy * coeff_t)
			scheme.modSwitchOneAndEqual(cgradit); // p * gamma * (xy * coeff_t) (-1)
			if(bit(t, 0)) {
				scheme.modEmbedAndEqual(cgradit, cpows[0].level);
				scheme.multModSwitchOneAndEqual(cgradit, cpows[0]); // p * gamma * (xy * coeff_t * ip^t)
			}
			if(bit(t, 1)) {
				scheme.modEmbedAndEqual(cgradit, cpows[1].level);
				scheme.multModSwitchOneAndEqual(cgradit, cpows[1]); // p * gamma * (xy * coeff_t * ip^t)
			}
			if(bit(t, 2)) {
				scheme.modEmbedAndEqual(cgradit, cpows[2].level);
				scheme.multModSwitchOneAndEqual(cgradit, cpows[2]); // p * gamma * (xy * coeff_t * ip^t)
			}
			scheme.modEmbedAndEqual(cgrad[i], cgradit.level);
			scheme.addAndEqual(cgrad[i], cgradit); // p * gamma * (xy * sigmoid(ip)) (-4)
		}
		NTL_EXEC_RANGE_END;
	}

	long logslots = log2(slots);

	NTL_EXEC_RANGE(cnum, first, last);
	for (long i = first; i < last; ++i) {
		for (long l = xybatchBits; l < logslots; ++l) {
			Cipher rot = scheme.leftRotateByPo2(cgrad[i], l);
			scheme.addAndEqual(cgrad[i], rot); // p * gamma * sum(xy * sigmoid(ip)) (-4)
		}
	}
	NTL_EXEC_RANGE_END;

	cnst = to_RR(eta);
	pcnst = pmult(cnst);

	NTL_EXEC_RANGE(cnum, first, last);
	for (long i = first; i < last; ++i) {
		scheme.modEmbedAndEqual(cwData[i], cgrad[i].level);
		cgrad[i] = scheme.sub(cwData[i], cgrad[i]); // tmp = w - gamma * grad (-4)
		scheme.modEmbedAndEqual(cvData[i], cgrad[i].level);
		cwData[i] = scheme.sub(cvData[i], cgrad[i]); // w = v - tmp (-4)
		scheme.multByConstAndEqual(cwData[i], pcnst); // w = p * eta * (v - tmp) (-4)
		scheme.modSwitchOneAndEqual(cwData[i]); // w = eta * (v - tmp) (-5)
		cvData[i] = scheme.modEmbed(cgrad[i], cwData[i].level); // v = tmp (-5)
		scheme.addAndEqual(cwData[i], cvData[i]); //  w = tmp + eta * (v - tmp) (-5)
	}
	NTL_EXEC_RANGE_END;
}

void CipherGD::encStepNLGD7XYB5(Cipher*& cxyData, Cipher*& cwData, Cipher*& cvData, ZZX& poly, long& slots, long& learnDim, long learnDimPo2, long& xybatchBits, long& xyBatch, long& cnum, double& gamma, double& eta, double& etaprev) {

	Cipher* cprod = new Cipher[cnum];

	NTL_EXEC_RANGE(cnum, first, last);
	for (long i = first; i < last; ++i) {
		cprod[i] = scheme.modEmbed(cxyData[i], cwData[i].level);
		scheme.multAndEqual(cprod[i], cwData[i]);
		for (long l = 0; l < xybatchBits; ++l) {
			Cipher rot = scheme.leftRotateByPo2(cprod[i], l);
			scheme.addAndEqual(cprod[i], rot);
		}
	}
	NTL_EXEC_RANGE_END;
	Cipher cip = algo.sum(cprod, cnum);
	scheme.modSwitchOneAndEqual(cip);

	scheme.multByPolyAndEqual(cip, poly);
	for (long l = 0; l < xybatchBits; ++l) {
		Cipher rot = scheme.rightRotateByPo2(cip, l);
		scheme.addAndEqual(cip, rot);
	}
	scheme.modSwitchOneAndEqual(cip);

 	double* coeffs = scheme.aux.taylorCoeffsMap.at(SIGMOIDPRIMEGOOD7);
	Cipher* cpows = algo.powerOf2Extended(cip, 2); // ip (-1), ip^2 (-2), ip^4 (-3)

	RR cnst = to_RR(gamma) * coeffs[0] * (1. - eta); // gamma * coeff_0
	ZZ pcnst = pmult(cnst); // p * gamma * coeff_0

	Cipher* cgrad = new Cipher[cnum];
	NTL_EXEC_RANGE(cnum, first, last);
	for (long i = first; i < last; ++i) {
		cgrad[i] = scheme.multByConst(cxyData[i], pcnst); // p * p * gamma * (xy * coeff_0)
		scheme.modSwitchOneAndEqual(cgrad[i]); // p * gamma * (xy * coeff_0) (-1)
	}
	NTL_EXEC_RANGE_END;

	for (long t = 1; t < 8; t=t+2) {

		cnst = to_RR(gamma) * coeffs[t] * (1. - eta); // gamma * coeff_t
		pcnst = pmult(cnst); // p * gamma * coeff_t

		NTL_EXEC_RANGE(cnum, first, last);
		for (long i = first; i < last; ++i) {
			Cipher cgradit = scheme.multByConst(cxyData[i], pcnst); // p * p * gamma * (xy * coeff_t)
			scheme.modSwitchOneAndEqual(cgradit); // p * gamma * (xy * coeff_t) (-1)
			if(bit(t, 0)) {
				scheme.modEmbedAndEqual(cgradit, cpows[0].level);
				scheme.multModSwitchOneAndEqual(cgradit, cpows[0]); // p * gamma * (xy * coeff_t * ip^t)
			}
			if(bit(t, 1)) {
				scheme.modEmbedAndEqual(cgradit, cpows[1].level);
				scheme.multModSwitchOneAndEqual(cgradit, cpows[1]); // p * gamma * (xy * coeff_t * ip^t)
			}
			if(bit(t, 2)) {
				scheme.modEmbedAndEqual(cgradit, cpows[2].level);
				scheme.multModSwitchOneAndEqual(cgradit, cpows[2]); // p * gamma * (xy * coeff_t * ip^t)
			}
			scheme.modEmbedAndEqual(cgrad[i], cgradit.level);
			scheme.addAndEqual(cgrad[i], cgradit); // p * gamma * (xy * sigmoid(ip)) (-4)
		}
		NTL_EXEC_RANGE_END;
	}

	long logslots = log2(slots);

	NTL_EXEC_RANGE(cnum, first, last);
	for (long i = first; i < last; ++i) {
		for (long l = xybatchBits; l < logslots; ++l) {
			Cipher rot = scheme.leftRotateByPo2(cgrad[i], l);
			scheme.addAndEqual(cgrad[i], rot); // p * gamma * sum(xy * sigmoid(ip)) (-4)
		}
	}
	NTL_EXEC_RANGE_END;
	cnst = to_RR(1. - eta);
	pcnst = pmult(cnst);
	if(etaprev < 0) {
		RR cnst2 = to_RR(eta) / (1 - etaprev);
		ZZ pcnst2 = pmult(cnst2);

		NTL_EXEC_RANGE(cnum, first, last);
		for (long i = first; i < last; ++i) {
			scheme.multByConstAndEqual(cwData[i], pcnst); // w = (1 - eta) * w
			scheme.modSwitchOneAndEqual(cwData[i]);
			scheme.modEmbedAndEqual(cwData[i], cgrad[i].level);
			cgrad[i] = scheme.sub(cwData[i], cgrad[i]); // tmp = (1 - eta) * w - (1 - eta) * gamma * grad

			scheme.multByConstAndEqual(cvData[i], pcnst2);
			scheme.modSwitchOneAndEqual(cvData[i]); // v = eta / (1 + etaprev) * v
			scheme.modEmbedAndEqual(cvData[i], cgrad[i].level);

			cwData[i] = scheme.add(cvData[i], cgrad[i]); // w = (1 - eta) * w - (1 - eta) * gamma * grad + eta / (1 - etaprev) * v
			cvData[i] = cgrad[i]; // v = (1 - eta) * w - (1 - eta) * gamma * grad
		}
		NTL_EXEC_RANGE_END;
	} else {
		NTL_EXEC_RANGE(cnum, first, last);
		for (long i = first; i < last; ++i) {
			cvData[i] = scheme.multByConst(cwData[i], pcnst); // tmp2 = (1 - eta) * w
			scheme.modSwitchOneAndEqual(cvData[i]);
			scheme.modEmbedAndEqual(cvData[i], cgrad[i].level);
			scheme.subAndEqual(cvData[i], cgrad[i]); // v = (1 - eta) * w - (1 - eta) * gamma * grad
			scheme.modEmbedAndEqual(cwData[i], cgrad[i].level);
			scheme.subAndEqual(cwData[i], cgrad[i]);
		}
		NTL_EXEC_RANGE_END;
	}
}

void CipherGD::encStepNLGD3XYB4(Cipher*& cxyData, Cipher*& cwData, Cipher*& cvData, ZZX& poly, long& slots, long& learnDim, long learnDimPo2, long& xybatchBits, long& xyBatch, long& cnum, double& gamma, double& eta, double& etaprev) {

	Cipher* cprod = new Cipher[cnum];

	NTL_EXEC_RANGE(cnum, first, last);
	for (long i = first; i < last; ++i) {
		cprod[i] = scheme.modEmbed(cxyData[i], cwData[i].level);
		scheme.multAndEqual(cprod[i], cwData[i]);
		for (long l = 0; l < xybatchBits; ++l) {
			Cipher rot = scheme.leftRotateByPo2(cprod[i], l);
			scheme.addAndEqual(cprod[i], rot);
		}
	}
	NTL_EXEC_RANGE_END;
	Cipher cip = algo.sum(cprod, cnum);
	scheme.modSwitchOneAndEqual(cip);

	scheme.multByPolyAndEqual(cip, poly);
	for (long l = 0; l < xybatchBits; ++l) {
		Cipher rot = scheme.rightRotateByPo2(cip, l);
		scheme.addAndEqual(cip, rot);
	}
	scheme.modSwitchOneAndEqual(cip);

 	double* coeffs = scheme.aux.taylorCoeffsMap.at(SIGMOIDPRIMEGOOD3);
	Cipher* cpows = algo.powerOf2Extended(cip, 1); // ip (-1), ip^2 (-2)

	RR cnst = to_RR(gamma) * coeffs[0] * (1. - eta); // gamma * coeff_0
	ZZ pcnst = pmult(cnst); // p * gamma * coeff_0

	Cipher* cgrad = new Cipher[cnum];
	NTL_EXEC_RANGE(cnum, first, last);
	for (long i = first; i < last; ++i) {
		cgrad[i] = scheme.multByConst(cxyData[i], pcnst); // p * p * gamma * (xy * coeff_0)
		scheme.modSwitchOneAndEqual(cgrad[i]); // p * gamma * (xy * coeff_0) (-1)
	}
	NTL_EXEC_RANGE_END;

	for (long t = 1; t < 4; t=t+2) {

		cnst = to_RR(gamma) * coeffs[t] * (1. - eta); // gamma * coeff_t
		pcnst = pmult(cnst); // p * gamma * coeff_t

		NTL_EXEC_RANGE(cnum, first, last);
		for (long i = first; i < last; ++i) {
			Cipher cgradit = scheme.multByConst(cxyData[i], pcnst); // p * p * gamma * (xy * coeff_t)
			scheme.modSwitchOneAndEqual(cgradit); // p * gamma * (xy * coeff_t) (-1)
			if(bit(t, 0)) {
				scheme.modEmbedAndEqual(cgradit, cpows[0].level);
				scheme.multModSwitchOneAndEqual(cgradit, cpows[0]); // p * gamma * (xy * coeff_t * ip^t)
			}
			if(bit(t, 1)) {
				scheme.modEmbedAndEqual(cgradit, cpows[1].level);
				scheme.multModSwitchOneAndEqual(cgradit, cpows[1]); // p * gamma * (xy * coeff_t * ip^t)
			}
			scheme.modEmbedAndEqual(cgrad[i], cgradit.level);
			scheme.addAndEqual(cgrad[i], cgradit); // p * gamma * (xy * sigmoid(ip)) (-3)
		}
		NTL_EXEC_RANGE_END;
	}

	long logslots = log2(slots);

	NTL_EXEC_RANGE(cnum, first, last);
	for (long i = first; i < last; ++i) {
		for (long l = xybatchBits; l < logslots; ++l) {
			Cipher rot = scheme.leftRotateByPo2(cgrad[i], l);
			scheme.addAndEqual(cgrad[i], rot); // p * gamma * sum(xy * sigmoid(ip)) (-4)
		}
	}
	NTL_EXEC_RANGE_END;
	cnst = to_RR(1. - eta);
	pcnst = pmult(cnst);
	if(etaprev < 0) {
		RR cnst2 = to_RR(eta) / (1 - etaprev);
		ZZ pcnst2 = pmult(cnst2);

		NTL_EXEC_RANGE(cnum, first, last);
		for (long i = first; i < last; ++i) {
			scheme.multByConstAndEqual(cwData[i], pcnst); // w = (1 - eta) * w
			scheme.modSwitchOneAndEqual(cwData[i]);
			scheme.modEmbedAndEqual(cwData[i], cgrad[i].level);
			cgrad[i] = scheme.sub(cwData[i], cgrad[i]); // tmp = (1 - eta) * w - (1 - eta) * gamma * grad

			scheme.multByConstAndEqual(cvData[i], pcnst2);
			scheme.modSwitchOneAndEqual(cvData[i]); // v = eta / (1 + etaprev) * v
			scheme.modEmbedAndEqual(cvData[i], cgrad[i].level);

			cwData[i] = scheme.add(cvData[i], cgrad[i]); // w = (1 - eta) * w - (1 - eta) * gamma * grad + eta / (1 - etaprev) * v
			cvData[i] = cgrad[i]; // v = (1 - eta) * w - (1 - eta) * gamma * grad
		}
		NTL_EXEC_RANGE_END;
	} else {
		NTL_EXEC_RANGE(cnum, first, last);
		for (long i = first; i < last; ++i) {
			cvData[i] = scheme.multByConst(cwData[i], pcnst); // tmp2 = (1 - eta) * w
			scheme.modSwitchOneAndEqual(cvData[i]);
			scheme.modEmbedAndEqual(cvData[i], cgrad[i].level);
			scheme.subAndEqual(cvData[i], cgrad[i]); // v = (1 - eta) * w - (1 - eta) * gamma * grad
			scheme.modEmbedAndEqual(cwData[i], cgrad[i].level);
			scheme.subAndEqual(cwData[i], cgrad[i]);
		}
		NTL_EXEC_RANGE_END;
	}
}

void CipherGD::encStepNLGD3XYB4(Cipher*& cxyData, Cipher*& cwData, Cipher*& cvData, ZZX& poly, long& slots, long& learnDim, long learnDimPo2, long& xybatchBits, long& xyBatch, long& cnum, double& gamma, double& eta, double& etaprev) {

	Cipher* cprod = new Cipher[cnum];

	NTL_EXEC_RANGE(cnum, first, last);
	for (long i = first; i < last; ++i) {
		cprod[i] = scheme.modEmbed(cxyData[i], cwData[i].level);
		scheme.multAndEqual(cprod[i], cwData[i]);
		for (long l = 0; l < xybatchBits; ++l) {
			Cipher rot = scheme.leftRotateByPo2(cprod[i], l);
			scheme.addAndEqual(cprod[i], rot);
		}
	}
	NTL_EXEC_RANGE_END;
	Cipher cip = algo.sum(cprod, cnum);
	scheme.modSwitchOneAndEqual(cip);

	scheme.multByPolyAndEqual(cip, poly);
	for (long l = 0; l < xybatchBits; ++l) {
		Cipher rot = scheme.rightRotateByPo2(cip, l);
		scheme.addAndEqual(cip, rot);
	}
	scheme.modSwitchOneAndEqual(cip);

	RR cnst = to_RR(gamma) * 0.5 * (1. - eta); // gamma * coeff_0
	ZZ pcnst = pmult(cnst); // p * gamma * coeff_0

	Cipher* cgrad = new Cipher[cnum];
	NTL_EXEC_RANGE(cnum, first, last);
	for (long i = first; i < last; ++i) {
		cgrad[i] = scheme.multByConst(cxyData[i], pcnst); // p * p * gamma * (xy * coeff_0)
		scheme.modSwitchOneAndEqual(cgrad[i]); // p * gamma * (xy * coeff_0) (-1)
	}
	NTL_EXEC_RANGE_END;

	cnst = to_RR(gamma) * 0.20 * (1. - eta); // gamma * coeff_t
	pcnst = pmult(cnst); // p * gamma * coeff_t

	NTL_EXEC_RANGE(cnum, first, last);
	for (long i = first; i < last; ++i) {
		Cipher cgradit = scheme.multByConst(cxyData[i], pcnst); // p * p * gamma * (xy * coeff_t)
		scheme.modSwitchOneAndEqual(cgradit); // p * gamma * (xy * coeff_t) (-1)
		scheme.modEmbedAndEqual(cgradit, cip.level);
		scheme.multModSwitchOneAndEqual(cgradit, cip); // p * gamma * (xy * coeff_t * ip^t)
		scheme.modEmbedAndEqual(cgrad[i], cgradit.level);
		scheme.addAndEqual(cgrad[i], cgradit); // p * gamma * (xy * sigmoid(ip)) (-3)
	}
	NTL_EXEC_RANGE_END;

	long logslots = log2(slots);

	NTL_EXEC_RANGE(cnum, first, last);
	for (long i = first; i < last; ++i) {
		for (long l = xybatchBits; l < logslots; ++l) {
			Cipher rot = scheme.leftRotateByPo2(cgrad[i], l);
			scheme.addAndEqual(cgrad[i], rot); // p * gamma * sum(xy * sigmoid(ip)) (-4)
		}
	}
	NTL_EXEC_RANGE_END;
	cnst = to_RR(1. - eta);
	pcnst = pmult(cnst);
	if(etaprev < 0) {
		RR cnst2 = to_RR(eta) / (1 - etaprev);
		ZZ pcnst2 = pmult(cnst2);

		NTL_EXEC_RANGE(cnum, first, last);
		for (long i = first; i < last; ++i) {
			scheme.multByConstAndEqual(cwData[i], pcnst); // w = (1 - eta) * w
			scheme.modSwitchOneAndEqual(cwData[i]);
			scheme.modEmbedAndEqual(cwData[i], cgrad[i].level);
			cgrad[i] = scheme.sub(cwData[i], cgrad[i]); // tmp = (1 - eta) * w - (1 - eta) * gamma * grad

			scheme.multByConstAndEqual(cvData[i], pcnst2);
			scheme.modSwitchOneAndEqual(cvData[i]); // v = eta / (1 + etaprev) * v
			scheme.modEmbedAndEqual(cvData[i], cgrad[i].level);

			cwData[i] = scheme.add(cvData[i], cgrad[i]); // w = (1 - eta) * w - (1 - eta) * gamma * grad + eta / (1 - etaprev) * v
			cvData[i] = cgrad[i]; // v = (1 - eta) * w - (1 - eta) * gamma * grad
		}
		NTL_EXEC_RANGE_END;
	} else {
		NTL_EXEC_RANGE(cnum, first, last);
		for (long i = first; i < last; ++i) {
			cvData[i] = scheme.multByConst(cwData[i], pcnst); // tmp2 = (1 - eta) * w
			scheme.modSwitchOneAndEqual(cvData[i]);
			scheme.modEmbedAndEqual(cvData[i], cgrad[i].level);
			scheme.subAndEqual(cvData[i], cgrad[i]); // v = (1 - eta) * w - (1 - eta) * gamma * grad
			scheme.modEmbedAndEqual(cwData[i], cgrad[i].level);
			scheme.subAndEqual(cwData[i], cgrad[i]);
		}
		NTL_EXEC_RANGE_END;
	}
}

void CipherGD::encwsumWB(Cipher*& cwData, long& factorDim, long& wBatch) {
	for (long i = 0; i < factorDim; ++i) {
		algo.partialSlotsSumAndEqual(cwData[i], wBatch);
	}
}

double* CipherGD::decWB(SecKey& secretKey, Cipher*& cw, long& factorDim) {
	double* w = new double[factorDim];
	for (long i = 0; i < factorDim; ++i) {
		CZZ* dcw = scheme.decrypt(secretKey, cw[i]);
		RR wi = to_RR(dcw[0].r);
		wi.e -= scheme.params.logp;
		w[i] = to_double(wi);
	}
	return w;
}

double* CipherGD::decXYB(SecKey& secretKey, Cipher*& cw, long& factorDim, long& xyBatch, long& cnum) {
	double* w = new double[factorDim];
	for (long i = 0; i < (cnum - 1); ++i) {
		CZZ* dcw = scheme.decrypt(secretKey, cw[i]);
		for (long l = 0; l < xyBatch; ++l) {
			RR wi = to_RR(dcw[l].r);
			wi.e -= scheme.params.logp;
			w[xyBatch * i + l] = to_double(wi);
		}
	}
	CZZ* dcw = scheme.decrypt(secretKey, cw[cnum-1]);
	long rest = factorDim - xyBatch * (cnum - 1);
	for (long l = 0; l < rest; ++l) {
		RR wi = to_RR(dcw[l].r);
		wi.e -= scheme.params.logp;
		w[xyBatch * (cnum - 1) + l] = to_double(wi);
	}
	return w;
}

ZZ CipherGD::pmult(RR val) {
	val.e += scheme.params.logp;
	return to_ZZ(val);
}

void CipherGD::debugcheck(string prefix, SecKey& secretKey, Cipher*& ciphers, long factorCheck, long slotCheck) {
	cout << prefix << " " << ciphers[0].level << endl;
	for (long i = 0; i < factorCheck; ++i) {
		CZZ* deci = scheme.decrypt(secretKey, ciphers[i]);
		for (int j = 0; j < slotCheck; ++j) {
			RR wir = to_RR(deci[j].r);
			wir.e -= scheme.params.logp;
			double wi = to_double(wir);
			cout << wi << ",";
		}
		cout << endl;
	}
	cout << endl;
}

void CipherGD::debugcheck(string prefix, SecKey& secretKey, Cipher& cipher, long slotCheck) {
	cout << prefix << " " << cipher.level << endl;
	CZZ* dec = scheme.decrypt(secretKey, cipher);
	for (long j = 0; j < slotCheck; ++j) {
		RR wr = to_RR(dec[j].r);
		wr.e -= scheme.params.logp;
		double w = to_double(wr);
		cout << w << ",";
	}
	cout << endl;
}
