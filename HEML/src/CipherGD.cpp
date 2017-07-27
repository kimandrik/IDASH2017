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

Cipher* CipherGD::encxyData(long**& xyData, long& slots, long& factorDim, long& learnDim, long& wBatch) {
	Cipher* cxyData = new Cipher[factorDim];
	NTL_EXEC_RANGE(factorDim, first, last);
	for (long i = first; i < last; ++i) {
		CZZ* pxyData = new CZZ[slots];
		for (long j = 0; j < learnDim; ++j) {
			if(xyData[i][j] == -1) {
				for (long l = 0; l < wBatch; ++l) {
					pxyData[wBatch * j + l] = CZZ(-scheme.params.p);
				}
			} else if(xyData[i][j] == 1) {
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

Cipher* CipherGD::encwData(double**& wData, long& slots, long& factorDim, long& learnDim, long& wBatch) {
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

void CipherGD::encStepQGD(Cipher*& cxyData, Cipher*& cwData, long& slots, long& factorDim, long& learnDim, long& wBatch, double& gamma, double& lambda) {

	long dimcheck = 5;
	long slotscheck = 10;

	debugcheck("c xyData: ", secretKey, cxyData, dimcheck, slotscheck);

	debugcheck("c wData: ", secretKey, cwData, dimcheck, slotscheck);

	Cipher* cprod = new Cipher[factorDim];
	NTL_EXEC_RANGE(factorDim, first, last);
	for (long i = first; i < last; ++i) {
		cout << i << endl;
		cprod[i] = scheme.modEmbed(cxyData[i], cwData[i].level);
		scheme.multAndEqual(cprod[i], cwData[i]);
	}
	NTL_EXEC_RANGE_END;

	Cipher cip = algo.sum(cprod, factorDim);

	scheme.modSwitchOneAndEqual(cip);

	debugcheck("c inner prod:", secretKey, cip, slotscheck);

	RR cnst = 1.0 - to_RR(gamma) * lambda; // 1 - gamma * lambda
	ZZ pcnst = pmult(cnst); // p * (1 - gamma * lambda)
	for (long i = 0; i < factorDim; ++i) {
		scheme.multByConstAndEqual(cwData[i], pcnst); // p * p * (1 - gamma * lambda) * w
		scheme.modSwitchOneAndEqual(cwData[i]); // p * (1 - gamma * lambda) * w  (-1)
	}

	debugcheck("c (1-gamma * lambda) * wdata: ", secretKey, cwData, dimcheck, slotscheck);

	ZZ minusp = -scheme.params.p; // -p
	scheme.addConstAndEqual(cip, minusp); // p * (ip - 1) (-1)
	scheme.doubleAndEqual(cip); // p * (2 * ip - 2) (-1)

	debugcheck("c (2ip - 2): ", secretKey, cip, slotscheck);

	Cipher* cgrad = new Cipher[factorDim];

	cnst = to_RR(gamma) / learnDim; // gamma / learnDim
	pcnst = pmult(cnst); // p * gamma / learnDim

	NTL_EXEC_RANGE(factorDim, first, last);
	for (long i = first; i < last; ++i) {
		cgrad[i] = scheme.multByConst(cxyData[i], pcnst); // p * p * gamma * (z / m)
		scheme.modSwitchOneAndEqual(cgrad[i]); // p * gamma * (z / m) (-1)
		scheme.modEmbedAndEqual(cgrad[i], cip.level);
		scheme.multModSwitchOneAndEqual(cgrad[i], cip); // p * gamma * (z * ip / m) (-2)
	}
	NTL_EXEC_RANGE_END;

	debugcheck("c grad: ", secretKey, cgrad, dimcheck, slotscheck);

	long logslots = log2(slots);
	long logwnum = log2(wBatch);

	NTL_EXEC_RANGE(factorDim, first, last);
	for (long i = first; i < last; ++i) {
		for (long l = logwnum; l < logslots; ++l) {
			Cipher rot = scheme.leftRotateByPo2(cgrad[i], l);
			scheme.addAndEqual(cgrad[i], rot); // p * gamma * sum(z * (ip - 2)) (-2)
		}
		scheme.modEmbedAndEqual(cwData[i], cgrad[i].level);
		scheme.subAndEqual(cwData[i], cgrad[i]); // w - gamma * (sum(z * (ip - 2)) - lambda * w) (-2)
	}
	NTL_EXEC_RANGE_END;

	debugcheck("c grad: ", secretKey, cgrad, dimcheck, slotscheck);

	debugcheck("c wdata: ", secretKey, cwData, dimcheck, slotscheck);
}

void CipherGD::encStepLGD(Cipher*& cxyData, Cipher*& cwData, long& slots, long& factorDim, long& learnDim, long& wBatch, double& lambda, double& gamma) {

	long dimcheck = 5;
	long slotscheck = 10;
	debugcheck("c xyData: ", secretKey, cxyData, dimcheck, slotscheck);

	debugcheck("c wData: ", secretKey, cwData, dimcheck, slotscheck);

	Cipher* cprod = new Cipher[factorDim];

	NTL_EXEC_RANGE(factorDim, first, last);
	for (long i = first; i < last; ++i) {
		cprod[i] = scheme.modEmbed(cxyData[i], cwData[i].level);
		scheme.multAndEqual(cprod[i], cwData[i]);
	}
	NTL_EXEC_RANGE_END;

	Cipher cip = algo.sum(cprod, factorDim);

	scheme.modSwitchOneAndEqual(cip); // cip (-1)

	debugcheck("c inner prod:", secretKey, cip, slotscheck);


	RR cnst = 1.0 - to_RR(gamma) * lambda; // 1 - gamma * lambda
	ZZ pcnst = pmult(cnst); // p * (1 - gamma * lambda)

	for (long i = 0; i < factorDim; ++i) {
		scheme.multByConstAndEqual(cwData[i], pcnst); // (1 - gamma * lambda) * w
		scheme.modSwitchOneAndEqual(cwData[i]); // (1 - gamma * lambda) * w  (-1)
	}

	debugcheck("c (1-gamma * lambda) * wdata: ", secretKey, cwData, dimcheck, slotscheck);

	double* coeffs = scheme.aux.taylorCoeffsMap.at(SIGMOIDPRIMEGOOD);
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

	debugcheck("c grad: ", secretKey, cgrad, dimcheck, 20);

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

	debugcheck("c grad: ", secretKey, cgrad, dimcheck, 20);

	debugcheck("c wData: ", secretKey, cwData, dimcheck, slotscheck);
}

void CipherGD::encStepMLGD(Cipher*& cxyData, Cipher*& cwData, Cipher*& cvData, long& slots, long& factorDim, long& learnDim, long& wBatch, double& lambda, double& gamma, double& eta) {

	long dimcheck = 5;
	long slotscheck = 10;
	debugcheck("c xyData: ", secretKey, cxyData, dimcheck, slotscheck);

	debugcheck("c wData: ", secretKey, cwData, dimcheck, slotscheck);

	Cipher* cprod = new Cipher[factorDim];

	NTL_EXEC_RANGE(factorDim, first, last);
	for (long i = first; i < last; ++i) {
		cprod[i] = scheme.modEmbed(cxyData[i], cwData[i].level);
		scheme.multAndEqual(cprod[i], cwData[i]);
	}
	NTL_EXEC_RANGE_END;

	Cipher cip = algo.sum(cprod, factorDim);

	scheme.modSwitchOneAndEqual(cip); // cip (-1)

	debugcheck("c inner prod:", secretKey, cip, slotscheck);


	RR cnst = 1.0 - to_RR(gamma) * lambda; // 1 - gamma * lambda
	ZZ pcnst = pmult(cnst); // p * (1 - gamma * lambda)

	for (long i = 0; i < factorDim; ++i) {
		scheme.multByConstAndEqual(cwData[i], pcnst); // (1 - gamma * lambda) * w
		scheme.modSwitchOneAndEqual(cwData[i]); // (1 - gamma * lambda) * w  (-1)
	}

	debugcheck("c (1-gamma * lambda) * wdata: ", secretKey, cwData, dimcheck, slotscheck);

	double* coeffs = scheme.aux.taylorCoeffsMap.at(SIGMOIDPRIMEGOOD);
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

	debugcheck("c grad: ", secretKey, cgrad, dimcheck, 20);

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

	debugcheck("c grad: ", secretKey, cgrad, dimcheck, 20);

	debugcheck("c wData: ", secretKey, cwData, dimcheck, slotscheck);
	//TODO: kimandrik
}

void CipherGD::encStepNLGD(Cipher*& cxyData, Cipher*& cwData, Cipher*& cvData, long& slots, long& factorDim, long& learnDim, long& wBatch, double& lambda, double& gamma, double& eta) {

	Cipher* cprod = new Cipher[factorDim];

	NTL_EXEC_RANGE(factorDim, first, last);
	for (long i = first; i < last; ++i) {
		cprod[i] = scheme.modEmbed(cxyData[i], cwData[i].level);
		scheme.multAndEqual(cprod[i], cwData[i]);
	}
	NTL_EXEC_RANGE_END;

	debugcheck("c prod: ", secretKey, cprod, 5, 7);

	Cipher cip = algo.sum(cprod, factorDim);


	scheme.modSwitchOneAndEqual(cip); // cip (-1)
	debugcheck("c ip: ", secretKey, cip, 7);

	double* coeffs = scheme.aux.taylorCoeffsMap.at(SIGMOIDPRIMEGOOD);
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

	debugcheck("c grad: ", secretKey, cgrad, 5, 7);

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

	debugcheck("c grad: ", secretKey, cgrad, 5, 7);

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

	debugcheck("c grad: ", secretKey, cgrad, 5, 7);

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
	debugcheck("c wData: ", secretKey, cwData, 5, 1);
	debugcheck("c vData: ", secretKey, cvData, 5, 1);
}

Cipher* CipherGD::encwsum(Cipher*& cwData, long& factorDim, long& wBatch) {
	Cipher* cw = new Cipher[factorDim];
	for (long i = 0; i < factorDim; ++i) {
		cw[i] = algo.partialSlotsSum(cwData[i], wBatch);
	}
	return cw;
}

double* CipherGD::decw(SecKey& secretKey, Cipher*& cw, long& factorDim) {
	double* w = new double[factorDim];
	for (long i = 0; i < factorDim; ++i) {
		CZZ* dcw = scheme.decrypt(secretKey, cw[i]);
		RR wi = to_RR(dcw[0].r);
		wi.e -= scheme.params.logp;
		w[i] = to_double(wi);
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
