#include "CipherSGD.h"

#include <NTL/ZZ.h>
#include <NTL/BasicThreadPool.h>

#include <Cipher.h>
#include <CZZ.h>
#include <EvaluatorUtils.h>
#include <NTL/RR.h>
#include <NTL/ZZ.h>
#include <Params.h>
#include <SchemeAux.h>
#include <cmath>
#include <map>

Cipher* CipherSGD::encxyData(long**& xyData, long& slots, long& wBatch, long& factorDim, long& learnDim) {
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

Cipher* CipherSGD::encwData(double**& wData, long& slots, long& wBatch, long& factorDim, long& learnDim) {
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

void CipherSGD::encStepQuadraticRegress(Cipher*& cxyData, Cipher*& cwData, ZZ& pgamma, double& lambda, long& slots, long& wBatch, long& factorDim, long& learnDim) {

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

	ZZ wcnst = scheme.params.p - to_ZZ(to_RR(pgamma) * lambda); // p * (1 - gamma * lambda)
	for (long i = 0; i < factorDim; ++i) {
		scheme.multByConstAndEqual(cwData[i], wcnst); // p * p * (1 - gamma * lambda) * w
		scheme.modSwitchOneAndEqual(cwData[i]); // p * (1 - gamma * lambda) * w  (-1)
	}

	debugcheck("c (1-gamma * lambda) * wdata: ", secretKey, cwData, dimcheck, slotscheck);

	ZZ minusp = -scheme.params.p; // -p
	scheme.addConstAndEqual(cip, minusp); // p * (ip - 1) (-1)
	scheme.doubleAndEqual(cip); // p * (2 * ip - 2) (-1)

	debugcheck("c (2ip - 2): ", secretKey, cip, slotscheck);

	Cipher* cgrad = new Cipher[factorDim];
	ZZ cnst = pgamma / learnDim; // p * gamma * (1 / m)
	NTL_EXEC_RANGE(factorDim, first, last);
	for (long i = first; i < last; ++i) {
		cgrad[i] = scheme.multByConst(cxyData[i], cnst); // p * p * gamma * (z / m)
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

void CipherSGD::encStepLogRegress(Cipher*& cxyData, Cipher*& cwData, ZZ& pgamma, double& lambda, long& slots, long& wBatch, long& factorDim, long& learnDim) {

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

	cout << "pgamma=" << pgamma << endl;
	cout << "lamda=" << lambda << endl;
	ZZ wcnst = scheme.params.p - to_ZZ(to_RR(pgamma) * lambda);

	cout << "p(1 - lambda * gamma)=" << wcnst << endl;

	for (long i = 0; i < factorDim; ++i) {
		scheme.multByConstAndEqual(cwData[i], wcnst); // (1 - gamma * lambda) * w
		scheme.modSwitchOneAndEqual(cwData[i]); // (1 - gamma * lambda) * w  (-1)
	}

	debugcheck("c (1-gamma * lambda) * wdata: ", secretKey, cwData, dimcheck, slotscheck);

	ZZ* pows = scheme.aux.taylorPowsMap.at(SIGMOIDPRIMEGOOD);
	Cipher* cpows = algo.powerOf2Extended(cip, 2); // ip (-1), ip^2 (-2), ip^4 (-3)

	Cipher* cgrad = new Cipher[factorDim];

	ZZ cnst = pgamma * pows[0] / scheme.params.p; // p * gamma * (alpha_t)

	cout << "gamma * alpha0=" << cnst << endl;
	NTL_EXEC_RANGE(factorDim, first, last);
	for (long i = first; i < last; ++i) {
		cgrad[i] = scheme.multByConst(cxyData[i], cnst); // p * p * gamma * (z * alpha_t)
		scheme.modSwitchOneAndEqual(cgrad[i]); // p * gamma * (z * alpha_t) (-1)
		scheme.modEmbedAndEqual(cgrad[i], cpows[0].level);
		scheme.multModSwitchOneAndEqual(cgrad[i], cpows[0]); // p * gamma * (z * alpha_t * ip^t)
	}
	NTL_EXEC_RANGE_END;

	for (long t = 1; t < 8; t=t+2) {
		cnst = pgamma * pows[t] / scheme.params.p; // p * gamma * (alpha_t)
		cout << "gamma * alphat=" << cnst << endl;
		NTL_EXEC_RANGE(factorDim, first, last);
		for (long i = first; i < last; ++i) {
			Cipher cgradit = scheme.multByConst(cxyData[i], cnst); // p * p * gamma * (z * alpha_t)
			scheme.modSwitchOneAndEqual(cgradit); // p * gamma * (z * alpha_t) (-1)
			if(bit(t, 0)) {
				scheme.modEmbedAndEqual(cgradit, cpows[0].level);
				scheme.multModSwitchOneAndEqual(cgradit, cpows[0]); // p * gamma * (z * alpha_t * ip^t)
			}
			if(bit(t, 1)) {
				scheme.modEmbedAndEqual(cgradit, cpows[1].level);
				scheme.multModSwitchOneAndEqual(cgradit, cpows[1]); // p * gamma * (z * alpha_t * ip^t)
			}
			if(bit(t, 2)) {
				scheme.modEmbedAndEqual(cgradit, cpows[2].level);
				scheme.multModSwitchOneAndEqual(cgradit, cpows[2]); // p * gamma * (z * alpha_t * ip^t)
			}
			scheme.modEmbedAndEqual(cgrad[i], cgradit.level);
			scheme.addAndEqual(cgrad[i], cgradit); // p * gamma * (z * sigmoid(ip)) (-4)
		}
		NTL_EXEC_RANGE_END;
	}

	debugcheck("c grad: ", secretKey, cgrad, dimcheck, 20);

	long logslots = log2(slots);
	long logwnum = log2(wBatch);

	cout << "logslots=" << logslots << endl;
	cout << "logwnum=" << logwnum << endl;

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

Cipher* CipherSGD::encwaverage(Cipher*& cwData, long& wBatch, long& factorDim) {
	Cipher* cw = new Cipher[factorDim];
	for (long i = 0; i < factorDim; ++i) {
		cw[i] = algo.partialSlotsSum(cwData[i], wBatch);
	}
	return cw;
}

double* CipherSGD::decw(SecKey& secretKey, Cipher*& cw, long& factorDim) {
	double* w = new double[factorDim];
	for (long i = 0; i < factorDim; ++i) {
		CZZ* dcw = scheme.decrypt(secretKey, cw[i]);
		RR wi = to_RR(dcw[0].r);
		wi.e -= scheme.params.logp;
		w[i] = to_double(wi);
	}

	return w;
}

void CipherSGD::debugcheck(string prefix, SecKey& secretKey, Cipher*& ciphers, long factorDim, long slots) {
	cout << prefix << " " << ciphers[0].level << endl;
	for (long i = 0; i < factorDim; ++i) {
		CZZ* deci = scheme.decrypt(secretKey, ciphers[i]);
		for (int j = 0; j < factorDim; ++j) {
			RR wi = to_RR(deci[j].r);
			wi.e -= scheme.params.logp;
			double w = to_double(wi);
			cout << w << ",";
		}
		cout << endl;
	}
	cout << endl;
}

void CipherSGD::debugcheck(string prefix, SecKey& secretKey, Cipher& cipher, long slots) {
	cout << prefix << " " << cipher.level << endl;
	CZZ* dec = scheme.decrypt(secretKey, cipher);
	for (long j = 0; j < slots; ++j) {
		RR wrr = to_RR(dec[j].r);
		wrr.e -= scheme.params.logp;
		double w = to_double(wrr);
		cout << w << ",";
	}
	cout << endl;
}
