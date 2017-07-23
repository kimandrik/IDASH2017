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

Cipher* CipherSGD::enczdata(long**& zdata, long& slots, long& wnum, long& dim, long& learndim, ZZ& p) {
	Cipher* czdata = new Cipher[dim];
	NTL_EXEC_RANGE(dim, first, last);
	for (long i = first; i < last; ++i) {
		CZZ* pzdata = new CZZ[slots];
		for (long j = 0; j < learndim; ++j) {
			if(zdata[i][j] == -1) {
				for (long l = 0; l < wnum; ++l) {
					pzdata[wnum * j + l] = CZZ(-p);
				}
			} else if(zdata[i][j] == 1) {
				for (long l = 0; l < wnum; ++l) {
					pzdata[wnum * j + l] = CZZ(p);
				}
			}
		}
		czdata[i] = scheme.encrypt(pzdata, slots);
	}
	NTL_EXEC_INDEX_END;
	return czdata;
}

Cipher* CipherSGD::encwdata(double**& wdata, long& slots, long& wnum, long& dim, long& learndim, long& logp) {
	Cipher* cwdata = new Cipher[dim];
	NTL_EXEC_RANGE(dim, first, last);
	for (long i = first; i < last; ++i) {
		CZZ* pwdata = new CZZ[slots];
		for (long l = 0; l < wnum; ++l) {
			CZZ tmp = EvaluatorUtils::evaluateVal(wdata[l][i], 0.0, logp);
			for (long j = 0; j < learndim; ++j) {
				pwdata[wnum * j + l] = tmp;
			}
		}
		cwdata[i] = scheme.encrypt(pwdata, slots);
	}
	NTL_EXEC_RANGE_END;
	return cwdata;
}

ZZ* CipherSGD::pgammagen(double*& alpha, long& iter, long& logp) {
	ZZ* palpha = new ZZ[iter];
	for (long k = 0; k < iter; ++k) {
		RR ralpha = to_RR(alpha[k]);
		RR pralpha = MakeRR(ralpha.x, ralpha.e + logp);
		palpha[k] = to_ZZ(pralpha);
	}
	return palpha;
}

void CipherSGD::encSteplogregress(Cipher*& czdata, Cipher*& cwdata, ZZ& pgamma, double& lambda, long& slots, long& wnum, long& dim, long& sampledim) {

	Cipher cip = scheme.modEmbed(czdata[dim-1], cwdata[dim-1].level);
	scheme.multAndEqual(cip, cwdata[dim-1]);

	NTL_EXEC_RANGE(dim-1, first, last);
	for (long i = first; i < last; ++i) {
		Cipher cprodi = scheme.modEmbed(czdata[i], cwdata[i].level);
		scheme.multAndEqual(cprodi, cwdata[i]);
		scheme.addAndEqual(cip, cprodi);
	}
	NTL_EXEC_RANGE_END;

	scheme.modSwitchOneAndEqual(cip);

	Cipher* cpows = algo.powerOf2Extended(cip, 2); // ip (-1), ip^2 (-2), ip^4 (-3)

	ZZ* pows = scheme.aux.taylorPowsMap.at(SIGMOIDPRIMEGOOD);
	ZZ wcnst = scheme.params.p - to_ZZ(to_RR(pgamma) * lambda);
	for (long i = 0; i < dim; ++i) {
		scheme.multByConstAndEqual(cwdata[i], wcnst); // (1 - gamma * lambda) * w
		scheme.modSwitchOneAndEqual(cwdata[i]); // (1 - gamma * lambda) * w  (-1)
	}

	Cipher* cgrad = new Cipher[dim];
	for (long t = 0; t < 8; ++t) {
		ZZ cnst = pgamma * pows[t] / scheme.params.p; // p * gamma * (alpha_t)
		NTL_EXEC_RANGE(dim, first, last);
		for (long i = first; i < last; ++i) {
			if(cnst != ZZ::zero()) {
				Cipher cgradit = scheme.multByConst(czdata[i], cnst); // p * p * gamma * (z * alpha_t)
				scheme.modSwitchOneAndEqual(cgradit); // p * gamma * (z * alpha_t) (-1)
				for (int b = 0; b < 3; ++b) {
					if(bit(t, b)) {
						scheme.modEmbedAndEqual(cgradit, cpows[b].level);
						scheme.multModSwitchOneAndEqual(cgradit, cpows[b]); // p * gamma * (z * alpha_t * ip^t)
					}
				}
				scheme.modEmbedAndEqual(cgrad[i], cgradit.level);
				scheme.addAndEqual(cgrad[i], cgradit); // p * gamma * (z * sigmoid(ip)) (-4)
			}
		}
		NTL_EXEC_RANGE_END;
	}

	long logslots = log2(slots);
	long logwnum = log2(wnum);

	NTL_EXEC_RANGE(dim, first, last);
	for (long i = first; i < last; ++i) {
		for (long i = logwnum; i < logslots; ++i) {
			Cipher rot = scheme.leftRotateByPo2(cgrad[i], i);
			scheme.addAndEqual(cgrad[i], rot); // p * gamma * sum(z * sigmoid(ip)) (-4)
		}
		scheme.modEmbedAndEqual(cwdata[i], cgrad[i].level);
		scheme.subAndEqual(cwdata[i], cgrad[i]); // w - gamma * grad(w) - gamma * lambda * w (-4)
	}
	NTL_EXEC_RANGE_END;
}

void CipherSGD::encStepsimpleregress(Cipher*& czdata, Cipher*& cwdata, ZZ& pgamma, double& lambda, long& slots, long& wnum, long& dim, long& learndim) {

	debugcheck("c zdata: ", secretKey, czdata, dim);
	debugcheck("c wdata: ", secretKey, cwdata, dim);

	Cipher cip = scheme.modEmbed(czdata[dim-1], cwdata[dim-1].level);
	scheme.multAndEqual(cip, cwdata[dim-1]);

	NTL_EXEC_RANGE(dim-1, first, last);
	for (long i = first; i < last; ++i) {
		Cipher cprodi = scheme.modEmbed(czdata[i], cwdata[i].level);
		scheme.multAndEqual(cprodi, cwdata[i]);
		scheme.addAndEqual(cip, cprodi);
	}
	NTL_EXEC_RANGE_END;

	scheme.modSwitchOneAndEqual(cip);

	debugcheck("c inner prod:", secretKey, cip);
	ZZ wcnst = scheme.params.p - to_ZZ(to_RR(pgamma) * lambda); // p * (1 - gamma * lambda)
	for (long i = 0; i < dim; ++i) {
		scheme.multByConstAndEqual(cwdata[i], wcnst); // p * p * (1 - gamma * lambda) * w
		scheme.modSwitchOneAndEqual(cwdata[i]); // p * (1 - gamma * lambda) * w  (-1)
	}
	debugcheck("c (1-gamma * lambda) * wdata: ", secretKey, cwdata, dim);

	ZZ minusp = -scheme.params.p; // -p
	scheme.addConstAndEqual(cip, minusp); // p * (ip - 1) (-1)
	scheme.doubleAndEqual(cip); // p * (2 * ip - 2) (-1)

	debugcheck("c (2ip - 2): ", secretKey, cip);

	Cipher* cgrad = new Cipher[dim];
	ZZ cnst = pgamma / learndim; // p * gamma * (1 / m)
	NTL_EXEC_RANGE(dim, first, last);
	for (long i = first; i < last; ++i) {
		cgrad[i] = scheme.multByConst(czdata[i], cnst); // p * p * gamma * (z / m)
		scheme.modSwitchOneAndEqual(cgrad[i]); // p * gamma * (z / m) (-1)
		scheme.modEmbedAndEqual(cgrad[i], cip.level);
		scheme.multModSwitchOneAndEqual(cgrad[i], cip); // p * gamma * (z * ip / m) (-2)
	}
	NTL_EXEC_RANGE_END;

	debugcheck("c grad: ", secretKey, cgrad, dim);

	long logslots = log2(slots);
	long logwnum = log2(wnum);

	NTL_EXEC_RANGE(dim, first, last);
	for (long i = first; i < last; ++i) {
		for (long i = logwnum; i < logslots; ++i) {
			Cipher rot = scheme.leftRotateByPo2(cgrad[i], i);
			scheme.addAndEqual(cgrad[i], rot); // p * gamma * sum(z * (ip - 2)) (-2)
		}
		scheme.modEmbedAndEqual(cwdata[i], cgrad[i].level);
		scheme.subAndEqual(cwdata[i], cgrad[i]); // w - gamma * (sum(z * (ip - 2)) - lambda * w) (-2)
	}
	NTL_EXEC_RANGE_END;
	debugcheck("c grad: ", secretKey, cgrad, dim);
	debugcheck("c wdata: ", secretKey, cwdata, dim);
}

Cipher* CipherSGD::encwout(Cipher*& cwdata, long& wnum, long& dim) {
	Cipher* cw = new Cipher[dim];
	for (long i = 0; i < dim; ++i) {
		cw[i] = algo.partialSlotsSum(cwdata[i], wnum);
	}
	return cw;
}

double* CipherSGD::decw(SecKey& secretKey, Cipher*& cw, long& dim) {
	double* w = new double[dim];
	for (long i = 0; i < dim; ++i) {
		CZZ* dcw = scheme.decrypt(secretKey, cw[i]);
		RR wi = to_RR(dcw[0].r);
		wi.e -= scheme.params.logp;
		w[i] = to_double(wi);
	}

	return w;
}

void CipherSGD::debugcheck(string prefix, SecKey& secretKey, Cipher*& ciphers, long& dim) {
	cout << prefix;
	for (long i = 0; i < dim; ++i) {
		CZZ* deci = scheme.decrypt(secretKey, ciphers[i]);
		RR wi = to_RR(deci[0].r);
		wi.e -= scheme.params.logp;
		double w = to_double(wi);
		cout << w << ",";
	}
	cout << endl;
}

void CipherSGD::debugcheck(string prefix, SecKey& secretKey, Cipher& cipher) {
	cout << prefix;
	CZZ* deci = scheme.decrypt(secretKey, cipher);
	RR wi = to_RR(deci[0].r);
	wi.e -= scheme.params.logp;
	double w = to_double(wi);
	cout << w << "," << endl;
}
