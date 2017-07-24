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

ZZ* CipherSGD::pgammagen(double*& gamma, long& iter, long& logp) {
	ZZ* pgamma = new ZZ[iter];
	for (long k = 0; k < iter; ++k) {
		RR rgamma = to_RR(gamma[k]);
		RR prgamma = MakeRR(rgamma.x, rgamma.e + logp);
		pgamma[k] = to_ZZ(prgamma);
	}
	return pgamma;
}

void CipherSGD::encSteplogregress(Cipher*& czdata, Cipher*& cwdata, ZZ& pgamma, double& lambda, long& slots, long& wnum, long& dim, long& sampledim) {

	long dimcheck = 5;
	debugcheck("c zdata: ", secretKey, czdata, dimcheck);

	debugcheck("c wdata: ", secretKey, cwdata, dimcheck);

	Cipher* cprod = new Cipher[dim];

	NTL_EXEC_RANGE(dim, first, last);
	for (long i = first; i < last; ++i) {
		cprod[i] = scheme.modEmbed(czdata[i], cwdata[i].level);
		scheme.multAndEqual(cprod[i], cwdata[i]);
	}
	NTL_EXEC_RANGE_END;

	Cipher cip = algo.sum(cprod, dim);

	scheme.modSwitchOneAndEqual(cip); // cip (-1)

	debugcheck("c inner prod:", secretKey, cip);

	cout << "pgamma=" << pgamma << endl;
	cout << "lamda=" << lambda << endl;
	ZZ wcnst = scheme.params.p - to_ZZ(to_RR(pgamma) * lambda);

	cout << "p(1 - lambda * gamma)=" << wcnst << endl;

	for (long i = 0; i < dim; ++i) {
		scheme.multByConstAndEqual(cwdata[i], wcnst); // (1 - gamma * lambda) * w
		scheme.modSwitchOneAndEqual(cwdata[i]); // (1 - gamma * lambda) * w  (-1)
	}

	debugcheck("c (1-gamma * lambda) * wdata: ", secretKey, cwdata, dimcheck);

	ZZ* pows = scheme.aux.taylorPowsMap.at(SIGMOIDPRIMEGOOD);
	Cipher* cpows = algo.powerOf2Extended(cip, 2); // ip (-1), ip^2 (-2), ip^4 (-3)

	Cipher* cgrad = new Cipher[dim];

	ZZ cnst = pgamma * pows[0] / scheme.params.p; // p * gamma * (alpha_t)

	NTL_EXEC_RANGE(dim, first, last);
	for (long i = first; i < last; ++i) {
		cgrad[i] = scheme.multByConst(czdata[i], cnst); // p * p * gamma * (z * alpha_t)
		scheme.modSwitchOneAndEqual(cgrad[i]); // p * gamma * (z * alpha_t) (-1)
		scheme.modEmbedAndEqual(cgrad[i], cpows[0].level);
		scheme.multModSwitchOneAndEqual(cgrad[i], cpows[0]); // p * gamma * (z * alpha_t * ip^t)
	}
	NTL_EXEC_RANGE_END;

	for (long t = 1; t < 8; t=t+2) {
		cnst = pgamma * pows[t] / scheme.params.p; // p * gamma * (alpha_t)
		NTL_EXEC_RANGE(dim, first, last);
		for (long i = first; i < last; ++i) {
			Cipher cgradit = scheme.multByConst(czdata[i], cnst); // p * p * gamma * (z * alpha_t)
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

	debugcheck("c grad: ", secretKey, cgrad, dimcheck);

	long logslots = log2(slots);
	long logwnum = log2(wnum);

	NTL_EXEC_RANGE(dim, first, last);
	for (long i = first; i < last; ++i) {
		for (long l = logwnum; l < logslots; ++l) {
			Cipher rot = scheme.leftRotateByPo2(cgrad[i], l);
			scheme.addAndEqual(cgrad[i], rot); // p * gamma * sum(z * sigmoid(ip)) (-4)
		}
		scheme.modEmbedAndEqual(cwdata[i], cgrad[i].level);
		scheme.subAndEqual(cwdata[i], cgrad[i]); // w - gamma * grad(w) - gamma * lambda * w (-4)
	}
	NTL_EXEC_RANGE_END;

	debugcheck("c grad: ", secretKey, cgrad, dimcheck);

	debugcheck("c wdata: ", secretKey, cwdata, dimcheck);
}

void CipherSGD::encStepsimpleregress(Cipher*& czdata, Cipher*& cwdata, ZZ& pgamma, double& lambda, long& slots, long& wnum, long& dim, long& learndim) {

	long dimcheck = 10;
	debugcheck("c zdata: ", secretKey, czdata, dimcheck);

	debugcheck("c wdata: ", secretKey, cwdata, dimcheck);

	Cipher* cprod = new Cipher[dim];
	NTL_EXEC_RANGE(dim, first, last);
	for (long i = first; i < last; ++i) {
		cout << i << endl;
		cprod[i] = scheme.modEmbed(czdata[i], cwdata[i].level);
		scheme.multAndEqual(cprod[i], cwdata[i]);
	}
	NTL_EXEC_RANGE_END;

	Cipher cip = algo.sum(cprod, dim);

	scheme.modSwitchOneAndEqual(cip);

	debugcheck("c inner prod:", secretKey, cip);

	ZZ wcnst = scheme.params.p - to_ZZ(to_RR(pgamma) * lambda); // p * (1 - gamma * lambda)
	for (long i = 0; i < dim; ++i) {
		scheme.multByConstAndEqual(cwdata[i], wcnst); // p * p * (1 - gamma * lambda) * w
		scheme.modSwitchOneAndEqual(cwdata[i]); // p * (1 - gamma * lambda) * w  (-1)
	}

	debugcheck("c (1-gamma * lambda) * wdata: ", secretKey, cwdata, dimcheck);

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

	debugcheck("c grad: ", secretKey, cgrad, dimcheck);

	long logslots = log2(slots);
	long logwnum = log2(wnum);

	NTL_EXEC_RANGE(dim, first, last);
	for (long i = first; i < last; ++i) {
		for (long l = logwnum; l < logslots; ++l) {
			Cipher rot = scheme.leftRotateByPo2(cgrad[i], l);
			scheme.addAndEqual(cgrad[i], rot); // p * gamma * sum(z * (ip - 2)) (-2)
		}
		scheme.modEmbedAndEqual(cwdata[i], cgrad[i].level);
		scheme.subAndEqual(cwdata[i], cgrad[i]); // w - gamma * (sum(z * (ip - 2)) - lambda * w) (-2)
	}
	NTL_EXEC_RANGE_END;

	debugcheck("c grad: ", secretKey, cgrad, dimcheck);

	debugcheck("c wdata: ", secretKey, cwdata, dimcheck);
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
	cout << prefix << " " << ciphers[0].level << endl;
	for (long i = 0; i < dim; ++i) {
//		cout << i << endl;
//		cout << "start dec" << endl;
		CZZ* deci = scheme.decrypt(secretKey, ciphers[i]);
//		cout << "stop dec" << endl;
		RR wi = to_RR(deci[0].r);
		wi.e -= scheme.params.logp;
		double w = to_double(wi);
		cout << w << ",";
	}
	cout << endl;
}

void CipherSGD::debugcheck(string prefix, SecKey& secretKey, Cipher& cipher) {
	cout << prefix << " " << cipher.level << endl;
	CZZ* dec = scheme.decrypt(secretKey, cipher);
	RR wrr = to_RR(dec[0].r);
	wrr.e -= scheme.params.logp;
	double w = to_double(wrr);
	cout << w << endl;
}
