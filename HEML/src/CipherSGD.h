#ifndef CIPHERSGD_H_
#define CIPHERSGD_H_

#include <Scheme.h>
#include <SchemeAlgo.h>
#include <SecKey.h>

using namespace NTL;

class CipherSGD {
public:
	Scheme scheme;
	SchemeAlgo algo;

	CipherSGD(Scheme& scheme, SchemeAlgo& algo) : scheme(scheme), algo(algo) {}

	Cipher* enczdata(long**& zdata, long& slots, long& wnum, long& dim, long& sampledim, ZZ& p);
	Cipher* encwdata(double**& wdata, long& slots, long& wnum, long& dim, long& sampledim, long& logp);

	ZZ* pgammagen(double*& alpha, long& iter, long& logp);

	void encSteplogregress(Cipher*& czdata, Cipher*& cwdata, ZZ& pgamma, double& lambda, long& slots, long& wnum, long& dim);

	void encStepsimpleregress(Cipher*& czdata, Cipher*& cwdata, ZZ& pgamma, double& lambda, long& slots, long& wnum, long& dim);

	Cipher* encwout(Cipher*& cwdata, long& wnum, long& dim);
	double* decw(SecKey& secretKey, Cipher*& cw, long& dim);

};

#endif /* CIPHERSGD_H_ */
