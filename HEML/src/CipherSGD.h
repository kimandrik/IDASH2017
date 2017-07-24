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
	SecKey secretKey;

	CipherSGD(Scheme& scheme, SchemeAlgo& algo, SecKey& secretKey) : scheme(scheme), algo(algo), secretKey(secretKey) {}

	Cipher* enczdata(long**& zdata, long& slots, long& wnum, long& dim, long& learndim, ZZ& p);
	Cipher* encwdata(double**& wdata, long& slots, long& wnum, long& dim, long& learndim, long& logp);

	void encStepQuadraticRegress(Cipher*& czdata, Cipher*& cwdata, ZZ& pgamma, double& lambda, long& slots, long& wnum, long& dim, long& learndim);
	void encStepLogRegress(Cipher*& czdata, Cipher*& cwdata, ZZ& pgamma, double& lambda, long& slots, long& wnum, long& dim, long& learndim);

	Cipher* encwout(Cipher*& cwdata, long& wnum, long& dim);
	double* decw(SecKey& secretKey, Cipher*& cw, long& dim);

	void debugcheck(string prefix, SecKey& secretKey, Cipher*& ciphers, long dim, long slots);
	void debugcheck(string prefix, SecKey& secretKey, Cipher& cipher, long slots);

};

#endif /* CIPHERSGD_H_ */
