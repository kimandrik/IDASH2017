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

	Cipher* encxyData(long**& xyData, long& slots, long& wnum, long& factorDim, long& learnDim);
	Cipher* encwData(double**& wData, long& slots, long& wnum, long& factorDim, long& learnDim);

	void encStepQGD(Cipher*& cxyData, Cipher*& cwData, ZZ& pgamma, double& lambda, long& slots, long& wBatch, long& factorDim, long& learnDim);
	void encStepLGD(Cipher*& cxyData, Cipher*& cwData, ZZ& pgamma, double& lambda, long& slots, long& wBatch, long& factorDim, long& learnDim);

	void encStepMLGD(Cipher*& cxyData, Cipher*& cwData, ZZ& pgamma, double& lambda, double& eta, long& slots, long& wBatch, long& factorDim, long& learnDim);
	void encStepNLGD(Cipher*& cxyData, Cipher*& cwData, ZZ& pgamma, double& lambda, double& eta, long& slots, long& wBatch, long& factorDim, long& learnDim);

	Cipher* encwaverage(Cipher*& cwData, long& wBatch, long& factorDim);
	double* decw(SecKey& secretKey, Cipher*& cw, long& factorDim);

	void debugcheck(string prefix, SecKey& secretKey, Cipher*& ciphers, long factorData, long slots);
	void debugcheck(string prefix, SecKey& secretKey, Cipher& cipher, long slots);

};

#endif /* CIPHERSGD_H_ */
