#ifndef CIPHERGD_H_
#define CIPHERGD_H_

#include <Scheme.h>
#include <SchemeAlgo.h>
#include <SecKey.h>

using namespace NTL;

class CipherGD {
public:
	Scheme scheme;
	SchemeAlgo algo;
	SecKey secretKey;

	CipherGD(Scheme& scheme, SchemeAlgo& algo, SecKey& secretKey) : scheme(scheme), algo(algo), secretKey(secretKey) {}

	Cipher* encxyData(long**& xyData, long& slots, long& factorDim, long& learnDim, long& wBatch);
	Cipher* encwData(double**& wData, long& slots, long& factorDim, long& learnDim, long& wBatch);

	void encStepQGD(Cipher*& cxyData, Cipher*& cwData, long& slots, long& factorDim, long& learnDim, long& wBatch, double& lambda, double& gamma);
	void encStepLGD(Cipher*& cxyData, Cipher*& cwData, long& slots, long& factorDim, long& learnDim, long& wBatch, double& lambda, double& gamma);

	void encStepMLGD(Cipher*& cxyData, Cipher*& cwData, Cipher*& cvData, long& slots, long& factorDim, long& learnDim, long& wBatch, double& lambda, double& gamma, double& eta);

	void encStepNLGD(Cipher*& cxyData, Cipher*& cwData, Cipher*& cvData, long& slots, long& factorDim, long& learnDim, long& wBatch, double& gamma, double& eta);

	Cipher* encwsum(Cipher*& cwData, long& factorDim, long& wBatch);
	double* decw(SecKey& secretKey, Cipher*& cw, long& factorDim);

	ZZ pmult(RR val);

	void debugcheck(string prefix, SecKey& secretKey, Cipher*& ciphers, long factorCheck, long slotCheck);
	void debugcheck(string prefix, SecKey& secretKey, Cipher& cipher, long slotCheck);

};

#endif /* CIPHERGD_H_ */
