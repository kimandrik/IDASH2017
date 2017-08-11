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

	Cipher* encxyData(long**& xyData, long& slots, long& factorDim, long& learnDim, long& learnDimPo2, long& xyBatch, long& cnum, long& xyBits);
	Cipher* encwData(Cipher*& cxyData, long& slotBits, long& ldimBits, long& xyBatchBits, long& cnum, long& xyBits, long& wBits);

	void encStepNLGD5(Cipher*& cxyData, Cipher*& cwData, Cipher*& cvData, ZZX& poly, long& slots, long& learnDim, long learnDimPo2, long& xybatchBits, long& xyBatch, long& cnum, double& gamma, double& eta, double& etaprev, long& xyBits, long& wBits, long& pBits);

	void encStepNLGD3(Cipher*& cxyData, Cipher*& cwData, Cipher*& cvData, ZZX& poly, long& slots, long& learnDim, long learnDimPo2, long& xybatchBits, long& xyBatch, long& cnum, double& gamma, double& eta, double& etaprev, long& xyBits, long& wBits, long& pBits);

	double* decwData(SecKey& secretKey, Cipher*& cw, long& factorDim, long& xyBatch, long& cnum, long& wBits);

	void debugcheck(string prefix, SecKey& secretKey, Cipher*& ciphers, long factorCheck, long slotCheck);
	void debugcheck(string prefix, SecKey& secretKey, Cipher& cipher, long slotCheck);

};

#endif /* CIPHERGD_H_ */
