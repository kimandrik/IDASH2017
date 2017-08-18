#ifndef CIPHERGD_H_
#define CIPHERGD_H_

#include <Scheme.h>
#include <SchemeAlgo.h>
#include <SecKey.h>

using namespace NTL;

class CipherGD {
public:
	Scheme& scheme;
	SecKey& secretKey;

	CipherGD(Scheme& scheme, SecKey& secretKey) : scheme(scheme), secretKey(secretKey) {}

	void encxyData(Cipher*& cxyData, long**& xyData, long& slots, long& factorDim, long& learnDim, long& xyBatch, long& cnum, long& xyBits);
	void encwData(Cipher*& cwData, Cipher*& cxyData, long& slotBits, long& ldimBits, long& xyBatchBits, long& cnum, long& xyBits, long& wBits);

	void encStepNLGD5(Cipher*& cxyData, Cipher*& cwData, Cipher*& cvData, ZZX& poly, long& slots, long& learnDim, long& xybatchBits, long& xyBatch, long& cnum, double& gamma, double& eta, double& etaprev, long& xyBits, long& wBits, long& pBits);
	void encStepNLGD3(Cipher*& cxyData, Cipher*& cwData, Cipher*& cvData, ZZX& poly, long& slots, long& learnDim, long& xybatchBits, long& xyBatch, long& cnum, double& gamma, double& eta, double& etaprev, long& xyBits, long& wBits, long& pBits);

	double* decwData(SecKey& secretKey, Cipher*& cw, long& factorDim, long& xyBatch, long& cnum, long& wBits);
};

#endif /* CIPHERGD_H_ */
