#ifndef CIPHERGD_H_
#define CIPHERGD_H_

#include <Scheme.h>
#include <SecKey.h>

using namespace NTL;

class CipherGD {
public:
	Scheme& scheme;
	SecKey& secretKey;

	double degree3[4] = {-0.5,0.19,0,-0.0035};
	double degree5[6] = {-0.5,0.2166,0,-0.0077,0,0.00011};

	CipherGD(Scheme& scheme, SecKey& secretKey) : scheme(scheme), secretKey(secretKey) {}

	void encxyData(Cipher*& cxyData, long**& xyData, long& slots, long& factorDim, long& learnDim, long& batch, long& cnum, long& xyBits);
	void encwData(Cipher*& cwData, Cipher*& cxyData, long& cnum, long& sBits, long& ldimBits, long& bBits, long& xyBits, long& wBits);

	void encStepNLGD3(Cipher*& cxyData, Cipher*& cwData, Cipher*& cvData, ZZX& poly, long& cnum, double& gamma, double& eta, double& etaprev, long& sBits, long& bBits, long& xyBits, long& wBits, long& pBits, long& aBits);
	void encStepNLGD5(Cipher*& cxyData, Cipher*& cwData, Cipher*& cvData, ZZX& poly, long& cnum, double& gamma, double& eta, double& etaprev, long& sBits, long& bBits, long& xyBits, long& wBits, long& pBits, long& aBits);

	void decwData(double*& wData, Cipher*& cwData, long& factorDim, long& batch, long& cnum, long& wBits);
};

#endif /* CIPHERGD_H_ */
