#ifndef CIPHERGD_H_
#define CIPHERGD_H_

#include <Scheme.h>
#include <SecKey.h>

using namespace NTL;

class CipherGD {
public:
	Scheme& scheme;
	SecKey& secretKey;

	double degree3[3] = {-0.5,0.19,-0.0035};
	double degree5[4] = {-0.5,0.2166,-0.0077,0.00011};
	double degree7[5] = {-0.5,0.216884,-0.00819276,0.000165861,-0.00000119581};

	CipherGD(Scheme& scheme, SecKey& secretKey) : scheme(scheme), secretKey(secretKey) {}

	void encxyData(Cipher*& cxyData, long**& xyData, long& slots, long& factorDim, long& learnDim, long& batch, long& cnum, long& xyBits);
	void encwData(Cipher*& cwData, Cipher*& cxyData, long& cnum, long& sBits, long& ldimBits, long& bBits, long& xyBits, long& wBits);
	ZZX generateAuxPoly(long& slots, long& batch, long& pBits);

	Cipher encIP(Cipher*& cxyData, Cipher*& cwData, Cipher*& cgrad, Cipher*& cprod, ZZX& poly, long& cnum, long& bBits, long& xyBits, long& pBits, long& aBits);

	void encSigmoid(long& approxDeg, Cipher*& cxyData, Cipher*& cgrad, Cipher*& cprod, Cipher& cip, long& cnum, double& gamma, long& sBits, long& bBits, long& xyBits, long& wBits, long& pBits, long& aBits);

	void encLGDstep(Cipher*& cwData, Cipher*& cgrad, long& cnum, long& wBits, long& bitsDown);
	void encMLGDstep(Cipher*& cwData, Cipher*& cvData, Cipher*& cgrad, double& eta, long& cnum, long& wBits, long& bitsDown);
	void encNLGDstep(Cipher*& cwData, Cipher*& cvData, Cipher*& cgrad, double& eta, double& etaprev, long& cnum, long& wBits, long& bitsDown);

	void encLGDiteration(long& approxDeg, Cipher*& cxyData, Cipher*& cwData, ZZX& poly, long& cnum, double& gamma, long& sBits, long& bBits, long& xyBits, long& wBits, long& pBits, long& aBits);
	void encMLGDiteration(long& approxDeg, Cipher*& cxyData, Cipher*& cwData, Cipher*& cvData, ZZX& poly, long& cnum, double& gamma, double& eta, long& sBits, long& bBits, long& xyBits, long& wBits, long& pBits, long& aBits);
	void encNLGDiteration(long& approxDeg, Cipher*& cxyData, Cipher*& cwData, Cipher*& cvData, ZZX& poly, long& cnum, double& gamma, double& eta, double& etaprev, long& sBits, long& bBits, long& xyBits, long& wBits, long& pBits, long& aBits);

	void decwData(double*& wData, Cipher*& cwData, long& factorDim, long& batch, long& cnum, long& wBits);
};

#endif /* CIPHERGD_H_ */
