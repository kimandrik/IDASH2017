#ifndef HEML_CIPHERGD_H_
#define HEML_CIPHERGD_H_

#include <Scheme.h>
#include <SecretKey.h>

#include <complex>

using namespace std;
using namespace NTL;

class CipherGD {
public:
	Scheme& scheme;
	SecretKey& secretKey;

	CipherGD(Scheme& scheme, SecretKey& secretKey) : scheme(scheme), secretKey(secretKey) {}

	void encxyData(Ciphertext* cxyData, double** xyData, long slots, long factorDim, long sampleDim, long batch, long cnum, long wBits);

	void encwData(Ciphertext* cwData, Ciphertext* cxyData, long cnum, long sBits, long bBits);

	void encwData0(Ciphertext* cwData, long cnum, long slots, long wBits);

	ZZX generateAuxPoly(long slots, long batch, long pBits);

	Ciphertext encIP(Ciphertext* cxyData, Ciphertext* cwData, ZZX& poly, long cnum, long bBits, long wBits, long pBits);

	void encSigmoid(long approxDeg, Ciphertext* cxyData, Ciphertext* cgrad, Ciphertext& cip, long cnum, double gamma, long sBits, long bBits, long wBits, long aBits);

	void encLGDstep(Ciphertext* cwData, Ciphertext* cgrad, long cnum);
	void encMLGDstep(Ciphertext* cwData, Ciphertext* cvData, Ciphertext* cgrad, double eta, long cnum, long pBits);
	void encNLGDstep(Ciphertext* cwData, Ciphertext* cvData, Ciphertext* cgrad, double eta, long cnum, long pBits);

	void encLGDiteration(long approxDeg, Ciphertext* cxyData, Ciphertext* cwData, ZZX& poly, long cnum, double gamma, long sBits, long bBits, long wBits, long pBits, long aBits);
	void encMLGDiteration(long approxDeg, Ciphertext* cxyData, Ciphertext* cwData, Ciphertext* cvData, ZZX& poly, long cnum, double gamma, double eta, long sBits, long bBits, long wBits, long pBits, long aBits);
	void encNLGDiteration(long approxDeg, Ciphertext* cxyData, Ciphertext* cwData, Ciphertext* cvData, ZZX& poly, long cnum, double gamma, double eta, long sBits, long bBits, long wBits, long pBits, long aBits);

	void decwData(double* wData, Ciphertext* cwData, long factorDim, long batch, long cnum, long wBits);
};

#endif /* CIPHERGD_H_ */
