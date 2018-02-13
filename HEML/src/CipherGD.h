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

	void encZData(Ciphertext* encZData, double** zData, long slots, long factorDim, long sampleDim, long batch, long cnum, long wBits);

	void encWDataVDataAverage(Ciphertext* encWData, Ciphertext* encVData, Ciphertext* encZData, long cnum, long sBits, long bBits);

	void encWDataVDataZero(Ciphertext* encWData, Ciphertext* encVData, long cnum, long slots, long wBits);

	ZZX generateAuxPoly(long slots, long batch, long pBits);

	Ciphertext encInnerProduct(Ciphertext* encZData, Ciphertext* encWData, ZZX& poly, long cnum, long bBits, long wBits, long pBits);

	void encSigmoid(long kdeg, Ciphertext* encZData, Ciphertext* encGrad, Ciphertext& encIP, long cnum, double gamma, long sBits, long bBits, long wBits, long aBits);

	void encLGDstep(Ciphertext* encWData, Ciphertext* encGrad, long cnum);
	void encMLGDstep(Ciphertext* encWData, Ciphertext* encVData, Ciphertext* encGrad, double eta, long cnum, long pBits);
	void encNLGDstep(Ciphertext* encWData, Ciphertext* encVData, Ciphertext* encGrad, double eta, long cnum, long pBits);

	void encLGDiteration(long kdeg, Ciphertext* encZData, Ciphertext* encWData, ZZX& poly, long cnum, double gamma, long sBits, long bBits, long wBits, long pBits, long aBits);
	void encMLGDiteration(long kdeg, Ciphertext* encZData, Ciphertext* encWData, Ciphertext* encVData, ZZX& poly, long cnum, double gamma, double eta, long sBits, long bBits, long wBits, long pBits, long aBits);
	void encNLGDiteration(long kdeg, Ciphertext* encZData, Ciphertext* encWData, Ciphertext* encVData, ZZX& poly, long cnum, double gamma, double eta, long sBits, long bBits, long wBits, long pBits, long aBits);

	void decWData(double* wData, Ciphertext* encWData, long factorDim, long batch, long cnum, long wBits);
};

#endif /* CIPHERGD_H_ */
