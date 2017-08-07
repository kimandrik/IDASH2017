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

	Cipher* encxyDataWB(long**& xyData, long& slots, long& factorDim, long& learnDim, long& wBatch);
	Cipher* encwDataWB(double**& wData, long& slots, long& factorDim, long& learnDim, long& wBatch);

	Cipher* encxyDataXYB(long**& xyData, long& slots, long& factorDim, long& learnDim, long& learnDimPo2, long& xyBatch, long& cnum);
	Cipher* encwDataXYB(double*& wData, long& slots, long& factorDim, long& learnDim, long& learnDimPo2, long& xyBatch, long& cnum);

	void encwsumWB(Cipher*& cwData, long& factorDim, long& wBatch);

	void encStepLGDWB(Cipher*& cxyData, Cipher*& cwData, long& slots, long& factorDim, long& learnDim, long& wBatch, double& lambda, double& gamma);

	void encStepNLGD7WB(Cipher*& cxyData, Cipher*& cwData, Cipher*& cvData, long& slots, long& factorDim, long& learnDim, long& wBatch, double& gamma, double& eta);

	void encStepNLGD7XYB6(Cipher*& cxyData, Cipher*& cwData, Cipher*& cvData, ZZX& kkk, long& slots, long& learnDim, long learnDimPo2, long& xybatchBits, long& xyBatch, long& cnum, double& gamma, double& eta);
	void encStepNLGD3XYB5(Cipher*& cxyData, Cipher*& cwData, Cipher*& cvData, ZZX& kkk, long& slots, long& learnDim, long learnDimPo2, long& xybatchBits, long& xyBatch, long& cnum, double& gamma, double& eta);

	void encStepNLGD7XYBfast5(Cipher*& cxyData, Cipher*& cwData, Cipher*& cvData, ZZX& poly, long& slots, long& learnDim, long learnDimPo2, long& xybatchBits, long& xyBatch, long& cnum, double& gamma, double& eta, double& etaprev);
	void encStepNLGD3XYBfast4(Cipher*& cxyData, Cipher*& cwData, Cipher*& cvData, ZZX& poly, long& slots, long& learnDim, long learnDimPo2, long& xybatchBits, long& xyBatch, long& cnum, double& gamma, double& eta, double& etaprev);

	double* decWB(SecKey& secretKey, Cipher*& cw, long& factorDim);
	double* decXYB(SecKey& secretKey, Cipher*& cw, long& factorDim, long& xyBatch, long& cnum);

	ZZ pmult(RR val);

	void debugcheck(string prefix, SecKey& secretKey, Cipher*& ciphers, long factorCheck, long slotCheck);
	void debugcheck(string prefix, SecKey& secretKey, Cipher& cipher, long slotCheck);

};

#endif /* CIPHERGD_H_ */
