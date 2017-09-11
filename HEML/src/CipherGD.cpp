#include "CipherGD.h"

#include <Cipher.h>
#include <CZZ.h>
#include <EvaluatorUtils.h>
#include <NTL/BasicThreadPool.h>
#include <NTL/RR.h>
#include <NTL/ZZ.h>

void CipherGD::encxyData(Cipher*& cxyData, long**& xyData, long& slots, long& factorDim, long& learnDim, long& batch, long& cnum, long& xyBits) {
	ZZ precision = power2_ZZ(xyBits);
	CZZ* pxyData = new CZZ[slots];
	for (long i = 0; i < cnum - 1; ++i) {
		for (long j = 0; j < learnDim; ++j) {
			for (long l = 0; l < batch; ++l) {
				pxyData[batch * j + l] = xyData[j][batch * i + l] == -1 ? CZZ(-precision) :
						xyData[j][batch * i + l] == 1 ? CZZ(precision) : CZZ();
			}
		}
		cxyData[i] = scheme.encrypt(pxyData, slots);
	}

	long rest = factorDim - batch * (cnum - 1);
	for (long j = 0; j < learnDim; ++j) {
		for (long l = 0; l < rest; ++l) {
			pxyData[batch * j + l] = xyData[j][batch * (cnum - 1) + l] == -1 ? CZZ(-precision) :
					xyData[j][batch * (cnum - 1) + l] == 1 ? CZZ(precision) : CZZ();
		}
		for (long l = rest; l < batch; ++l) {
			pxyData[batch * j + l] = CZZ();
		}
	}
	cxyData[cnum - 1] = scheme.encrypt(pxyData, slots);

	delete[] pxyData;
}

void CipherGD::encwData(Cipher*& cwData, Cipher*& cxyData, long& cnum, long& sBits, long& ldimBits, long& bBits, long& xyBits, long& wBits) {
	long downBits = ldimBits + xyBits - wBits;
	NTL_EXEC_RANGE(cnum, first, last);
	for (long i = first; i < last; ++i) {
		cwData[i] = cxyData[i];
		for (long l = bBits; l < sBits; ++l) {
			Cipher rot = scheme.leftRotateByPo2(cwData[i], l);
			scheme.addAndEqual(cwData[i], rot);
		}
		scheme.modSwitchAndEqual(cwData[i], downBits);
	}
	NTL_EXEC_RANGE_END;
}

void CipherGD::encStepLGD3(Cipher*& cxyData, Cipher*& cwData, ZZX& poly, long& cnum, double& gamma, long& sBits, long& bBits, long& xyBits, long& wBits, long& pBits, long& aBits) {
	//TODO implement method
}

void CipherGD::encStepLGD5(Cipher*& cxyData, Cipher*& cwData, ZZX& poly, long& cnum, double& gamma, long& sBits, long& bBits, long& xyBits, long& wBits, long& pBits, long& aBits) {
	//TODO implement method
}

void CipherGD::encStepLGD7(Cipher*& cxyData, Cipher*& cwData, ZZX& poly, long& cnum, double& gamma, long& sBits, long& bBits, long& xyBits, long& wBits, long& pBits, long& aBits) {
	//TODO implement method
}

void CipherGD::encStepMLGD3(Cipher*& cxyData, Cipher*& cwData, Cipher*& cvData, ZZX& poly, long& cnum, double& gamma, double& eta, double& etaprev, long& sBits, long& bBits, long& xyBits, long& wBits, long& pBits, long& aBits) {
	//TODO implement method
}

void CipherGD::encStepMLGD5(Cipher*& cxyData, Cipher*& cwData, Cipher*& cvData, ZZX& poly, long& cnum, double& gamma, double& eta, double& etaprev, long& sBits, long& bBits, long& xyBits, long& wBits, long& pBits, long& aBits) {
	//TODO implement method
}

void CipherGD::encStepMLGD7(Cipher*& cxyData, Cipher*& cwData, Cipher*& cvData, ZZX& poly, long& cnum, double& gamma, double& eta, double& etaprev, long& sBits, long& bBits, long& xyBits, long& wBits, long& pBits, long& aBits) {
	//TODO implement method
}

void CipherGD::encStepNLGD3(Cipher*& cxyData, Cipher*& cwData, Cipher*& cvData, ZZX& poly, long& cnum, double& gamma, double& eta, double& etaprev, long& sBits, long& bBits, long& xyBits, long& wBits, long& pBits, long& aBits) {
	Cipher* cprod = new Cipher[cnum];
 	Cipher* cgrad = new Cipher[cnum];

	long bitsDown = cxyData[0].cbits - cwData[0].cbits;

	NTL_EXEC_RANGE(cnum, first, last);
	for (long i = first; i < last; ++i) {
		scheme.modEmbedAndEqual(cxyData[i], bitsDown);
		cprod[i] = scheme.mult(cxyData[i], cwData[i]);
		for (long l = 0; l < bBits; ++l) {
			cgrad[i] = scheme.leftRotateByPo2(cprod[i], l);
			scheme.addAndEqual(cprod[i], cgrad[i]);
		}
	}
	NTL_EXEC_RANGE_END;
	Cipher cip = cprod[0];
	for (long i = 1; i < cnum; ++i) {
		scheme.addAndEqual(cip, cprod[i]);
	}
	scheme.modSwitchAndEqual(cip, xyBits);

	scheme.multByPolyAndEqual(cip, poly);
	for (long l = 0; l < bBits; ++l) {
		cprod[0] = scheme.rightRotateByPo2(cip, l);
		scheme.addAndEqual(cip, cprod[0]);
	}
	scheme.modSwitchAndEqual(cip, pBits + aBits); // (xyBits + pBits + aBits, wBits - aBits)

 	Cipher cip2 = scheme.square(cip);
 	scheme.modSwitchAndEqual(cip2, wBits); // cip2 (xyBits + pBits + aBits + wBits, wBits - 2aBits)

 	RR tmp = to_RR(degree3[1]) / degree3[2];
 	ZZ tmpzz = EvaluatorUtils::evaluateVal(tmp, wBits - 2 * aBits);
 	scheme.addConstAndEqual(cip2, tmpzz); // cip2 (xyBits + pBits + aBits + wBits, wBits - 2aBits)

 	tmp = to_RR(gamma) * (1 - eta) * degree3[0];
 	tmpzz = EvaluatorUtils::evaluateVal(tmp, wBits);

 	tmp = to_RR(gamma) * (1 - eta) * degree3[2];
 	ZZ tmpzz1 = EvaluatorUtils::evaluateVal(tmp, wBits + 3 * aBits);

	NTL_EXEC_RANGE(cnum, first, last);
	for (long i = first; i < last; ++i) {
		cprod[i] = scheme.multByConst(cxyData[i], tmpzz); // cprod (0, xyBits + wBits)
		scheme.modSwitchAndEqual(cprod[i], xyBits); // cprod (xyData, wBits)
		scheme.modEmbedAndEqual(cprod[i], pBits  + aBits + 2 * wBits); // cprod (xyData + pBits + aBits + 2wBits, wBits)

		cgrad[i] = scheme.multByConst(cxyData[i], tmpzz1); // cgrad (0, xyData + wBits + 3aBits)
		scheme.modSwitchAndEqual(cgrad[i], xyBits); // cgrad (xyData, wBits + 3aBits)
		scheme.modEmbedAndEqual(cgrad[i], pBits + aBits); // cgrad (xyData + pBits + aBits, wBits + 3aBits)
		scheme.multAndEqual(cgrad[i], cip); // cgrad (xyData + pBits + aBits, 2wBits + 2aBits)
		scheme.modSwitchAndEqual(cgrad[i], wBits); // cgrad (xyData + pBits + aBits + wBits, wBits + 2aBits)
		scheme.multAndEqual(cgrad[i], cip2); // cgrad (xyData + pBits + aBits + 2wBits, 2wBits)
		scheme.modSwitchAndEqual(cgrad[i], wBits); // cgrad (xyData + pBits + aBits + 2wBits, wBits)
		scheme.addAndEqual(cgrad[i], cprod[i]); // cgrad (xyData + pBits + aBits + 2wBits, wBits)
	}
	NTL_EXEC_RANGE_END;

	NTL_EXEC_RANGE(cnum, first, last);
	for (long i = first; i < last; ++i) {
		for (long l = bBits; l < sBits; ++l) {
			cprod[i] = scheme.leftRotateByPo2(cgrad[i], l);
			scheme.addAndEqual(cgrad[i], cprod[i]);
		}
	}
	NTL_EXEC_RANGE_END;

	delete[] cprod;

	tmp = to_RR(1. - eta);
	tmpzz = EvaluatorUtils::evaluateVal(tmp, wBits);

	bitsDown = xyBits + pBits + aBits + wBits;

	if(etaprev < 0) {
		tmp = to_RR(eta / (1 - etaprev));
		tmpzz1 = EvaluatorUtils::evaluateVal(tmp, wBits);

		NTL_EXEC_RANGE(cnum, first, last);
		for (long i = first; i < last; ++i) {
			scheme.multByConstAndEqual(cwData[i], tmpzz);
			scheme.modSwitchAndEqual(cwData[i], wBits);
			scheme.modEmbedAndEqual(cwData[i], bitsDown);
			cgrad[i] = scheme.sub(cwData[i], cgrad[i]);

			scheme.multByConstAndEqual(cvData[i], tmpzz1);
			scheme.modSwitchAndEqual(cvData[i], wBits);
			scheme.modEmbedAndEqual(cvData[i], bitsDown);

			cwData[i] = scheme.add(cvData[i], cgrad[i]);
			cvData[i] = cgrad[i];
		}
		NTL_EXEC_RANGE_END;
	} else {
		NTL_EXEC_RANGE(cnum, first, last);
		for (long i = first; i < last; ++i) {
			cvData[i] = scheme.multByConst(cwData[i], tmpzz);
			scheme.modSwitchAndEqual(cvData[i], wBits);
			scheme.modEmbedAndEqual(cvData[i], bitsDown);
			scheme.subAndEqual(cvData[i], cgrad[i]);
			scheme.modEmbedAndEqual(cwData[i], bitsDown + wBits);
			scheme.subAndEqual(cwData[i], cgrad[i]);
		}
		NTL_EXEC_RANGE_END;
	}
	delete[] cgrad;
}

void CipherGD::encStepNLGD5(Cipher*& cxyData, Cipher*& cwData, Cipher*& cvData, ZZX& poly, long& cnum, double& gamma, double& eta, double& etaprev, long& sBits, long& bBits, long& xyBits, long& wBits, long& pBits, long& aBits) {
	Cipher* cprod = new Cipher[cnum];
 	Cipher* cgrad = new Cipher[cnum];

	long bitsDown = cxyData[0].cbits - cwData[0].cbits;

	NTL_EXEC_RANGE(cnum, first, last);
	for (long i = first; i < last; ++i) {
		scheme.modEmbedAndEqual(cxyData[i], bitsDown); // cxy(0, xyBits)
		cprod[i] = scheme.mult(cxyData[i], cwData[i]);
		for (long l = 0; l < bBits; ++l) {
			cgrad[i] = scheme.leftRotateByPo2(cprod[i], l);
			scheme.addAndEqual(cprod[i], cgrad[i]);
		}
	}
	NTL_EXEC_RANGE_END;

	Cipher cip = cprod[0]; // cip (0, xyBits + wBits)
	for (long i = 1; i < cnum; ++i) {
		scheme.addAndEqual(cip, cprod[i]);
	}
	scheme.modSwitchAndEqual(cip, xyBits); // cip (xyBits, wBits)

	scheme.multByPolyAndEqual(cip, poly);
	for (long l = 0; l < bBits; ++l) {
		cprod[0] = scheme.rightRotateByPo2(cip, l);
		scheme.addAndEqual(cip, cprod[0]);
	}
	scheme.modSwitchAndEqual(cip, pBits + aBits); // cip (xyBits + pBits + aBits, wBits - aBits)

 	Cipher cip2 = scheme.square(cip);
 	scheme.modSwitchAndEqual(cip2, wBits); // cip2 (xyBits + pBits + aBits + wBits, wBits - 2aBits)
 	Cipher cip4 = scheme.square(cip2);

 	RR tmp = to_RR(degree5[2]) / degree5[3];
 	ZZ tmpzz = EvaluatorUtils::evaluateVal(tmp, wBits - 2 * aBits);

 	scheme.multByConstAndEqual(cip2, tmpzz);  // cip2 (xyBits + pBits + aBits + wBits, 2wBits - 4aBits)
 	scheme.modSwitchAndEqual(cip2, wBits); // cip2 (xyBits + pBits + aBits + 2wBits, wBits - 4aBits)

 	scheme.modSwitchAndEqual(cip4, wBits); // cip4 (xyBits + pBits + aBits + 2wBits, wBits - 4aBits)
 	scheme.addAndEqual(cip4, cip2);

 	tmp = to_RR(degree5[1]) / degree5[3];
 	tmpzz = EvaluatorUtils::evaluateVal(tmp, wBits - 4 * aBits);
 	scheme.addConstAndEqual(cip4, tmpzz); // cip4 (xyBits + pBits + aBits + 2wBits, wBits - 4aBits)


 	tmp = to_RR(gamma) * (1 - eta) * degree5[0];
 	tmpzz = EvaluatorUtils::evaluateVal(tmp, wBits);

 	tmp = to_RR(gamma) * (1 - eta) * degree5[3];
 	ZZ tmpzz1 = EvaluatorUtils::evaluateVal(tmp, wBits + 5 * aBits);

	NTL_EXEC_RANGE(cnum, first, last);
	for (long i = first; i < last; ++i) {
		cprod[i] = scheme.multByConst(cxyData[i], tmpzz); // cprod (0, xyBits + wBits)
		scheme.modSwitchAndEqual(cprod[i], xyBits); // cprod (xyData, wBits)
		scheme.modEmbedAndEqual(cprod[i], pBits  + aBits + 3 * wBits); // cprod (xyData + pBits + aBits + 3wBits, wBits)

		cgrad[i] = scheme.multByConst(cxyData[i], tmpzz1); // cgrad (0, xyData + wBits + 5aBits)
		scheme.modSwitchAndEqual(cgrad[i], xyBits); // cgrad (xyData, wBits + 5aBits)
		scheme.modEmbedAndEqual(cgrad[i], pBits + aBits); // cgrad (xyData + pBits + aBits, wBits + 5aBits)
		scheme.multAndEqual(cgrad[i], cip); // cgrad (xyData + pBits + aBits, 2wBits + 4aBits)
		scheme.modSwitchAndEqual(cgrad[i], wBits); // cgrad (xyData + pBits + aBits + wBits, wBits + 4aBits)
		scheme.modEmbedAndEqual(cgrad[i], wBits); // cgrad (xyData + pBits + aBits + 2wBits, wBits + 4aBits)
		scheme.multAndEqual(cgrad[i], cip4); // cgrad (xyData + pBits + aBits + 2wBits, 2wBits)
		scheme.modSwitchAndEqual(cgrad[i], wBits); // cgrad (xyData + pBits + aBits + 3wBits, wBits)
		scheme.addAndEqual(cgrad[i], cprod[i]); // cgrad (xyData + pBits + aBits + 3wBits, wBits)
	}
	NTL_EXEC_RANGE_END;

	NTL_EXEC_RANGE(cnum, first, last);
	for (long i = first; i < last; ++i) {
		for (long l = bBits; l < sBits; ++l) {
			cprod[i] = scheme.leftRotateByPo2(cgrad[i], l);
			scheme.addAndEqual(cgrad[i], cprod[i]);
		}
	}
	NTL_EXEC_RANGE_END;
	delete[] cprod;

	tmp = to_RR(1. - eta);
	tmpzz = EvaluatorUtils::evaluateVal(tmp, wBits);

	bitsDown = xyBits + pBits + aBits + 2 * wBits;

	if(etaprev < 0) {
		tmp = to_RR(eta) / (1 - etaprev);
		tmpzz1 = EvaluatorUtils::evaluateVal(tmp, wBits);
		NTL_EXEC_RANGE(cnum, first, last);
		for (long i = first; i < last; ++i) {
			scheme.multByConstAndEqual(cwData[i], tmpzz);
			scheme.modSwitchAndEqual(cwData[i], wBits);
			scheme.modEmbedAndEqual(cwData[i], bitsDown);
			cgrad[i] = scheme.sub(cwData[i], cgrad[i]);

			scheme.multByConstAndEqual(cvData[i], tmpzz1);
			scheme.modSwitchAndEqual(cvData[i], wBits);
			scheme.modEmbedAndEqual(cvData[i], bitsDown);

			cwData[i] = scheme.add(cvData[i], cgrad[i]);
			cvData[i] = cgrad[i];
		}
		NTL_EXEC_RANGE_END;
	} else {
		NTL_EXEC_RANGE(cnum, first, last);
		for (long i = first; i < last; ++i) {
			cvData[i] = scheme.multByConst(cwData[i], tmpzz);
			scheme.modSwitchAndEqual(cvData[i], wBits);
			scheme.modEmbedAndEqual(cvData[i], bitsDown);
			scheme.subAndEqual(cvData[i], cgrad[i]);
			scheme.modEmbedAndEqual(cwData[i], bitsDown + wBits);
			scheme.subAndEqual(cwData[i], cgrad[i]);
		}
		NTL_EXEC_RANGE_END;
	}
	delete[] cgrad;
}

void CipherGD::encStepNLGD7(Cipher*& cxyData, Cipher*& cwData, Cipher*& cvData, ZZX& poly, long& cnum, double& gamma, double& eta, double& etaprev, long& sBits, long& bBits, long& xyBits, long& wBits, long& pBits, long& aBits) {
	Cipher* cprod = new Cipher[cnum];
 	Cipher* cgrad = new Cipher[cnum];

	long bitsDown = cxyData[0].cbits - cwData[0].cbits;

	NTL_EXEC_RANGE(cnum, first, last);
	for (long i = first; i < last; ++i) {
		scheme.modEmbedAndEqual(cxyData[i], bitsDown); // cxy (0, xyBits)
		cprod[i] = scheme.mult(cxyData[i], cwData[i]);
		for (long l = 0; l < bBits; ++l) {
			cgrad[i] = scheme.leftRotateByPo2(cprod[i], l);
			scheme.addAndEqual(cprod[i], cgrad[i]);
		}
	}
	NTL_EXEC_RANGE_END;

	Cipher cip = cprod[0]; // cip (0, xyBits + wBits)
	for (long i = 1; i < cnum; ++i) {
		scheme.addAndEqual(cip, cprod[i]);
	}
	scheme.modSwitchAndEqual(cip, xyBits); // cip (xyBits, wBits)

	scheme.multByPolyAndEqual(cip, poly);
	for (long l = 0; l < bBits; ++l) {
		cprod[0] = scheme.rightRotateByPo2(cip, l);
		scheme.addAndEqual(cip, cprod[0]);
	}
	scheme.modSwitchAndEqual(cip, pBits + aBits); // cip (xyBits + pBits + aBits, wBits - aBits)

 	Cipher cip2 = scheme.square(cip);
 	scheme.modSwitchAndEqual(cip2, wBits); // cip2 (xyBits + pBits + aBits + wBits, wBits - 2aBits)

 	Cipher cip4 = scheme.square(cip2);
 	scheme.modSwitchAndEqual(cip4, wBits); // cip4 (xyBits + pBits + aBits + 2wBits, wBits - 4aBits)

 	RR tmp = to_RR(degree7[3]) / degree7[4];
 	ZZ tmpzz = EvaluatorUtils::evaluateVal(tmp, wBits - 2 * aBits);
 	Cipher cip2c = scheme.multByConst(cip2, tmpzz);  // cip2c (xyBits + pBits + aBits + wBits, 2wBits - 4aBits)
 	scheme.modSwitchAndEqual(cip2c, wBits); // cip2c (xyBits + pBits + aBits + 2wBits, wBits - 4aBits)
 	scheme.addAndEqual(cip4, cip2c);

 	tmp = to_RR(degree7[2]) / degree7[4];
 	tmpzz = EvaluatorUtils::evaluateVal(tmp, wBits - 4 * aBits);
 	scheme.addConstAndEqual(cip4, tmpzz); // cip4 (xyBits + pBits + aBits + 2wBits, wBits - 4aBits)

 	tmp = to_RR(gamma) * (1 - eta) * degree7[0];
 	tmpzz = EvaluatorUtils::evaluateVal(tmp, wBits);

 	tmp = to_RR(gamma) * (1 - eta) * degree7[1];
 	ZZ tmpzz1 = EvaluatorUtils::evaluateVal(tmp, wBits + aBits);

 	tmp = to_RR(gamma) * (1 - eta) * degree7[4];
 	ZZ tmpzz2 = EvaluatorUtils::evaluateVal(tmp, wBits + 7 * aBits);

	NTL_EXEC_RANGE(cnum, first, last);
	for (long i = first; i < last; ++i) {
		cgrad[i] = scheme.multByConst(cxyData[i], tmpzz); // cgrad (0, xyBits + wBits)
		scheme.modSwitchAndEqual(cgrad[i], xyBits); // cgrad (xybits, wBits)
		scheme.modEmbedAndEqual(cgrad[i], pBits  + aBits + 3 * wBits); // cgrad (xyBits + pBits + aBits + 3wBits, wBits)

		cprod[i] = scheme.multByConst(cxyData[i], tmpzz1); // cprod (0, xyBits + wBits + aBits);
		scheme.modSwitchAndEqual(cprod[i], xyBits); // cprod (xyBits, wBits + aBits);
		scheme.modEmbedAndEqual(cprod[i], pBits  + aBits); // cprod (xyBits + pBits + aBits, wBits + aBits);
		scheme.multAndEqual(cprod[i], cip); // cprod (xyBits + pBits + aBits, 2wBits);
		scheme.modSwitchAndEqual(cprod[i], wBits); // cprod (xyBits + pBits + aBits + wBits, wBits);
		scheme.modEmbedAndEqual(cprod[i], 2 * wBits); // cprod (xyBits + pBits + aBits + 3wBits, wBits);
		scheme.addAndEqual(cprod[i], cgrad[i]); // cprod (xyBits + pBits + aBits + 3wBits, wBits);

		cgrad[i] = scheme.multByConst(cxyData[i], tmpzz2); // cgrad (0, xyData + wBits + 7aBits)
		scheme.modSwitchAndEqual(cgrad[i], xyBits); // cgrad (xyBits, wBits + 7aBits)
		scheme.modEmbedAndEqual(cgrad[i], pBits + aBits); // cgrad (xyBits + pBits + aBits, wBits + 7aBits)
		scheme.multAndEqual(cgrad[i], cip); // cgrad (xyData + pBits + aBits, 2wBits + 6aBits)
		scheme.modSwitchAndEqual(cgrad[i], wBits); // cgrad (xyBits + pBits + aBits + wBits, wBits + 6aBits)
		scheme.multAndEqual(cgrad[i], cip2); // cgrad (xyBits + pBits + aBits + wBits, 2wBits + 4aBits)
		scheme.modSwitchAndEqual(cgrad[i], wBits); // cgrad (xyData + pBits + aBits + 2wBits, wBits + 4aBits)
		scheme.multAndEqual(cgrad[i], cip4); // cgrad (xyData + pBits + aBits + 2wBits, 2wBits)
		scheme.modSwitchAndEqual(cgrad[i], wBits); // cgrad (xyData + pBits + aBits + 3wBits, wBits)
		scheme.addAndEqual(cgrad[i], cprod[i]); // cgrad (xyData + pBits + aBits + 3wBits, wBits)
	}
	NTL_EXEC_RANGE_END;

	NTL_EXEC_RANGE(cnum, first, last);
	for (long i = first; i < last; ++i) {
		for (long l = bBits; l < sBits; ++l) {
			cprod[i] = scheme.leftRotateByPo2(cgrad[i], l);
			scheme.addAndEqual(cgrad[i], cprod[i]);
		}
	}
	NTL_EXEC_RANGE_END;
	delete[] cprod;

	tmp = to_RR(1. - eta);
	tmpzz = EvaluatorUtils::evaluateVal(tmp, wBits);

	bitsDown = xyBits + pBits + aBits + 2 * wBits;

	if(etaprev < 0) {
		tmp = to_RR(eta) / (1 - etaprev);
		tmpzz1 = EvaluatorUtils::evaluateVal(tmp, wBits);
		NTL_EXEC_RANGE(cnum, first, last);
		for (long i = first; i < last; ++i) {
			scheme.multByConstAndEqual(cwData[i], tmpzz);
			scheme.modSwitchAndEqual(cwData[i], wBits);
			scheme.modEmbedAndEqual(cwData[i], bitsDown);
			cgrad[i] = scheme.sub(cwData[i], cgrad[i]);

			scheme.multByConstAndEqual(cvData[i], tmpzz1);
			scheme.modSwitchAndEqual(cvData[i], wBits);
			scheme.modEmbedAndEqual(cvData[i], bitsDown);

			cwData[i] = scheme.add(cvData[i], cgrad[i]);
			cvData[i] = cgrad[i];
		}
		NTL_EXEC_RANGE_END;
	} else {
		NTL_EXEC_RANGE(cnum, first, last);
		for (long i = first; i < last; ++i) {
			cvData[i] = scheme.multByConst(cwData[i], tmpzz);
			scheme.modSwitchAndEqual(cvData[i], wBits);
			scheme.modEmbedAndEqual(cvData[i], bitsDown);
			scheme.subAndEqual(cvData[i], cgrad[i]);
			scheme.modEmbedAndEqual(cwData[i], bitsDown + wBits);
			scheme.subAndEqual(cwData[i], cgrad[i]);
		}
		NTL_EXEC_RANGE_END;
	}
	delete[] cgrad;
}

void CipherGD::decwData(double*& wData, Cipher*& cwData, long& factorDim, long& batch, long& cnum, long& wBits) {
	for (long i = 0; i < (cnum - 1); ++i) {
		CZZ* dcw = scheme.decrypt(secretKey, cwData[i]);
		for (long l = 0; l < batch; ++l) {
			RR wi = to_RR(dcw[l].r);
			wi.e -= wBits;
			wData[batch * i + l] = to_double(wi);
		}
		delete[] dcw;
	}
	CZZ* dcw = scheme.decrypt(secretKey, cwData[cnum-1]);
	long rest = factorDim - batch * (cnum - 1);
	for (long l = 0; l < rest; ++l) {
		RR wi = to_RR(dcw[l].r);
		wi.e -= wBits;
		wData[batch * (cnum - 1) + l] = to_double(wi);
	}
	delete[] dcw;
}
