#ifndef BOOT_H_
#define BOOT_H_

#include "BootKey.h"
#include <NTL/ZZ.h>
#include <NTL/RR.h>
#include <Cipher.h>
#include <Scheme.h>
#include <SecKey.h>
#include <map>

using namespace NTL;

class Boot {
public:

	Scheme& scheme;
	SecKey& secretKey;
	map<long, BootKey> bootKeyMap;

	double sinDegree3[2] = {0.99749,-0.15652};
	double sinDegree5[3] = {0.99997860,-0.16649840,0.00799232};
	double sinDegree7[4] = {0.9999999401,-0.1666654105,0.008328769024,-0.000192297728};

	double cosDegree2[2] = {1,-0.0035};
	double cosDegree4[3] = {1,0.2166,-0.0077};
	double cosDegree6[4] = {1,0.216884,-0.00819276,0.000165861};

	ZZX* axBabyRot;
	ZZX* bxBabyRot;

	Boot(Scheme& scheme, SecKey& secretKey);

	void generateKey(long l, long pBits);

	virtual ~Boot();

	void normalize(Cipher& cipher, long N);

	void linearTransform(Cipher& enclin, Cipher& encx, long size);
	void linearTransformInv(Cipher& res, Cipher& enclin, long size);

	Cipher leftBabyRotate(Cipher& cipher, long i);
	void leftGiantRotateAndEqual(Cipher& cipher, long l, long k, long i);

	Cipher evaluateEncSin2pix7(Cipher& enclinx, long precisionBits);
	Cipher evaluateEncSin2pix3(Cipher& enclinx, long precisionBits);

	Cipher evaluateEncCos2pix6(Cipher& enclinx, long precisionBits);
	Cipher evaluateEncCos2pix2(Cipher& enclinx, long precisionBits);

	Cipher evaluateEncSin2x(Cipher& encSinx, Cipher& encCosx, long precisionBits);
	Cipher evaluateEncCos2x(Cipher& encSinx, Cipher& encCosx, long precisionBits);

	Cipher removeIpart(Cipher& cipher, long logq0, long logT, long logI);

	Cipher bootstrap(Cipher& encx, long logq0, long logT, long logI);
};

#endif /* BOOT_H_ */
