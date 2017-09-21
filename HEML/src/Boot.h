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

	void bootstrapAndEqual(Cipher& cipher, long logq0, long logT, long logI);
};

#endif /* BOOT_H_ */
