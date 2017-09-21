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

	void normalizeAndEqual(Cipher& cipher, long N);

	void linearTransform(Cipher& res, Cipher& cipher, long size);
	void linearTransformInv(Cipher& res, Cipher& cipher, long size);

	void linearTransformAndEqual(Cipher& cipher, long size);
	void linearTransformInvAndEqual(Cipher& cipher, long size);

	Cipher leftBabyRotate(Cipher& cipher, long i);
	void leftGiantRotateAndEqual(Cipher& cipher, long l, long k, long i);

	Cipher evaluateEncSin2pix7(Cipher& cipher, long precisionBits);
	Cipher evaluateEncSin2pix3(Cipher& cipher, long precisionBits);

	Cipher evaluateEncCos2pix6(Cipher& cipher, long precisionBits);
	Cipher evaluateEncCos2pix2(Cipher& cipher, long precisionBits);

	Cipher evaluateEncSin2x(Cipher& cipherSin, Cipher& cipherCos, long precisionBits);
	Cipher evaluateEncCos2x(Cipher& cipherSin, Cipher& cipherCos, long precisionBits);

	Cipher removeIpart(Cipher& cipher, long logq0, long logT, long logI);
	void removeIpartAndEqual(Cipher& cipher, long logq0, long logT, long logI);

	Cipher bootstrap(Cipher& encx, long logq0, long logT, long logI);
	void bootstrapAndEqual(Cipher& cipher, long logq0, long logT, long logI);
};

#endif /* BOOT_H_ */
