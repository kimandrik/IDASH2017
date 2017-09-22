#ifndef BOOTAUX_H_
#define BOOTAUX_H_

#include <NTL/RR.h>

using namespace NTL;

class BootAux {
public:

	RR* sinChebyshevDeg3;
	RR* sinChebyshevDeg5;
	RR* sinChebyshevDeg7;

	RR* cosChebyshevDeg2;
	RR* cosChebyshevDeg4;
	RR* cosChebyshevDeg6;

	RR* sinTaylorDeg7;

	RR* cosTaylorDeg6;

	RR oneOver2pi;

	BootAux();
	virtual ~BootAux();
};

#endif /* BOOTAUX_H_ */
