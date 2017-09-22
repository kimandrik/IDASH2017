#include "BootAux.h"

#include <EvaluatorUtils.h>

BootAux::BootAux() {
	sinChebyshevDeg3 = new RR[2];
	sinChebyshevDeg5 = new RR[3];
	sinChebyshevDeg7 = new RR[4];

	cosChebyshevDeg2 = new RR[2];
	cosChebyshevDeg4 = new RR[3];
	cosChebyshevDeg6 = new RR[4];

	sinTaylorDeg7 = new RR[4];

	cosTaylorDeg6 = new RR[4];

	sinTaylorDeg7[0] = -4*Pi*Pi*Pi/3; // c3
	sinTaylorDeg7[1] = -3/(2*Pi*Pi); // c1/c3
	sinTaylorDeg7[2] = -8*Pi*Pi*Pi*Pi*Pi/315; // c7
	sinTaylorDeg7[3] = -21/(2*Pi*Pi); // c5/c7

	cosTaylorDeg6[0] = -1/(2*Pi*Pi); // 1/c2
	cosTaylorDeg6[1] = -2*Pi*Pi; // c2
	cosTaylorDeg6[2] = -15/(2*Pi*Pi); // c4/c6
	cosTaylorDeg6[3] = -4*Pi*Pi*Pi*Pi*Pi*Pi/45; // c6

	oneOver2pi = 1/(2*Pi);
}

BootAux::~BootAux() {
}

