#include "BootAux.h"

BootAux::BootAux() {
	sinChebyshevDeg3 = new RR[2];
	sinChebyshevDeg5 = new RR[3];
	sinChebyshevDeg7 = new RR[4];

	cosChebyshevDeg2 = new RR[2];
	cosChebyshevDeg4 = new RR[3];
	cosChebyshevDeg6 = new RR[4];

	sinTaylorDeg3 = new RR[2];
	sinTaylorDeg5 = new RR[3];
	sinTaylorDeg7 = new RR[4];

	cosTaylorDeg2 = new RR[2];
	cosTaylorDeg4 = new RR[3];
	cosTaylorDeg6 = new RR[4];

	//TODO
}

BootAux::~BootAux() {
}

