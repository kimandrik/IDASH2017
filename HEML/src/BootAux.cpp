#include "BootAux.h"

//double sinChebyshevDegree3[2] = {0.99749,-0.15652};
//double sinChebyshevDegree5[3] = {0.99997860,-0.16649840,0.00799232};
//double sinChebyshevDegree7[4] = {0.9999999401,-0.1666654105,0.008328769024,-0.000192297728};
//double cosChebyshevDegree2[2] = {1,-0.0035};
//double cosChebyshevDegree4[3] = {1,0.2166,-0.0077};
//double cosChebyshevDegree6[4] = {1,0.216884,-0.00819276,0.000165861};

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

