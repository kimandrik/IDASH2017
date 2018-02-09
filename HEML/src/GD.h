#ifndef HEML_GD_H_
#define HEML_GD_H_

#include <iostream>

//static double degree3[3] = {-0.5,0.15012,-0.001593};
//static double degree5[4] = {-0.5,0.19131,-0.0045963, 0.0000412332};
//static double degree7[5] = {-0.5,0.216884,-0.00819276,0.000165861,-0.00000119581};

static double degree3[3] = {-0.5,0.19,-0.0035};
static double degree5[4] = {-0.5,0.2166,-0.0077,0.00011};
static double degree7[5] = {-0.5,0.216884,-0.00819276,0.000165861,-0.00000119581};

using namespace std;

class GD {

public:

	static double** xyDataFromFile(string& path, long& factorDim, long& sampleDim, bool isfirst = true);
	static double** RandomxyDataLearn(double** xyData, long sampleDimLearn, long sampleDim, long factorDim);
	static void normalizexyData(double** xyData, long factorDim, long sampleDim);
	static void normalizexyData2(double** xyDataLearn, double** xyDataTest, long factorDim, long sampleDimLearn, long sampleDimTest);

	static double* plainIP(double** a, double* b, long factorDim, long learnDim);

	static double* plainSigmoid(long approxDeg, double** xyData, double* ip, long factorDim, long sampleDim, double gamma);

	static void plainLGDstep(double* wData, double* grad, long factorDim);
	static void plainMLGDstep(double* wData, double* vData, double* grad, long factorDim, double eta);
	static void plainNLGDstep(double* wData, double* vData, double* grad, long factorDim, double eta);

	static void plainLGDiteration(long approxDeg, double** xyData, double* wData, long factorDim, long sampleDim, double gamma);
	static void plainMLGDiteration(long approxDeg, double** xyData, double* wData, double* vData, long factorDim, long sampleDim, double gamma, double eta);
	static void plainNLGDiteration(long approxDeg, double** xyData, double* wData, double* vData, long factorDim, long sampleDim, double gamma, double eta);

	//-----------------------------------------

	static double trueIP(double* a, double* b, long size);

	static void trueLGDiteration(double** xyData, double* wData, long factorDim, long sampleDim, double gamma);
	static void trueMLGDiteration(double** xyData, double* wData, double* vData, long factorDim, long sampleDim, double gamma, double eta);
	static void trueNLGDiteration(double** xyData, double* wData, double* vData, long factorDim, long sampleDim, double gamma, double eta);

	static void check(double** xyData, double* wData, long factorDim, long sampleDim);

	static double calculateAUC(double** xyData, double* wData, long factorDim, long sampleDim);

	static double calculateMSE(double* wData1, double* wData2, long factorDim);

	static double calculateNMSE(double* wData1, double* wData2, long factorDim);
};

#endif /* SGD_SGD_H_ */
