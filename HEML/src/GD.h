#ifndef SGD_SGD_H_
#define SGD_SGD_H_

#include <iostream>

using namespace std;

class GD {

public:

	static double** xyDataFromFile(string& path, long& factorDim, long& sampleDim, bool isfirst = true);

	static double** RandomxyDataLearn(double** xyData, long learnDim, long sampleDim, long factorDim);

	static double innerprod(double* w, double* xy, long size);

	static void stepLGD(double** xyData, double* wData, long factorDim, long learnDim, double& gamma);

	static void stepNLGD(double** xyData, double* wData, double* vData, long factorDim, long learnDim, double& gamma, double& eta);

	static void stepNLGDimitate(double** xyData, double* wData, double* vData, long factorDim, long learnDim, double& gamma, double& eta, double& etaprev);

	static void stepMLGD(double** xyData, double* wData, double* vData, long factorDim, long learnDim, double& gamma, double& eta);

	static void check(double** xyData, double* wData, long factorDim, long sampleDim);

	static double* calculateYtrueData(double** xyData, long sampleDim);

	static double* calculateYpredictData(double** xyData, double* wData, long factorDim, long sampleDim);

	static double calcuateAUC(double** xyData, double* wData, long factorDim, long sampleDim, long steps);
};

#endif /* SGD_SGD_H_ */
