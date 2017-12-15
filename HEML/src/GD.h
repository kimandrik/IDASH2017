#ifndef SGD_SGD_H_
#define SGD_SGD_H_

#include <iostream>

using namespace std;

class GD {

public:

	static long** xyDataFromFile(string& path, long& factorDim, long& sampleDim, bool isfirst = true);

	static long** RandomxyDataLearn(long** xyData, long& learnDim, long& sampleDim, long& factorDim);

	static double innerprod(double* w, long* xy, long& size);

	static void stepLGD(long** xyData, double* wData, long& factorDim, long& learnDim, double& gamma);

	static void stepNLGD(long** xyData, double* wData, double* vData, long& factorDim, long& learnDim, double& gamma, double& eta);

	static void stepMLGD(long** xyData, double* wData, double* vData, long& factorDim, long& learnDim, double& gamma, double& eta);

	static void check(long** xyData, double* wData, long& factorDim, long& sampleDim);

	static long* calculateYtrueData(long** xyData, long& sampleDim);

	static double* calculateYpredictData(long** xyData, double* wData, long& factorDim, long& sampleDim);

	static double calcuateAUC(long** xyData, double* wData, long& factorDim, long& sampleDim, long steps);
};

#endif /* SGD_SGD_H_ */
