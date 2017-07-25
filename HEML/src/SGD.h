#ifndef SGD_SGD_H_
#define SGD_SGD_H_

#include <Cipher.h>
#include <Scheme.h>
#include <SchemeAlgo.h>

using namespace std;
using namespace NTL;

class SGD {
public:

	//-----------------------------------------

	SGD() {}

	//-----------------------------------------

	long** xyDataFromFile(string& path, long& factordim, long& sampledim); // returns x_ij * y_i

	//-----------------------------------------

	double innerprod(double*& wdata, long*& x, long& size);

	void stepQGD(double*& wData, long**& xyData, double& gamma, double& lambda, long& factorDim, long& learnDim);
	void stepSQGD(double*& wData, long**& xyData, double& gamma, double& lambda, long& factorDim, long& learnDim, long& stochDim);

	void stepLGD(double*& wData, long**& xyData, double& gamma, double& lambda, long& factorDim, long& learnDim);
	void stepSLGD(double*& wData, long**& xyData, double& gamma, double& lambda, long& factorDim, long& learnDim, long& stochDim);

	void stepMLGD(double*& wData, double*& vData, long**& xyData, double& gamma, long& factorDim, long& learnDim, double& eta);
	void stepNLGD(double*& wData, double*& vData, long**& xyData, double& gamma, long& factorDim, long& learnDim, double& eta);

	double* waverage(double**& wdata, long& wnum, long& dim);

	void check(double*& w, long**& zdata, long& dim, long& sampledim);

	void debugcheck(string prefix, double*& w, long& dim);
	void debugcheck(string prefix, long*& z, long& dim);
	//-----------------------------------------

};

#endif /* SGD_SGD_H_ */
