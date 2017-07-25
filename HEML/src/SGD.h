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

	void stepQuadraticRegress(double*& wdata, long**& zdata, double& gamma, double& lambda, long& dim, long& learndim);
	void stepLogRegress(double*& wdata, long**& zdata, double& gamma, double& lambda, long& dim, long& learndim);
	void stepStochasticQuadraticRegress(double*& wdata, long**& zdata, double& gamma, double& lambda, long& dim, long& learndim, long& stochdim);
	void stepStochasticLogRegress(double*& wdata, long**& zdata, double& gamma, double& lambda, long& dim, long& learndim, long& stochdim);
	void stepMomentumLogRegress(double*& wdata, double*& vdata, long**& zdata, double& gamma, long& dim, long& learndim, double& eta);
	void stepNesterovLogRegress(double*& wdata, double*& vdata, long**& zdata, double& gamma, long& dim, long& learndim, double& eta);

	double* waverage(double**& wdata, long& wnum, long& dim);

	void check(double*& w, long**& zdata, long& dim, long& sampledim);

	void debugcheck(string prefix, double*& w, long& dim);
	void debugcheck(string prefix, long*& z, long& dim);
	//-----------------------------------------

};

#endif /* SGD_SGD_H_ */
