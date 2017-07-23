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

	long** zdataFromFile(string& path, long& dim, long& sampledim); // returns x_ij * y_i

	//-----------------------------------------

	double** wdatagen(long& wnum, long& dim);
	double* gammagen(long& iter);
	double innerprod(double*& wdata, long*& x, long& size);

	void steplogregress(double*& wdata, long**& zdata, double& gamma, double& lambda, long& dim, long& learndim);

	void stepsimpleregress(double*& wdata, long**& zdata, double& gamma, double& lambda, long& dim, long& learndim);

	double* wout(double**& wdata, long& wnum, long& dim);
	void check(double*& w, long**& zdata, long& dim, long& sampledim);

	//-----------------------------------------

};

#endif /* SGD_SGD_H_ */
