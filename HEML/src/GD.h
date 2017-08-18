#ifndef SGD_SGD_H_
#define SGD_SGD_H_

#include <Cipher.h>
#include <Scheme.h>
#include <SchemeAlgo.h>

using namespace std;
using namespace NTL;

class GD {
public:

	//-----------------------------------------

	GD() {}

	//-----------------------------------------

	long** xyDataFromFile(string& path, long& factorDim, long& sampleDim, bool isfirst = true);

	long** RandomxyDataLearn(long**& xyData, long& learnDim, long& sampleDim, long& factorDim);

	double innerprod(double*& w, long*& xy, long& size);

	void stepLGD(long**& xyData, double*& wData, long& factorDim, long& learnDim, double& lambda, double& gamma);

	void stepNLGD(long**& xyData, double*& wData, double*& vData, long& factorDim, long& learnDim, double& gamma, double& eta);

	double* wsum(double**& wdata, long& factorDim, long& wBatch);

	void check(long**& xyData, double*& w, long& factorDim, long& sampleDim);

};

#endif /* SGD_SGD_H_ */
