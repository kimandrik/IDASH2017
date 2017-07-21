#ifndef SGD_SGD_H_
#define SGD_SGD_H_

#include <Cipher.h>
#include <Scheme.h>
#include <SchemeAlgo.h>

using namespace std;
using namespace NTL;

class SGD {
public:
	Scheme scheme;
	SchemeAlgo algo;

	//-----------------------------------------

	SGD(Scheme& scheme, SchemeAlgo& algo) : scheme(scheme), algo(algo) {}

	//-----------------------------------------

	long** zdataFromFile(string& path, long& dim, long& sampledim); // returns x_ij * y_i

	//-----------------------------------------

	double** wdatagen(long& wnum, long& dim);
	double* gammagen(long& iter);
	double innerprod(double*& wdata, long*& x, long& size);

	void step(double**& wdata, long**& zdata, double& gamma, double& lambda, long& wnum, long& dim, long& sampledim);
	double* wout(double**& wdata, long& wnum, long& dim);
	void check(double*& w, long**& zdata, long& dim, long& sampledim);

	//-----------------------------------------

	Cipher* enczdata(long**& zdata, long& slots, long& wnum, long& dim, long& sampledim, ZZ& p);
	Cipher* encwdata(double**& wdata, long& slots, long& wnum, long& dim, long& sampledim, long& logp);
	ZZ* pgammagen(double*& alpha, long& iter, long& logp);
	void encStep(Cipher*& czdata, Cipher*& cwdata, ZZ& pgamma, double& lambda, long& slots, long& wnum, long& dim, long& sampledim);
	Cipher* encwout(Cipher*& cwdata, long& wnum, long& dim);
	double* decw(SecKey& secretKey, Cipher*& cw, long& dim);
};

#endif /* SGD_SGD_H_ */
