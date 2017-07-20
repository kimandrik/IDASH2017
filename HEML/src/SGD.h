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

	double plainphi(double& b);
	double plainphiprime(double& b);

	double plainip(double*& wdata, long*& x, long& size);
	double* plainGradient(double*& wdata, long**& zdata, long& dim, long& sampledim, double& lambda);
	double* sgd(long& iter, double**& wdata, long**& zdata, double*& alpha, double& lambda, long& wnum, long& dim, long& sampledim);

	void check(double*& w, long**& zdata, long& dim, long& sampledim);
	//-----------------------------------------

	Cipher* encryptzdata(long**& zdata, long& slots, long& wnum, long& dim, long& sampledim, ZZ& p);
	Cipher* encryptwdata(double**& wdata, long& slots, long& wnum, long& dim, long& sampledim, long& logp);

	Cipher* cipherGradient(Cipher*& zciphers, Cipher*& wciphers, long& slots, long& wnum, long& dim);

	Cipher* ciphersgd(long& iter, Cipher*& zciphers, Cipher*& wciphers, ZZ*& palpha, long& slots, long& wnum, long& dim);

};

#endif /* SGD_SGD_H_ */
