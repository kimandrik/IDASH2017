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

	long** dataFromFile(string& path, long& dim, long& sampledim);

	//-----------------------------------------

	double plaincost(double*& w, long**& data, long& dim, long& sampledim);
	double plainnorm(double*& w, long& size);
	double plainphi(double& b);
	double plainphiprime(double& b);
	double plainip(double*& wdata, long*& x, long& size);
	double* plainGradient(double*& wdata, long**& data, long& dim, long& sampledim, double& lambda);
	void plainsgd(long& iter, double*& wdata, long**& zdata, long& dim, long& sampledim);

	//-----------------------------------------

	Cipher* cipherGradient(Cipher*& zcipher, Cipher*& wcipher, const long& dim, const long& slots, const long& wnum);
};

#endif /* SGD_SGD_H_ */
