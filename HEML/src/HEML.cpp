#include <Params.h>
#include <PubKey.h>
#include <Scheme.h>
#include <SchemeAlgo.h>
#include <SchemeAux.h>
#include <SecKey.h>
#include <TimeUtils.h>
#include <TestScheme.h>

#include "TestAK.h"
#include "TestJH.h"
#include "TestKW.h"

#include <string>

using namespace std;

void test() {

}

int main() {

//	string filename = "data/data5x500.txt";     // false   415/500
//	string filename = "data/data9x1253.txt";    // false   775/1253 not good results
//	string filename = "data/data15x1500.txt";   // false   1270/1500
//	string filename = "data/data16x101.txt";    // false   101/101
//	string filename = "data/data27x148.txt";    // false   132/148
//	string filename = "data/data43x3247.txt";   // false   3182/3247
//	string filename = "data/data45x296.txt";    // false   257/296
//	string filename = "data/data51x653.txt";    // false   587/653
//	string filename = "data/data67x216.txt";    // false   216/216 slow convergence
	string filename = "data/data103x1579.txt";  // true    1086/1579 has many troubles with ip

//	14: logq <= 310 (secure parameters)
//	15: logq <= 620 (secure parameters)
//	16: logq <= 1241 (secure parameters)
//	17: logq <= 2483 (secure parameters)

	bool isYfirst = true;

	long iter = 10;

	long wBits = 30;
	long xyBits = 32;
	long pBits = 20;
	long lBits = 40;

	double gammaCnst = -4.; // if gammaCnst > 0 then gammaUpCnst / learndim / gammaCnst, else gammaUpCnst / learndim / (|gammaCnst| + i)
	double gammaUpCnst = 4.; // gammaUpCnst / ...
	bool is3approx = false; // if true then 3 degree approximation, else 7 degree approximation
	bool isAllsample = true; // if true then all sample is learned, else log2(sample) is learned
	bool isEncrypted = false; // if true then encrypted learn, else unecnrypted (for testing)


	TestAK::testNLGD(filename, iter, gammaCnst, gammaUpCnst, is3approx, isAllsample, isEncrypted, isYfirst, xyBits, wBits, pBits, lBits);

	return 0;
}
