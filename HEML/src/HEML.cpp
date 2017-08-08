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
		string filename = "data/data15x1500.txt";   // false   1270/1500
	//	string filename = "data/data16x101.txt";    // false   101/101
	//	string filename = "data/data27x148.txt";    // false   132/148
	//	string filename = "data/data43x3247.txt";   // false   3182/3247
	//	string filename = "data/data45x296.txt";    // false   257/296
	//	string filename = "data/data51x653.txt";    // false   587/653
	//	string filename = "data/data67x216.txt";    // false   216/216 slow convergence
	//	string filename = "data/data103x1579.txt";  // true    1086/1579 has many troubles with ip

	long iter = 10;
	long logl = 10;
	long logp = 32;
	double gammaCnst = 2.;
	bool is3approx = false;
	bool isAllsample = true;
	bool isEncrypted = true;

	TestAK::testNLGDXYB(filename, iter, logl, logp, gammaCnst, is3approx, isAllsample, isEncrypted);

//	TestScheme::testExponentBatch(13, 5, 30, 5, 7, 12);
	return 0;
}
