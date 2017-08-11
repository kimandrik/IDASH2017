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

//	string filename = "data/data5x500.txt";     // isYfirst = false   421/500
//	string filename = "data/data9x1253.txt";    // isYfirst = false   1147/1253
//	string filename = "data/data15x1500.txt";   // isYfirst = false   1277/1500
//	string filename = "data/data16x101.txt";    // isYfirst = false   101/101
//	string filename = "data/data27x148.txt";    // isYfirst = false   132/148
//	string filename = "data/data43x3247.txt";   // isYfirst = false   3182/3247
//	string filename = "data/data45x296.txt";    // isYfirst = false   257/296
//	string filename = "data/data51x653.txt";    // isYfirst = false   590/653
//	string filename = "data/data67x216.txt";    // isYfirst = false   216/216
	string filename = "data/data103x1579.txt";  // isYfirst = true    1086/1579

//	13: logq <= 155 (secure parameters)
//	14: logq <= 310 (secure parameters)
//	15: logq <= 620 (secure parameters)
//	16: logq <= 1241 (secure parameters)
//	17: logq <= 2483 (secure parameters)

	bool isYfirst = true;

	long iter = 10;

	long wBits = 32;
	long xyBits = 35;
	long pBits = 20;
	long lBits = 40;

	double gammaCnst = -3.; // if gammaCnst > 0 then gammaUpCnst / learndim / gammaCnst, else gammaUpCnst / learndim / (|gammaCnst| + i)
	double gammaUpCnst = 2.; // gammaUpCnst / ...
	double learnPortion = 0.9; // portion of learnPortion to be learned
	bool is3approx = false; // if true then 3 degree approximation, else 7 degree approximation
	bool isEncrypted = false; // if true then encrypted learn, else unecnrypted (for testing)


	TestAK::testNLGD(filename, iter, gammaCnst, gammaUpCnst, learnPortion, is3approx, isEncrypted, isYfirst, xyBits, wBits, pBits, lBits);
	return 0;
}
