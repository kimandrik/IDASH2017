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

void test() {

}

int main() {

//	test();

	//logN, logl, logp, L (some secure parameter sets)

	//15, 5, 25, 25
	//15, 30, 30, 20
	//15, 0, 35, 18
	//15, 30, 40, 15

	//16, 10, 25, 50
	//16, 60, 30, 40
	//16, 0, 35, 36
	//16, 60, 40, 30

	//17, 20, 25, 100
	//17, 120, 30, 80
	//17, 0, 35, 72
	//17, 120, 40, 60

	//18, 40, 25, 200
	//18, 90, 30, 165
	//18, 0, 35, 144
	//18, 40, 40, 125

//	TestAK::testAK(12, 60, 40, 20);
//	TestCH::testJH(16, 60, 30, 40);
//	TestKW::testKW(16, 60, 30, 40);

//	TestScheme::testEncodeBatch(13, 5, 30, 5, 12);
//	TestScheme::testExponentBatch(13, 5, 30, 5, 7, 12);
//	TestScheme::testSlotsSum(13, 5, 50, 5, 12);
//	TestScheme::testInverseBatch(13, 5, 30, 5, 4, 3);

	//-----------------------------------------

	return 0;
}
