#include <Params.h>
#include <PubKey.h>
#include <Scheme.h>
#include <SchemeAlgo.h>
#include <SchemeAux.h>
#include <SecKey.h>
#include <TimeUtils.h>
#include <TestScheme.h>

#include "TestSGD.h"

void test() {

}

int main() {

//	test();

	//logN, logl, logp, L

	//15, 30, 20, 30
	//16, 60, 20, 60
	//17, 20, 20, 125
	//18, 40, 20, 250
	//19, 80, 20, 500
	//20, 160, 20, 1000 (+)

	//15, 5, 25, 25
	//16, 10, 25, 50
	//17, 20, 25, 100
	//18, 40, 25, 200
	//19, 80, 25, 400
	//20, 160, 25, 800 (+)

	//15, 30, 30, 20
	//16, 60, 30, 40
	//17, 120, 30, 80
	//18, 90, 30, 165
	//19, 30, 30, 335
	//20, 60, 30, 670

//	TestSGD::testSimpleSGD(12, 10, 30, 10);
	TestSGD::testLogSGD(12, 10, 30, 10);
//	TestScheme::testEncodeBatch(13, 5, 30, 5, 12);
//	TestScheme::testExponentBatch(13, 5, 30, 5, 7, 12);
//	TestScheme::testSlotsSum(13, 5, 50, 5, 12);
	//-----------------------------------------

	return 0;
}
