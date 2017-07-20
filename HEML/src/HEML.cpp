#include <Params.h>
#include <PubKey.h>
#include <Scheme.h>
#include <SchemeAlgo.h>
#include <SchemeAux.h>
#include <SecKey.h>
#include <TimeUtils.h>

#include "TestBoot.h"
#include "TestHD.h"
#include "TestSGD.h"

void test() {

}

int main() {

//	test();

//	TestBoot::testBoot(6, 5, 30, 5, 5);
	TestSGD::testSGD(13, 5, 30, 5);
//	TestHD::testHD(13, 13, 30, 2, 11);
//	TestHD::testDNA(13, 13, 30, 2);
	//-----------------------------------------

	return 0;
}
