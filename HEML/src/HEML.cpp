#include <Params.h>
#include <PubKey.h>
#include <Scheme.h>
#include <SchemeAlgo.h>
#include <SchemeAux.h>
#include <SecKey.h>
#include <TimeUtils.h>

#include "TestSGD.h"

void test() {

}

int main() {

//	test();

	TestSGD::testSGD(16, 5, 20, 60);
	//-----------------------------------------

	return 0;
}
