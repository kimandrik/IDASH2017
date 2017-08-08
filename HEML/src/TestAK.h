#ifndef TEST_TESTSGD_H_
#define TEST_TESTSGD_H_

#include <string>

using namespace std;

class TestAK {
public:
	static void testNLGDWB();
	static void testNLGDXYB(string filename, long iter, long logl, long logp, double gammaCnst, bool is3approx, bool isAllsample, bool isEncrypted);
};

#endif /* TEST_TESTSGD_H_ */
