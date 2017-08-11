#ifndef TEST_TESTSGD_H_
#define TEST_TESTSGD_H_

#include <string>

using namespace std;

class TestAK {
public:
	static void testNLGD(string filename, long iter, double gammaCnst, double gammaUpCnst, bool is3approx, bool isAllsample, bool isEncrypted, bool isYfirst, long xyBits, long wBits, long pBits, long lBits);
};

#endif /* TEST_TESTSGD_H_ */
