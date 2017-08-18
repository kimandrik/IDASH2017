#ifndef TEST_TESTSGD_H_
#define TEST_TESTSGD_H_

#include <string>

using namespace std;

class TestAK {
public:
	static void testNLGD(string filename, long iter, double gammaDownCnst, double gammaUpCnst, double learnPortion, bool is3approx, bool isEncrypted, bool isYfirst);
};

#endif /* TEST_TESTSGD_H_ */
