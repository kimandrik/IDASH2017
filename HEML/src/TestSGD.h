#ifndef TEST_TESTSGD_H_
#define TEST_TESTSGD_H_

#include <vector>

using namespace std;

class TestSGD {
public:
	static void testSimpleSGD(long logN, long logl, long logp, long L);
	static void testLogSGD(long logN, long logl, long logp, long L);
};

#endif /* TEST_TESTSGD_H_ */
