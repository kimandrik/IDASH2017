#ifndef HEML_TESTGD_H_
#define HEML_TESTGD_H_

class TestGD {
public:

	static long suggestLogN(long lambda, long logQ);

	static void testEncNLGD(double** zDataTrain, double** zDataTest, long factorDim, long sampleDimTrain, long sampleDimTest,
			bool isYfirst, long numIter, long k, double gammaUp, double gammaDown, bool isInitZero);

	static void testEncNLGDFOLD(long fold, double** zData, long factorDim, long sampleDim,
			bool isYfirst, long numIter, long k, double gammaUp, double gammaDown, bool isInitZero);

	static void testPlainNLGD(double** zDataTrain, double** zDataTest, long factorDim, long sampleDimTrain, long sampleDimTest,
			bool isYfirst, long numIter, long k, double gammaUp, double gammaDown, bool isInitZero);

	static void testPlainNLGDFOLD(long fold, double** zData, long factorDim, long sampleDim,
			bool isYfirst, long numIter, long k, double gammaUp, double gammaDown, bool isInitZero);

};

#endif /* TESTGD_H_ */
