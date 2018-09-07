/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/

#ifndef IDASH2017_TESTGD_H_
#define IDASH2017_TESTGD_H_

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
