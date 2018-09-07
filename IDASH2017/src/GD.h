/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/

#ifndef IDASH2017_GD_H_
#define IDASH2017_GD_H_

#include <iostream>

static double degree3[3] = {-0.5,0.15012,-0.001593};
static double degree5[4] = {-0.5,0.19131,-0.0045963, 0.0000412332};
static double degree7[5] = {-0.5,0.216884,-0.00819276,0.000165861,-0.00000119581};

//static double degree3[3] = {-0.5,0.19,-0.0035};
//static double degree5[4] = {-0.5,0.2166,-0.0077,0.00011};
//static double degree7[5] = {-0.5,0.216884,-0.00819276,0.000165861,-0.00000119581};

using namespace std;

class GD {

public:

	static double** zDataFromFile(string& path, long& factorDim, long& sampleDim, bool isfirst = true);

	static void shuffleZData(double** zData, long factorDim, long sampleDim);

	static void normalizeZData(double** zData, long factorDim, long sampleDim);
	static void normalizezData2(double** zDataLearn, double** zDataTest, long factorDim, long sampleDimLearn, long sampleDimTest);

	static void initialWDataVDataAverage(double* wData, double* vData, double** zData, long factorDim, long sampleDim);
	static void initialWDataVDataZero(double* wData, double* vData, long factorDim);

	static double* plainIP(double** a, double* b, long factorDim, long sampleDim);
	static double* plainSigmoid(long approxDeg, double** zData, double* ip, long factorDim, long sampleDim, double gamma);

	static void plainLGDstep(double* wData, double* grad, long factorDim);
	static void plainMLGDstep(double* wData, double* vData, double* grad, long factorDim, double eta);
	static void plainNLGDstep(double* wData, double* vData, double* grad, long factorDim, double eta);

	static void plainLGDL2step(double* wData, double* grad, long factorDim, double lambda);
	static void plainMLGDL2step(double* wData, double* vData, double* grad, long factorDim, double eta, double lambda);
	static void plainNLGDL2step(double* wData, double* vData, double* grad, long factorDim, double eta, double lambda);

	static void plainLGDiteration(long approxDeg, double** zData, double* wData, long factorDim, long sampleDim, double gamma);
	static void plainMLGDiteration(long approxDeg, double** zData, double* wData, double* vData, long factorDim, long sampleDim, double gamma, double eta);
	static void plainNLGDiteration(long approxDeg, double** zData, double* wData, double* vData, long factorDim, long sampleDim, double gamma, double eta);

	static void plainLGDL2iteration(long approxDeg, double** zData, double* wData, long factorDim, long sampleDim, double gamma, double lambda);
	static void plainMLGDL2iteration(long approxDeg, double** zData, double* wData, double* vData, long factorDim, long sampleDim, double gamma, double eta, double lambda);
	static void plainNLGDL2iteration(long approxDeg, double** zData, double* wData, double* vData, long factorDim, long sampleDim, double gamma, double eta, double lambda);

	//-----------------------------------------

	static double trueIP(double* a, double* b, long size);

	static void trueLGDiteration(double** zData, double* wData, long factorDim, long sampleDim, double gamma);
	static void trueMLGDiteration(double** zData, double* wData, double* vData, long factorDim, long sampleDim, double gamma, double eta);
	static void trueNLGDiteration(double** zData, double* wData, double* vData, long factorDim, long sampleDim, double gamma, double eta);

	static void trueLGDL2iteration(double** zData, double* wData, long factorDim, long sampleDim, double gamma, double lambda);
	static void trueMLGDL2iteration(double** zData, double* wData, double* vData, long factorDim, long sampleDim, double gamma, double eta, double lambda);
	static void trueNLGDL2iteration(double** zData, double* wData, double* vData, long factorDim, long sampleDim, double gamma, double eta, double lambda);

	static void calculateAUC(double** zData, double* wData, long factorDim, long sampleDim, double& correctness, double& AUC);
	static double calculateMSE(double* wData1, double* wData2, long factorDim);
	static double calculateNMSE(double* wData1, double* wData2, long factorDim);
};

#endif /* SGD_SGD_H_ */
