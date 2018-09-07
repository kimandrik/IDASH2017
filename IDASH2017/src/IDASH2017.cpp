/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/

#include <NTL/BasicThreadPool.h>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>

#include "GD.h"
#include "TestGD.h"

using namespace std;
using namespace NTL;

/*
 * run: ./IDASH2017 trainfile isYfirst numIter k gammaUp gammaDown isInitZero fold isEncrypted testfile
 * ./HEML string bool long long double double bool long bool string
 * example: ./IDASH2017 "../data/data103x1579.txt" 1 7 5 1 -1 1 5 1
 * example: ./IDASH2017 "../data/1_training_data_csv" 1 7 5 1 -1 1 0 1 "../data/1_testing_data_csv"
 *
 * parameters:
 * trainfile - path to train file
 * isYfirst - {0,1} y parameter first OR last
 * numIter - number of iterations
 * kdeg - degree of sigmoid approximation function k in {3,5,7}
 * gammaUp - corresponds to learning rate
 * gammaDown - corresponds to learning rate
 * isInitZero - is initial weights zero or average
 * fold - folding method if arguments <= 8 we use folding method
 * isEncrypted - encrypted or plain
 * testfile - path to test file (checks if number of arguments > 8 then we use standard method
 *
 * current files that in data folder (filename isYfirst):
 * "../data/data5x500.txt" false
 * "../data/data9x1253.txt" false
 * "../data/data15x1500.txt" false
 * "../data/data16x101.txt" false
 * "../data/data27x148.txt" false
 * "../data/data51x653.txt" false
 * "../data/data67x216.txt" false
 * "../data/data103x1579.txt" true
 * "../data/1_training_data.csv" true
 * "../data/1_testing_data.csv" true
 *
 * FYI: approx 3 suggested iter: 4, 9, 18, 36, ...
 * FYI: approx 5 suggested iter: 3, 7, 14, 28, ...
 * FYI: approx 7 suggested iter: 3, 7, 14, 28, ...
 */


int main(int argc, char **argv) {
	//	size_t currentAfterSchemeSize = getCurrentRSS( ) >> 20;
	//	size_t peakAfterSchemeSize = getPeakRSS() >> 20;
	//	cout << "Current Memory Usage After Scheme Generation: " << currentAfterSchemeSize << "MB"<< endl;
	//	cout << "Peak Memory Usage After Scheme Generation: " << peakAfterSchemeSize << "MB"<< endl;
	SetNumThreads(8);
	string trainfile(argv[1]);
	bool isYfirst = atoi(argv[2]);
	long numIter = atol(argv[3]);
	long kdeg = atol(argv[4]);
	double gammaUp = atof(argv[5]);
	double gammaDown = atof(argv[6]);
	bool isInitZero = atoi(argv[7]);
	long fold = atol(argv[8]);
	bool isEncrypted = atoi(argv[9]);
	string testfile = argc > 10? string(argv[10]) : trainfile;
	//-----------------------------------------

	if(argc > 10) {
		long factorDimTrain = 0, sampleDimTrain = 0, factorDimTest = 0, sampleDimTest = 0;
		double** zDataTrain = GD::zDataFromFile(trainfile, factorDimTrain, sampleDimTrain, isYfirst);
		double** zDataTest = GD::zDataFromFile(testfile, factorDimTest, sampleDimTest, isYfirst);
		if(factorDimTest != factorDimTrain) {
			invalid_argument("factor dimensions of learn and test datasets do not match");
		}
		if(isEncrypted) {
			TestGD::testEncNLGD(zDataTrain, zDataTest, factorDimTrain, sampleDimTrain, sampleDimTest,
					isYfirst, numIter, kdeg, gammaUp, gammaDown, isInitZero);
		} else {
			TestGD::testPlainNLGD(zDataTrain, zDataTest, factorDimTrain, sampleDimTrain, sampleDimTest,
					isYfirst, numIter, kdeg, gammaUp, gammaDown, isInitZero);
		}
	} else {
		long sampleDim = 0, factorDim = 0;
		double** zData = GD::zDataFromFile(trainfile, factorDim, sampleDim, isYfirst);
		GD::shuffleZData(zData, factorDim, sampleDim);
		if(isEncrypted) {
			TestGD::testEncNLGDFOLD(fold, zData, factorDim, sampleDim, isYfirst, numIter, kdeg, gammaUp, gammaDown, isInitZero);
		} else {
			TestGD::testPlainNLGDFOLD(fold, zData, factorDim, sampleDim, isYfirst, numIter, kdeg, gammaUp, gammaDown, isInitZero);
		}
	}


	return 0;
}
