#include "GD.h"

#include <EvaluatorUtils.h>
#include <CZZ.h>
#include <NTL/BasicThreadPool.h>
#include <NTL/ZZ.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

long** GD::xyDataFromFile(string& path, long& factorDim, long& sampleDim) {
	vector<vector<long>> xdata;
	vector<long> ydata;
	factorDim = 0; 		// dimension of x
	sampleDim = 0;	// number of samples
	ifstream openFile(path.data());
	if(openFile.is_open()) {
		string line;
		getline(openFile, line);
		for(long i = 0; i < line.length(); ++i) if(line[i] == ',' ) factorDim++;
		while(getline(openFile, line)){
			if(line.length() != 2 * factorDim + 1 ) {
				cout << "Error: data format" << endl;
				break;
			}
			if(line[0] == '0') {
				ydata.push_back(-1);
			} else if(line[0] == '1') {
				ydata.push_back(1);
			} else {
				cout << "Error: data value" << endl;
				break;
			}
			vector<long> vecline;
			for(long i = 2; i < 2 * factorDim + 1; i += 2) {
				if(line[i] == '0') {
					vecline.push_back(0);
				} else if(line[i] == '1') {
					vecline.push_back(1);
				} else{
					cout << "Error: data value" << endl;
					break;
				}
			}
			xdata.push_back(vecline);
			sampleDim++;
		}
		openFile.close();
	} else {
		cout << "Error: cannot read file" << endl;
	}
	long** zdata = new long*[sampleDim];
	for(int j = 0; j < sampleDim; ++j){
		long* zj = new long[factorDim];
		for(int i = 0; i < factorDim; ++i){
			zj[i] = ydata[j] * xdata[j][i];
		}
		zdata[j] = zj;
	}
	return zdata;
}

double GD::innerprod(double*& w, long*& xy, long& size){
	double res = 0.0;
	for(int i = 0; i < size; ++i) {
		res += w[i] * xy[i];
	}
	return res;
}

void GD::stepQGD(long**& xyData, double*& wData, long& factorDim, long& learnDim, double& lambda, double& gamma) {
	double* grad = new double[factorDim];
	for(int i = 0; i < factorDim; ++i) {
		grad[i] = lambda * wData[i];
	}

	for(int j = 0; j < learnDim; ++j) {
		double ip = innerprod(wData, xyData[j], factorDim);
		double tmp = 2.0 * (ip - 1.0) / learnDim;
		for(int i = 0; i < factorDim; ++i) {
			grad[i] += tmp * (double) xyData[j][i];
		}
	}
	for (int i = 0; i < factorDim; ++i) {
		wData[i] -= gamma * grad[i];
	}
}

void GD::stepSQGD(long**& xyData, double*& wData, long& factorDim, long& learnDim, double& lambda, double& gamma, long& stochDim) {
	double* grad = new double[factorDim];
	for(int i = 0; i < factorDim; ++i) {
		grad[i] = lambda * wData[i];
	}

	for(int j = 0; j < stochDim; ++j) {
		long rnd = RandomBnd(learnDim);
		double ip = innerprod(wData, xyData[rnd], factorDim);
		double tmp = 2.0 * (ip - 1.0) / stochDim;
		for(int i = 0; i < factorDim; ++i) {
			grad[i] += tmp * (double) xyData[rnd][i];
		}
	}
	for (int i = 0; i < factorDim; ++i) {
		wData[i] -= gamma * grad[i];
	}
}

void GD::stepLGD(long**& xyData, double*& wData, long& factorDim, long& learnDim, double& lambda, double& gamma) {
	double* grad = new double[factorDim];
	for(int i = 0; i < factorDim; ++i) {
		grad[i] = lambda * wData[i];
	}

	for(int j = 0; j < learnDim; ++j) {
		double ip = innerprod(wData, xyData[j], factorDim);
		double tmp = (ip > 15.0) ? 0 : (ip < -15.0) ? -1.0 : - 1. / (1. + exp(ip));
		tmp /= (double)learnDim;
		for(int i = 0; i < factorDim; ++i) {
			grad[i] += tmp * (double) xyData[j][i];
		}
	}
	for (int i = 0; i < factorDim; ++i) {
		wData[i] -= gamma * grad[i];
	}
}

void GD::stepSLGD(long**& xyData, double*& wData, long& factorDim, long& learnDim, double& lambda, double& gamma, long& stochDim) {
	double* grad = new double[factorDim];
	for(int i = 0; i < factorDim; ++i) {
		grad[i] = lambda * wData[i];
	}

	for(int j = 0; j < stochDim; ++j) {
		long rnd = RandomBnd(learnDim);
		double ip = innerprod(wData, xyData[rnd], factorDim);
		double tmp = (ip > 15.0) ? 0 : (ip < -15.0) ? -1.0 : - 1. / (1. + exp(ip));
		for(int i = 0; i < factorDim; ++i) {
			grad[i] += tmp * (double) xyData[rnd][i];
		}
	}

	for (int i = 0; i < factorDim; ++i) {
		wData[i] -= gamma * grad[i];
	}
}

void GD::stepMLGD(long**& xyData, double*& wData, double*& vData, long& factorDim, long& learnDim, double& lambda, double& gamma, double& eta) {
	double* grad = new double[factorDim];
	for(int i = 0; i < factorDim; ++i) {
		grad[i] = lambda * wData[i];
	}

	for(int j = 0; j < learnDim; ++j) {
		double ip = innerprod(wData, xyData[j], factorDim);
		double tmp = (ip > 15.0) ? 0 : (ip < -15.0) ? -1.0 : - 1. / (1. + exp(ip));
		for(int i = 0; i < factorDim; ++i) {
			grad[i] += tmp * (double) xyData[j][i];
		}
	}

	for (int i = 0; i < factorDim; ++i) {
		vData[i] *= eta;
		vData[i] += gamma * grad[i];
		wData[i] -= vData[i];
	}
}

void GD::stepNLGD(long**& xyData, double*& wData, double*& vData, long& factorDim, long& learnDim, double& lambda, double& gamma, double& eta) {
	double* grad = new double[factorDim];
	for(int i = 0; i < factorDim; ++i) {
		grad[i] = lambda * wData[i];
	}

	for(int j = 0; j < learnDim; ++j) {
		double ip = innerprod(wData, xyData[j], factorDim);
//		double tmp = 0;
//		if(ip > 5) {
//			cout << ip << endl;
//			tmp = 0.0;
//		} else if(tmp < -5)	{
//			cout << ip << endl;
//			tmp = -1.0;
//		} else {
//			tmp = - 1. / (1. + exp(ip));
//		}

		double tmp = (ip > 15.0) ? 0 : (ip < -15.0) ? -1.0 : - 1. / (1. + exp(ip));
		for(int i = 0; i < factorDim; ++i) {
			grad[i] += tmp * (double) xyData[j][i];
		}
	}

	// Nesterov steps
	for (int i = 0; i < factorDim; ++i) {
		double tmpv = wData[i] - gamma * grad[i];
		wData[i] = (1.0 - eta) * tmpv + eta * vData[i];
		vData[i] = tmpv;
	}
}

double* GD::wsum(double**& wData, long& factorDim, long& wBatch) {
	double* w = new double[factorDim];

	for (long i = 0; i < factorDim; ++i) {
		for (int l = 0; l < wBatch; ++l) {
			w[i] += wData[l][i];
		}
	}
	return w;
}

void GD::check(long**& xyData, double*& w, long& factorDim, long& sampleDim) {
	cout << "w:";
	for (long i = 0; i < factorDim; ++i) {
		cout << w[i] << ",";
	}
	cout << endl;

	long num = 0;
	for(long j = 0; j < sampleDim; ++j){
		if(innerprod(w, xyData[j], factorDim) > 0) num++;
	}
	cout << "Correctness: " << num << "/" << sampleDim << endl;

}

void GD::debugcheck(string prefix, double*& w, long& factorDim) {
	cout << prefix;
	for (long i = 0; i < factorDim; ++i) {
		cout << w[i] << ",";
	}
	cout << endl;
}

void GD:: debugcheck(string prefix, long*& xy, long& factorDim) {
	cout << prefix;
	for (long i = 0; i < factorDim; ++i) {
		cout << xy[i] << ",";
	}
	cout << endl;
}

