#include "SGD.h"

#include <EvaluatorUtils.h>
#include <CZZ.h>
#include <NTL/BasicThreadPool.h>
#include <NTL/ZZ.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

long** SGD::xyDataFromFile(string& filename, long& factorDim, long& sampleDim) {
	vector<vector<long>> xdata;
	vector<long> ydata;
	factorDim = 0; 		// dimension of x
	sampleDim = 0;	// number of samples
	ifstream openFile(filename.data());
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

double SGD::innerprod(double*& wdata, long*& x, long& size){
	double res = 0.0;
	for(int i = 0; i < size; ++i) {
		res += wdata[i] * x[i];
	}
	return res;
}

void SGD::stepQuadraticRegress(double*& wdata, long**& zdata, double& gamma, double& lambda, long& dim, long& learndim) {
	double* grad = new double[dim];
	for(int i = 0; i < dim; ++i) {
		grad[i] = lambda * wdata[i];
	}

	for(int j = 0; j < learndim; ++j) {
		double ip = innerprod(wdata, zdata[j], dim);
		double tmp = 2.0 * (ip - 1.0) / learndim;
		for(int i = 0; i < dim; ++i) {
			grad[i] += tmp * (double) zdata[j][i];
		}
	}
	for (int i = 0; i < dim; ++i) {
		wdata[i] -= gamma * grad[i];
	}
}


void SGD::stepLogRegress(double*& wdata, long**& zdata, double& gamma, double& lambda, long& dim, long& learndim) {
	double* grad = new double[dim];
	for(int i = 0; i < dim; ++i) {
		grad[i] = lambda * wdata[i];
	}

	for(int j = 0; j < learndim; ++j) {
		double ip = innerprod(wdata, zdata[j], dim);
		double tmp = - 1. / (1. + exp(ip));
		tmp /= (double)learndim;
		for(int i = 0; i < dim; ++i) {
			grad[i] += tmp * (double) zdata[j][i];
		}
	}
	for (int i = 0; i < dim; ++i) {
		wdata[i] -= gamma * grad[i];
	}
}

void SGD::stepStochasticQuadraticRegress(double*& wdata, long**& zdata, double& gamma, double& lambda, long& dim, long& learndim, long& stochdim) {
	double* grad = new double[dim];
	for(int i = 0; i < dim; ++i) {
		grad[i] = lambda * wdata[i];
	}

	for(int j = 0; j < stochdim; ++j) {
		long rnd = RandomBnd(learndim);
		double ip = innerprod(wdata, zdata[rnd], dim);
		double tmp = 2.0 * (ip - 1.0) / stochdim;
		for(int i = 0; i < dim; ++i) {
			grad[i] += tmp * (double) zdata[rnd][i];
		}
	}
	for (int i = 0; i < dim; ++i) {
		wdata[i] -= gamma * grad[i];
	}
}

void SGD::stepStochasticLogRegress(double*& wdata, long**& zdata, double& gamma, double& lambda, long& dim, long& learndim, long& stochdim) {
	double* grad = new double[dim];
	for(int i = 0; i < dim; ++i) {
		grad[i] = lambda * wdata[i];
	}

	for(int j = 0; j < stochdim; ++j) {
		long rnd = RandomBnd(learndim);
		double ip = innerprod(wdata, zdata[rnd], dim);
		double tmp = - 1. / (1. + exp(ip)) / stochdim;
		for(int i = 0; i < dim; ++i) {
			grad[i] += tmp * (double) zdata[rnd][i];
		}
	}

	for (int i = 0; i < dim; ++i) {
		wdata[i] -= gamma * grad[i];
	}
}

void SGD::stepMomentumLogRegress(double*& wdata, double*& vdata, long**& zdata, double& gamma, long& dim, long& learndim, double& eta) {
	double* grad = new double[dim]();

	for(int j = 0; j < learndim; ++j) {
		double ip = innerprod(wdata, zdata[j], dim);
		double tmp = - 1. / (1. + exp(ip));
		for(int i = 0; i < dim; ++i) {
			grad[i] += tmp * (double) zdata[j][i];
		}
	}

	for (int i = 0; i < dim; ++i) {
		vdata[i] *= eta;
		vdata[i] += gamma * grad[i];
		wdata[i] -= vdata[i];
	}
}

void SGD::stepNesterovLogRegress(double*& wdata, double*& vdata, long**& xtData, double& gamma, long& dim, long& learndim, double& eta) {
	double* grad = new double[dim]();

	for(int j = 0; j < learndim; ++j) {
		double ip = innerprod(wdata, xtData[j], dim);
		double tmp = - 1. / (1. + exp(ip));
		for(int i = 0; i < dim; ++i) {
			grad[i] += tmp * (double) xtData[j][i];
		}
	}

	for (int i = 0; i < dim; ++i) {
		double tmpv = wdata[i] - gamma * grad[i];
		wdata[i] = (1.0 - eta) * tmpv + eta * vdata[i];
		vdata[i] = tmpv;
	}
}

double* SGD::waverage(double**& wdata, long& wnum, long& dim) {
	double* w = new double[dim];

	for (long i = 0; i < dim; ++i) {
		for (int l = 0; l < wnum; ++l) {
			w[i] += wdata[l][i];
		}
		w[i] /= wnum;
	}
	return w;
}

void SGD::check(double*& w, long**& zdata, long& dim, long& sampledim) {
	cout << "w:";
	for (long i = 0; i < dim; ++i) {
		cout << w[i] << ",";
	}
	cout << endl;

	long num = 0;
	for(long j = 0; j < sampledim; ++j){
		if(innerprod(w, zdata[j], dim) > 0) num++;
	}
	cout << "Correctness: " << num << "/" << sampledim << endl;

}

void SGD::debugcheck(string prefix, double*& w, long& dim) {
	cout << prefix;
	for (long i = 0; i < dim; ++i) {
		cout << w[i] << ",";
	}
	cout << endl;
}

void SGD:: debugcheck(string prefix, long*& z, long& dim) {
	cout << prefix;
	for (long i = 0; i < dim; ++i) {
		cout << z[i] << ",";
	}
	cout << endl;
}

