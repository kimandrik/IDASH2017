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

long** SGD::zdataFromFile(string& path, long& dim, long& sampledim) {
	vector<vector<long>> xdata;
	vector<long> ydata;
	dim = 0; 		// dimension of x
	sampledim = 0;	// number of samples
	ifstream openFile(path.data());
	if(openFile.is_open()) {
		string line;
		getline(openFile, line);
		for(long i = 0; i < line.length(); ++i) if(line[i] == ',' ) dim++;
		while(getline(openFile, line)){
			if(line.length() != 2 * dim + 1 ) {
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
			for(long i = 2; i < 2 * dim + 1; i += 2) {
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
			sampledim++;
		}
		openFile.close();
	} else {
		cout << "Error: cannot read file" << endl;
	}
	long** zdata = new long*[sampledim];
	for(int j = 0; j < sampledim; ++j){
		long* zj = new long[dim];
		for(int i = 0; i < dim; ++i){
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

double** SGD::wdatagen(long& wnum, long& dim) {
	double** wdata = new double*[wnum];
	for (long l = 0; l < wnum; ++l) {
		wdata[l] = new double[dim];
		for (long i = 0; i < dim; ++i) {
			wdata[l][i] = (1.0 - 2.0 * (double)rand() / RAND_MAX) / 32.0; // change to good initial w choice
		}
	}
	return wdata;
}

double* SGD::gammagen(long& iter) {
	double* gamma = new double[iter];
	for (long k = 0; k < iter; ++k) {
		gamma[k] = 0.1 / (k + 1);
	}
	return gamma;
}

void SGD::steplogregress(double*& wdata, long**& zdata, double& gamma, double& lambda, long& dim, long& sampledim) {
	double* grad = new double[dim];
	for(int i = 0; i < dim; ++i) {
		grad[i] = lambda * wdata[i];
	}

	for(int j = 0; j < sampledim; ++j) {
		double ip = innerprod(wdata, zdata[j], dim);
		double tmp = - 1. / (1. + exp(ip));
		for(int i = 0; i < dim; ++i) {
			grad[i] += tmp * (double) zdata[j][i];
		}
	}
	for (int i = 0; i < dim; ++i) {
		wdata[i] -= gamma * grad[i];
	}
}

void SGD::stepsimpleregress(double*& wdata, long**& zdata, double& gamma, double& lambda, long& dim, long& sampledim) {
	double* grad = new double[dim];
	for(int i = 0; i < dim; ++i) {
		grad[i] = lambda * wdata[i];
	}

	for(int j = 0; j < sampledim; ++j) {
		double ip = innerprod(wdata, zdata[j], dim);
		double tmp = (ip - 2) / dim;
		for(int i = 0; i < dim; ++i) {
			grad[i] += tmp * (double) zdata[j][i];
		}
	}
	for (int i = 0; i < dim; ++i) {
		wdata[i] -= gamma * grad[i];
	}
}

double* SGD::wout(double**& wdata, long& wnum, long& dim) {
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

