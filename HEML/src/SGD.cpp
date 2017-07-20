#include "SGD.h"

#include <CZZ.h>
#include <NTL/BasicThreadPool.h>
#include <NTL/ZZ.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

long** SGD::dataFromFile(string& path, long& dim, long& sampledim) {
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
	for(int i = 0; i < sampledim; ++i){
		long* zi = new long[dim];
		for(int j = 0; j < dim; ++j){
			zi[j] = ydata[i] * xdata[i][j];
		}
		zdata[i] = zi;
	}
	return zdata;
}

double SGD::plaincost(double*& w, long**& data, long& dim, long& sampledim){
	double res = 0;
	for(int i = 0; i < dim; ++i){
		res += 0.5 * w[i] * w[i];
	}
	for(int i = 0; i < sampledim; ++i){
		double tmp = plainip(w, data[i], dim);
		res += plainphi(tmp);
	}
	return res;
}

double SGD::plainnorm(double*& w, long& size){
	double res = 0;
	for(int i = 0; i < size; ++i) res += w[i] * w[i];
	return sqrt(res);
}

double SGD::plainphi(double& b){
	double res = log(1.0 + exp(-b));
	return res;
}

double SGD::plainphiprime(double& b) {
	double res = - 1. / (1. + exp(b));
	return res;
}

double SGD::plainip(double*& wdata, long*& x, long& size){
	double res = 0.0;
	for(int i = 0; i < size; ++i) {
		res += wdata[i] * x[i];
	}
	return res;
}

double* SGD::plainGradient(double*& wdata, long**& zdata, long& dim, long& sampledim, double& lambda) {
	double* grad = new double[dim];
	for(int i = 0; i < dim; ++i) {
		grad[i] = lambda * wdata[i];
	}

	for(int i = 0; i < sampledim; ++i) {
		double ip = plainip(wdata, zdata[i], dim);
		double tmp = plainphiprime(ip);
		for(int j = 0; j < dim; ++j) {
			grad[j] += tmp * (double) zdata[i][j];
		}
	}
	return grad;
}

void SGD::plainsgd(long& iter, double*& wdata, long**& zdata, long& dim, long& sampledim) {

	for (long k = 0; k < iter; ++k) {
		double alpha = 1.0 / (k+1);
		double lambda = 2.0;
		double* grad = plainGradient(wdata, zdata, dim, sampledim, lambda);
		for (int i = 0; i < dim; ++i) {
			wdata[i] -= alpha * grad[i];
		}
		if((k+1) % (iter/5) == 0) {
			double c = plaincost(wdata, zdata, dim, sampledim);
			cout << k + 1  << "-th. ||w|| = " << plainnorm(wdata, dim) << ", ||grad|| = " << plainnorm(grad, dim) <<".\n";
			cout << "cost: " << c << "\n\n";
			cout << "w: ";
			for (long i = 0; i < dim; ++i) {
				cout << wdata[i] << ",";
			}
			cout << endl;
		}
	}

	long num = 0;
	for(long i = 0; i < sampledim; ++i){
		if(plainip(wdata, zdata[i], dim) > 0) num++;
	}
	cout << "Correctness: " << num << "/" << sampledim << endl;

}

Cipher* SGD::cipherGradient(Cipher*& zciphers, Cipher*& wciphers, const long& dim, const long& slots, const long& wnum) {
	Cipher ip = algo.innerProd(zciphers, wciphers, dim);
	scheme.doubleAndEqual(ip);
	Cipher sig =  algo.function(ip, SIGMOIDBARGOOD, 7); // 3 levels sig_i = sigmoid_i (z_i1 * w_1 + z_i2 * w_2 + ... + z_in * w_n) * p

	Cipher* grad = new Cipher[dim];

	NTL_EXEC_RANGE(dim, first, last);
	for (long i = first; i < last; ++i) {
		grad[i] = scheme.modEmbed(zciphers[i], sig.level);
		scheme.multModSwitchOneAndEqual(grad[i], sig); // 1 level grad_i = sigmoid_j ( (xj1w1 + xj2w2 + ... + xjnwn) * y_j ) * z_ij * p
		long logslots = log2(sampledim);
		long logwnum = log2(wnum);
		for (long i = logwnum; i < logslots; ++i) {
			Cipher rot = scheme.leftRotateByPo2(grad[i], i);
			scheme.addAndEqual(grad[i], rot);
		}
	}
	NTL_EXEC_RANGE_END;

	return grad; // 7 levels total
}
