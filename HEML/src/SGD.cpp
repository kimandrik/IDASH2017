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
	for(int i = 0; i < sampledim; ++i){
		long* zi = new long[dim];
		for(int j = 0; j < dim; ++j){
			zi[j] = ydata[i] * xdata[i][j];
		}
		zdata[i] = zi;
	}
	return zdata;
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

double* SGD::sgd(long& iter, double**& wdata, long**& zdata, double*& alpha, double& lambda, long& wnum, long& dim, long& sampledim) {
	for (long l = 0; l < wnum; ++l) {
		for (long k = 0; k < iter; ++k) {
			double* grad = plainGradient(wdata[l], zdata, dim, sampledim, lambda);
			for (int i = 0; i < dim; ++i) {
				wdata[l][i] -= alpha[k] * grad[i];
			}
		}
	}
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
	long num = 0;
	for(long i = 0; i < sampledim; ++i){
		if(plainip(w, zdata[i], dim) > 0) num++;
	}
	cout << "Correctness: " << num << "/" << sampledim << endl;

}

Cipher* SGD::cipherGradient(Cipher*& czdata, Cipher*& cwdata, long& slots, long& wnum, long& dim) {
	Cipher* grad = new Cipher[dim];

	Cipher ip = algo.innerProd(czdata, cwdata, dim);
	scheme.doubleAndEqual(ip);
	Cipher sig =  algo.function(ip, SIGMOIDBARGOOD, 7); // 3 levels sig_i = sigmoid_i (z_i1 * w_1 + z_i2 * w_2 + ... + z_in * w_n) * p
	NTL_EXEC_RANGE(dim, first, last);
	for (long i = first; i < last; ++i) {
		grad[i] = scheme.modEmbed(czdata[i], sig.level);
		scheme.multModSwitchOneAndEqual(grad[i], sig); // 1 level grad_i = sigmoid_j ( (xj1w1 + xj2w2 + ... + xjnwn) * y_j ) * z_ij * p
		long logslots = log2(slots);
		long logwnum = log2(wnum);
		for (long i = logwnum; i < logslots; ++i) {
			Cipher rot = scheme.leftRotateByPo2(grad[i], i);
			scheme.addAndEqual(grad[i], rot);
		}
	}
	NTL_EXEC_RANGE_END;

	return grad; // 7 levels total
}

Cipher* SGD::encryptzdata(long**& zdata, long& slots, long& wnum, long& dim, long& sampledim, ZZ& p) {
	Cipher* czdata = new Cipher[dim];
	for (long i = 0; i < dim; ++i) {
		CZZ* pzdata = new CZZ[slots];
		for (long j = 0; j < sampledim; ++j) {
			if(zdata[i][j] == -1) {
				for (long l = 0; l < wnum; ++l) {
					pzdata[wnum * j + l] = CZZ(-p);
				}
			} else if(zdata[i][j] == 1) {
				for (long l = 0; l < wnum; ++l) {
					pzdata[wnum * j + l] = CZZ(p);
				}
			}
		}
		cout << i << endl;
		czdata[i] = scheme.encrypt(pzdata, slots);
		cout << i << endl;
	}
	return czdata;
}

Cipher* SGD::encryptwdata(double**& wdata, long& slots, long& wnum, long& dim, long& sampledim, long& logp) {
	Cipher* cwdata = new Cipher[dim];
	for (long i = 0; i < dim; ++i) {
		CZZ* pwdata = new CZZ[slots];
		for (long j = 0; j < sampledim; ++j) {
			for (long l = 0; l < wnum; ++l) {
				pwdata[wnum * j + l] = EvaluatorUtils::evaluateVal(wdata[l][i], 0.0, logp);
			}
		}
		cwdata[i] = scheme.encrypt(pwdata, slots);
	}
	return cwdata;
}

Cipher* SGD::ciphersgd(long& iter, Cipher*& czdata, Cipher*& cwdata, ZZ*& palpha, long& slots, long& wnum, long& dim) {
	for (long k = 0; k < iter; ++k) {
		Cipher* cgrad = cipherGradient(czdata, cwdata, slots, wnum, dim);
		NTL_EXEC_RANGE(dim, first, last);
		for (long i = first; i < last; ++i) {
			scheme.multByConstAndEqual(cgrad[i], palpha[k]);
			scheme.modSwitchOneAndEqual(cgrad[i]);
			scheme.modEmbedAndEqual(cwdata[i], cgrad[i].level);
			scheme.subAndEqual(cwdata[i], cgrad[i]);
		}
		NTL_EXEC_RANGE_END;
	}
	Cipher* cw = new Cipher[dim];

	for (long i = 0; i < dim; ++i) {
		cw[i] = algo.partialSlotsSum(cwdata[i], wnum);
	}
	return cw;
}
