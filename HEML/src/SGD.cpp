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

double SGD::innerprod(double*& wdata, long*& x, long& size){
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
		double ip = ip(wdata, zdata[i], dim);
		double tmp = plainphiprime(ip);
		for(int j = 0; j < dim; ++j) {
			grad[j] += tmp * (double) zdata[i][j];
		}
	}
	return grad;
}

double** SGD::wdatagen(long& wnum, long& dim) {
	double** wdata = new double*[wnum];
	for (long l = 0; l < wnum; ++l) {
		wdata[l] = new double[dim];
		for (long i = 0; i < dim; ++i) {
			wdata[l][i] = 1.0 - (double)rand() / RAND_MAX; // change to good initial w choice
		}
	}
	return wdata;
}

double* SGD::gammagen(long& iter) {
	double* gamma = new double[iter];
	for (long k = 0; k < iter; ++k) {
		gamma[k] = 1.0 / (k + 1);
	}
	return gamma;
}

void SGD::step(double**& wdata, long**& zdata, double& gamma, double& lambda, long& wnum, long& dim, long& sampledim) {
	for (long l = 0; l < wnum; ++l) {
		double* grad = new double[dim];
		for(int i = 0; i < dim; ++i) {
			grad[i] = lambda * wdata[i];
		}

		for(int i = 0; i < sampledim; ++i) {
			double ip = ip(wdata, zdata[i], dim);
			double tmp = - 1. / (1. + exp(ip));
			for(int j = 0; j < dim; ++j) {
				grad[j] += tmp * (double) zdata[i][j];
			}
		}
		for (int i = 0; i < dim; ++i) {
			wdata[l][i] -= gamma * grad[i];
		}
	}
}

double* SGD::wgen(double**& wdata, long& wnum, long& dim) {
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
	for(long i = 0; i < sampledim; ++i){
		if(innerprod(w, zdata[i], dim) > 0) num++;
	}
	cout << "Correctness: " << num << "/" << sampledim << endl;

}

Cipher* SGD::enczdata(long**& zdata, long& slots, long& wnum, long& dim, long& sampledim, ZZ& p) {
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

Cipher* SGD::encwdata(double**& wdata, long& slots, long& wnum, long& dim, long& sampledim, long& logp) {
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

ZZ* SGD::pgammagen(double*& alpha, long& iter, long& logp) {
	ZZ* palpha = new ZZ[iter];
	for (long k = 0; k < iter; ++k) {
		RR ralpha = to_RR(alpha[k]);
		RR pralpha = MakeRR(ralpha.x, ralpha.e + logp);
		palpha[k] = to_ZZ(pralpha);
	}
	return palpha;
}

void SGD::encStep(Cipher*& czdata, Cipher*& cwdata, ZZ& pgamma, long& lambda, long& slots, long& wnum, long& dim, long& sampledim) {

	Cipher cip = algo.innerProd(czdata, cwdata, dim); // ip (-1)

	Cipher* cpows = algo.powerOf2Extended(cip, 2); // ip (-1), ip^2 (-2), ip^4 (-3)

	ZZ* pows = scheme.aux.taylorPowsMap.at(SIGMOIDBARGOOD);
	ZZ* wcnst = scheme.params.p - pgamma * lambda;
	for (long i = 0; i < dim; ++i) {
		scheme.multByConstAndEqual(cwdata[i], wcnst); // (1 - gamma * lambda) w
		scheme.modSwitchOneAndEqual(cwdata[i]); // (1 - gamma * lambda) w  (-1)
	}

	Cipher* cgrad = new Cipher[dim];
	for (long t = 0; t < 8; ++t) {
		ZZ* cnst = pgamma * pows[t] / sampledim / scheme.params.p; // p * (gamma * alpha_t / m)
		NTL_EXEC_RANGE(dim, first, last);
		for (long i = first; i < last; ++i) {
			if(cnst != ZZ::zero()) {
				Cipher cgradit = scheme.multByConst(czdata[i], cnst); // p * p * (z * gamma * alpha_t / m)
				scheme.modSwitchOneAndEqual(cgradit); // p * gamma * 1/m * (z * alpha_t) (-1)
				for (int b = 0; b < 3; ++b) {
					if(bit(t, b)) {
						scheme.modEmbedAndEqual(cgradit, cpows[b].level);
						scheme.multModSwitchOneAndEqual(cgradit, cpows[b]);
					}
				}
				scheme.modEmbedAndEqual(cgrad[i], cgradit.level);
				scheme.addAndEqual(cgrad[i], cgradit); // p * gamma * 1/m * (z * sigmoid(ip)) (-4)
			}
		}
		NTL_EXEC_RANGE_END;
	}

	long logslots = log2(slots);
	long logwnum = log2(wnum);

	NTL_EXEC_RANGE(dim, first, last);
	for (long i = first; i < last; ++i) {
		for (long i = logwnum; i < logslots; ++i) {
			Cipher rot = scheme.leftRotateByPo2(cgrad[i], i);
			scheme.addAndEqual(cgrad[i], rot); // p * gamma * 1/m * sum (z * sigmoid(ip)) (-4)
		}
		scheme.modEmbedAndEqual(cwdata[i], cgrad[i].level);
		scheme.subAndEqual(cwdata[i], cgrad[i]); // w - gamma * grad(w) - gamma * lambda * w (-4)
	}
	NTL_EXEC_RANGE_END;
}

Cipher* SGD::encwgen(Cipher*& cwdata, long& wnum, long& dim) {
	Cipher* cw = new Cipher[dim];
	for (long i = 0; i < dim; ++i) {
		cw[i] = algo.partialSlotsSum(cwdata[i], wnum);
	}
	return cw;
}

double* SGD::decw(SecKey& secretKey, Cipher*& cw, long& dim) {
	double* w = new double[dim];
	for (long i = 0; i < dim; ++i) {
		CZZ* dcw = scheme.decrypt(secretKey, cw[i]);
		w[i] = dcw[0];
	}
	return w;
}
