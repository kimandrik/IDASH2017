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

long** GD::xyDataFromFile(string& path, long& factorDim, long& sampleDim, bool isfirst) {
	vector<vector<long>> xdata;
	vector<long> ydata;
	factorDim = 0; 		// dimension of x
	sampleDim = 0;	// number of samples
	ifstream openFile(path.data());
	if(openFile.is_open()) {
		string line;
		getline(openFile, line);
		for(long i = 0; i < line.length(); ++i) if(line[i] == ',' ) factorDim++;
		if(isfirst) {
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
			while(getline(openFile, line)){
				if(line.length() != 2 * factorDim + 1 ) {
					cout << "Error: data format" << endl;
					break;
				}
				if(line[2 * factorDim] == '0') {
					ydata.push_back(-1);
				} else if(line[2 * factorDim] == '1') {
					ydata.push_back(1);
				} else {
					cout << "Error: data value" << endl;
					break;
				}
				vector<long> vecline;
				for(long i = 0; i < 2 * factorDim; i += 2) {
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
		}
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

void GD::stepNLGD(long**& xyData, double*& wData, double*& vData, long& factorDim, long& learnDim, double& gamma, double& eta) {
	double* grad = new double[factorDim]();

	for(int j = 0; j < learnDim; ++j) {
		double ip = innerprod(wData, xyData[j], factorDim);
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
	debugcheck("d wData: ", wData, 5);
}

void GD::decStepNLGD(long**& xyData, double*& wData, double*& vData, long& factorDim, long& learnDim, double& gamma, double& eta) {

	double** dprod = new double*[learnDim];
	for (long j = 0; j < learnDim; ++j) {
		dprod[j] = new double[factorDim];
	}

	for (long j = 0; j < learnDim; ++j) {
		for (long i = 0; i < factorDim; ++i) {
			dprod[j][i] = xyData[j][i] * wData[i];
		}
	}

	double* dip = new double[learnDim]();
	for (long j = 0; j < learnDim; ++j) {
		for (long i = 0; i < factorDim; ++i) {
			dip[j] += dprod[j][i];
		}
	}

	double** dpows = new double*[3];
	for (long l = 0; l < 3; ++l) {
		dpows[l] = new double[learnDim];
		for (long j = 0; j < learnDim; ++j) {
			dpows[l][j] = pow(dip[j], (double)(1<<l));
		}
	}
	double* coeffs = new double[8]{-0.5,0.216884,0,0.00819276,0,0.000165861,0,-0.00000119581};

	double** dgrad = new double*[learnDim];
	for (long j = 0; j < learnDim; ++j) {
		dgrad[j] = new double[factorDim];
		for (long i = 0; i < factorDim; ++i) {
			dgrad[j][i] = xyData[j][i] * gamma * coeffs[0];
		}
	}

	for (long t = 1; t < 8; t=t+2) {
		for (long j = 0; j < learnDim; ++j) {
			for (long i = 0; i < factorDim; ++i) {
				double dgradit = xyData[j][i] * gamma * coeffs[t];
				if(bit(t, 0)) {
					dgradit *= dpows[0][j];
				}
				if(bit(t, 1)) {
					dgradit *= dpows[1][j];
				}
				if(bit(t, 2)) {
					dgradit *= dpows[2][j];
				}
				dgrad[j][i] += dgradit;
			}
		}
	}

	double* dgradsum = new double[factorDim]();
	for (long j = 0; j < learnDim; ++j) {
		for (long i = 0; i < factorDim; ++i) {
			dgradsum[i] += dgrad[j][i];
		}
	}


	for (long i = 0; i < factorDim; ++i) {
		dgradsum[i] = vData[i] - dgradsum[i];
		wData[i] = vData[i] - dgradsum[i];
		wData[i] *= eta;
		vData[i] = dgradsum[i];
		wData[i] += vData[i];
	}
	debugcheck("d wData: ", wData, 5);
}


double* GD::wsum(double**& wData, long& factorDim, long& wBatch) {
	double* w = new double[factorDim]();

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

void GD::debugcheck(string prefix, double*& w, long factorCheck) {
	cout << prefix;
	for (long i = 0; i < factorCheck; ++i) {
		cout << w[i] << ",";
	}
	cout << endl;
	cout << endl;
}

void GD:: debugcheck(string prefix, double**& w, long factorCheck, long slotCheck) {
	cout << prefix;
	for (long i = 0; i < factorCheck; ++i) {
		for (long j = 0; j < slotCheck; ++j) {
			cout << w[j][i] << ",";
		}
		cout << endl;
	}
	cout << endl;
}

