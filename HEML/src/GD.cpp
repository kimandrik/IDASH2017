#include "GD.h"

#include <NTL/ZZ.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using namespace NTL;

double** GD::xyDataFromFile(string& path, long& factorDim, long& sampleDim, bool isfirst) {
	vector<vector<double>> xyline;
	factorDim = 1; 	// dimension of x
	sampleDim = 0;	// number of samples
	ifstream openFile(path.data());
	if(openFile.is_open()) {
		string line, temp;
		getline(openFile, line);
		long i;
		size_t start, end;
		for(i = 0; i < line.length(); ++i) if(line[i] == ',' ) factorDim++;

		while(getline(openFile, line)){
			vector<double> vecline;
			do {
				end = line.find_first_of (',', start);
				temp = line.substr(start,end);
				vecline.push_back(atof(temp.c_str()));
				start = end + 1;
			} while(start);
			xyline.push_back(vecline);
			sampleDim++;
		}
	} else {
		cout << "Error: cannot read file" << endl;
	}
	double** xyData = new double*[sampleDim];
	if(isfirst) {
		for(long j = 0; j < sampleDim; ++j){
			double* xyj = new double[factorDim];
			xyj[0] = 2 * xyline[j][0] - 1;
			for(long i = 1; i < factorDim; ++i){
				xyj[i] = xyj[0] * xyline[j][i];
			}
			xyData[j] = xyj;
		}
	} else {
		for(long j = 0; j < sampleDim; ++j){
			double* xyj = new double[factorDim];
			xyj[0] = 2 * xyline[j][factorDim - 1] - 1;
			for(long i = 1; i < factorDim; ++i){
				xyj[i] = xyj[0] * xyline[j][i-1];
			}
			xyData[j] = xyj;
		}
	}
	return xyData;
}

double** GD::RandomxyDataLearn(double** xyData, long learnDim, long sampleDim, long factorDim) {
	double** res = new double*[learnDim];

	bool* notTaken = new bool[sampleDim];
	fill_n(notTaken, sampleDim, false);
	long r;
	for (long i = 0; i < sampleDim - learnDim; ++i) {
		do {
			r = RandomBnd(sampleDim);
		} while(notTaken[r]);
		notTaken[r] = true;
	}
	long idx = 0;
	for (long i = 0; i < sampleDim; ++i) {
		if(!notTaken[i]) {
			res[idx++] = xyData[i];
		}
	}
	delete[] notTaken;
	return res;
}

double GD::innerprod(double* w, double* xy, long size){
	double res = 0.0;
	for(long i = 0; i < size; ++i) {
		res += w[i] * xy[i];
	}
	return res;
}

void GD::stepLGD(double** xyData, double* wData, long factorDim, long learnDim, double& gamma) {
	double* grad = new double[factorDim]();

	for(long j = 0; j < learnDim; ++j) {
		double ip = innerprod(wData, xyData[j], factorDim);
		double tmp = (ip > 15.0) ? 0 : (ip < -15.0) ? -1.0 : - 1. / (1. + exp(ip));
		if(ip > 6) {
			cout << "too big ip: " << ip << endl;
		} else if(ip < -6) {
			cout << "too small ip: " << ip << endl;
		}
		for(int i = 0; i < factorDim; ++i) {
			grad[i] += tmp * (double) xyData[j][i];
		}
	}

	for (int i = 0; i < factorDim; ++i) {
		wData[i] -= gamma * grad[i];
	}
	delete[] grad;
}

void GD::stepMLGD(double** xyData, double* wData, double* vData, long factorDim, long learnDim, double& gamma, double& eta) {
	double* grad = new double[factorDim]();

	for(long j = 0; j < learnDim; ++j) {
		double ip = innerprod(wData, xyData[j], factorDim);
		double tmp = (ip > 15.0) ? 0 : (ip < -15.0) ? -1.0 : - 1. / (1. + exp(ip));
		if(ip > 6) {
			cout << "too big ip: " << ip << endl;
		} else if(ip < -6) {
			cout << "too small ip: " << ip << endl;
		}
		for(long i = 0; i < factorDim; ++i) {
			grad[i] += tmp * (double) xyData[j][i];
		}
	}

	for (int i = 0; i < factorDim; ++i) {
		vData[i] = eta * vData[i] + gamma * grad[i];
		wData[i] -= vData[i];
	}
	delete[] grad;
}

void GD::stepNLGD(double** xyData, double* wData, double* vData, long factorDim, long learnDim, double& gamma, double& eta) {
	double* grad = new double[factorDim]();

	for(long j = 0; j < learnDim; ++j) {
		double ip = innerprod(vData, xyData[j], factorDim);
		double tmp = (ip > 15.0) ? 0 : (ip < -15.0) ? -1.0 : - 1. / (1. + exp(ip));
		if(ip > 6) {
			cout << "too big ip: " << ip << endl;
		} else if(ip < -6) {
			cout << "too small ip: " << ip << endl;
		}
		for(long i = 0; i < factorDim; ++i) {
			grad[i] += tmp * (double) xyData[j][i];
		}
	}

	// Nesterov steps
	for (long i = 0; i < factorDim; ++i) {
		double tmpw = vData[i] - gamma * grad[i];
		vData[i] = (1.0 - eta) * tmpw + eta * wData[i];
		wData[i] = tmpw;
	}
	delete[] grad;
}

void GD::check(double** xyData, double* wData, long factorDim, long sampleDim) {
	cout << "w:";
	for (long i = 0; i < factorDim; ++i) {
		cout << wData[i] << ",";
	}
	cout << endl;

	long num = 0;
	for(long j = 0; j < sampleDim; ++j){
		if(innerprod(wData, xyData[j], factorDim) > 0) num++;
	}
	cout << "Correctness: " << num << "/" << sampleDim << endl;

}

double* GD::calculateYtrueData(double** xyData, long sampleDim) {
	double* res = new double[sampleDim];
	for (long i = 0; i < sampleDim; ++i) {
		res[i] = xyData[i][0] < 0 ? 0 : 1;
	}
	return res;
}

double* GD::calculateYpredictData(double** xyData, double* wData, long factorDim, long sampleDim) {
	double* res = new double[sampleDim];
	for(long i = 0; i < sampleDim; ++i){
		res[i] = innerprod(wData, xyData[i], factorDim) * xyData[i][0] / 2. + 0.5;
	}
	return res;
}

double GD::calcuateAUC(double** xyData, double* wData, long factorDim, long sampleDim, long steps) {

	double* yTrueData = calculateYtrueData(xyData, sampleDim);
	double* yPredictData = calculateYpredictData(xyData, wData, factorDim, sampleDim);

	double* TPR = new double[steps + 1];
	double* FPR = new double[steps + 1];

	for (long i = 0; i < steps + 1; ++i) {
		double threshold =  (double)i / steps;

		long FP = 0;
		long TP = 0;
		long TN = 0;
		long FN = 0;

		for (long j = 0; j < sampleDim; ++j) {
			if(yTrueData[j] == 0) {
				if(yPredictData[j] > threshold) FP++;
				else TN++;
			} else {
				if(yPredictData[j] > threshold) TP++;
				else FN++;
			}
		}
		double TPFN = TP + FN;
		double FPTN = FP + TN;

		TPR[i] = TP / TPFN;
		FPR[i] = FP / FPTN;
	}
	double auc = 0.0;

	for (long i = 0; i < steps; ++i) {
		auc += (TPR[i] + TPR[i + 1]) * (FPR[i] - FPR[i + 1]) / 2.;
	}
	delete[] TPR;
	delete[] FPR;
	delete[] yTrueData;
	delete[] yPredictData;
	return auc;
}
