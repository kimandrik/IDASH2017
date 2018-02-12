#include "GD.h"

#include <NTL/ZZ.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <StringUtils.h>

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

double** GD::RandomxyDataLearn(double** xyData, long sampleDimLearn, long sampleDim, long factorDim) {
	double** res = new double*[sampleDimLearn];

	bool* notTaken = new bool[sampleDim];
	fill_n(notTaken, sampleDim, false);
	long r;
	for (long j = 0; j < sampleDim - sampleDimLearn; ++j) {
		do {
			r = RandomBnd(sampleDim);
		} while(notTaken[r]);
		notTaken[r] = true;
	}
	long idx = 0;
	for (long j = 0; j < sampleDim; ++j) {
		if(!notTaken[j]) {
			res[idx] = new double[factorDim];
			for (long i = 0; i < factorDim; ++i) {
				res[idx][i] = xyData[j][i];
			}
			idx++;
		}
	}
	delete[] notTaken;
	return res;
}

void GD::normalizexyData(double** xyData, long factorDim, long sampleDim) {
	for (long i = 0; i < factorDim; ++i) {
		double m = 0;
		for (long j = 0; j < sampleDim; ++i) {
			m = max(m, abs(xyData[j][i]));
		}
		if(m < 1e-10) continue;
		for (long j = 0; j < sampleDim; ++j) {
			xyData[j][i] /= m;
		}
	}
}

void GD::normalizexyData2(double** xyDataLearn, double** xyDataTest, long factorDim, long sampleDimLearn, long sampleDimTest) {
	long i, j;
	double m;
	for (i = 0; i < factorDim; ++i) {
		m = 0.0;
		for (j = 0; j < sampleDimLearn; ++j) {
			m = max(m, abs(xyDataLearn[j][i]));
		}
		for (j = 0; j < sampleDimTest; ++j) {
			m = max(m, abs(xyDataTest[j][i]));
		}
		if(m < 1e-10) continue;
		for (j = 0; j < sampleDimLearn; ++j) {
			xyDataLearn[j][i] /= m;
		}
		for (j = 0; j < sampleDimTest; ++j) {
			xyDataTest[j][i] /= m;
		}
	}
}

double* GD::plainIP(double** a, double* b, long factorDim, long sampleDim) {
	double* res = new double[sampleDim]();
	for (long j = 0; j < sampleDim; ++j) {
		for(long i = 0; i < factorDim; ++i) {
			res[j] += a[j][i] * b[i];
		}
	}
	return res;
}

double* GD::plainSigmoid(long approxDeg, double** xyData, double* ip, long factorDim, long sampleDim, double gamma) {
	double* grad = new double[factorDim]();
	if(approxDeg == 3) {
		for (long i = 0; i < factorDim; ++i) {
			for (long j = 0; j < sampleDim; ++j) {
				grad[i] += (degree3[0] + ip[j] * degree3[1] + pow(ip[j], 3) * degree3[2]) * xyData[j][i];
			}
			grad[i] *= gamma;
		}
	} else if(approxDeg == 5) {
		for (long i = 0; i < factorDim; ++i) {
			for (long j = 0; j < sampleDim; ++j) {
				grad[i] += (degree5[0] + ip[j] * degree5[1] + pow(ip[j], 3) * degree5[2] + pow(ip[j], 5) * degree5[3]) * xyData[j][i];
			}
			grad[i] *= gamma;
		}
	} else {
		for (long i = 0; i < factorDim; ++i) {
			for (long j = 0; j < sampleDim; ++j) {
				grad[i] += (degree7[0] + ip[j] * degree7[1] + pow(ip[j], 3) * degree7[2] + pow(ip[j], 5) * degree7[3] + pow(ip[j], 7) * degree7[4]) * xyData[j][i];
			}
			grad[i] *= gamma;
		}
	}
	return grad;
}

void GD::plainLGDstep(double* wData, double* grad, long factorDim) {
	for (long i = 0; i < factorDim; ++i) {
		wData[i] -= grad[i];
	}
}

void GD::plainMLGDstep(double* wData, double* vData, double* grad, long factorDim, double eta) {
	for (long i = 0; i < factorDim; ++i) {
		vData[i] = eta * vData[i] + grad[i];
		wData[i] -= vData[i];
	}
}

void GD::plainNLGDstep(double* wData, double* vData, double* grad, long factorDim, double eta) {
	for (long i = 0; i < factorDim; ++i) {
		double tmpw = vData[i] - grad[i];
		vData[i] = (1.0 - eta) * tmpw + eta * wData[i];
		wData[i] = tmpw;
	}
}

void GD::plainLGDL2step(double* wData, double* grad, long factorDim, double lambda) {
	//TODO: implement method
}

void GD::plainMLGDL2step(double* wData, double* vData, double* grad, long factorDim, double eta, double lambda) {
	//TODO: implement method
}

void GD::plainNLGDL2step(double* wData, double* vData, double* grad, long factorDim, double eta, double lambda) {
	//TODO: implement method
}

void GD::plainLGDiteration(long approxDeg, double** xyData, double* wData, long factorDim, long sampleDim, double gamma) {
	double* ip = plainIP(xyData, wData, factorDim, sampleDim);
	double* grad = plainSigmoid(approxDeg, xyData, ip, factorDim, sampleDim, gamma);
	plainLGDstep(wData, grad, factorDim);
	delete[] ip;
	delete[] grad;
}

void GD::plainMLGDiteration(long approxDeg, double** xyData, double* wData, double* vData, long factorDim, long sampleDim, double gamma, double eta) {
	double* ip = plainIP(xyData, wData, factorDim, sampleDim);
	double* grad = plainSigmoid(approxDeg, xyData, ip, factorDim, sampleDim, gamma);
	plainMLGDstep(wData, vData, grad, factorDim, eta);
	delete[] ip;
	delete[] grad;
}


void GD::plainNLGDiteration(long approxDeg, double** xyData, double* wData, double* vData, long factorDim, long sampleDim, double gamma, double eta) {
	double* ip = plainIP(xyData, vData, factorDim, sampleDim);
	double* grad = plainSigmoid(approxDeg, xyData, ip, factorDim, sampleDim, gamma);
	plainNLGDstep(wData, vData, grad, factorDim, eta);

	delete[] ip;
	delete[] grad;
}

void GD::plainLGDL2iteration(long approxDeg, double** xyData, double* wData, long factorDim, long sampleDim, double gamma, double lambda) {
	//TODO: implement method
}

void GD::plainMLGDL2iteration(long approxDeg, double** xyData, double* wData, double* vData, long factorDim, long sampleDim, double gamma, double eta, double lambda) {
	//TODO: implement method
}

void GD::plainNLGDL2iteration(long approxDeg, double** xyData, double* wData, double* vData, long factorDim, long sampleDim, double gamma, double eta, double lambda) {
	//TODO: implement method
}

//-----------------------------------------

double GD::trueIP(double* a, double* b, long size) {
	double res = 0.0;
	for(long i = 0; i < size; ++i) {
		res += a[i] * b[i];
	}
	return res;
}

void GD::trueLGDiteration(double** xyData, double* wData, long factorDim, long sampleDim, double gamma) {
	double* grad = new double[factorDim]();

	for(long j = 0; j < sampleDim; ++j) {
		double ip = trueIP(wData, xyData[j], factorDim);
		double tmp = (ip > 15.0) ? 0 : (ip < -15.0) ? -1.0 : - 1. / (1. + exp(ip));
		for(int i = 0; i < factorDim; ++i) {
			grad[i] += tmp * (double) xyData[j][i];
		}
	}

	for (int i = 0; i < factorDim; ++i) {
		wData[i] -= gamma * grad[i];
	}
	delete[] grad;
}

void GD::trueMLGDiteration(double** xyData, double* wData, double* vData, long factorDim, long sampleDim, double gamma, double eta) {
	double* grad = new double[factorDim]();

	for(long j = 0; j < sampleDim; ++j) {
		double ip = trueIP(wData, xyData[j], factorDim);
		double tmp = (ip > 15.0) ? 0 : (ip < -15.0) ? -1.0 : - 1. / (1. + exp(ip));
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

void GD::trueNLGDiteration(double** xyData, double* wData, double* vData, long factorDim, long sampleDim, double gamma, double eta) {
	double* grad = new double[factorDim]();

	for(long j = 0; j < sampleDim; ++j) {
		double ip = trueIP(vData, xyData[j], factorDim);
		double tmp = (ip > 15.0) ? 0 : (ip < -15.0) ? -1.0 : - 1. / (1. + exp(ip));
		for(long i = 0; i < factorDim; ++i) {
			grad[i] += tmp * (double) xyData[j][i];
		}
	}
	for (long i = 0; i < factorDim; ++i) {
		double tmpw = vData[i] - gamma * grad[i];
		vData[i] = (1.0 - eta) * tmpw + eta * wData[i];
		wData[i] = tmpw;
	}
	delete[] grad;
}

void GD::trueLGDiteration(double** xyData, double* wData, long factorDim, long sampleDim, double gamma, double lambda) {
	//TODO: implement method
}

void GD::trueMLGDiteration(double** xyData, double* wData, double* vData, long factorDim, long sampleDim, double gamma, double eta, double lambda) {
	//TODO: implement method
}

void GD::trueNLGDiteration(double** xyData, double* wData, double* vData, long factorDim, long sampleDim, double gamma, double eta, double lambda) {
	//TODO: implement method
}

void GD::check(double** xyData, double* wData, long factorDim, long sampleDim) {
	cout << "w:";
	for (long i = 0; i < factorDim; ++i) {
		cout << wData[i] << ",";
	}
	cout << endl;
	long num = 0;
	for(long j = 0; j < sampleDim; ++j){
		if(trueIP(wData, xyData[j], factorDim) > 0) num++;
	}
	cout << "Correctness: " << num << "/" << sampleDim << endl;

}

double GD::calculateAUC(double** xyData, double* wData, long factorDim, long sampleDim) {
    long n_fail_y1 = 0;
    long n_fail_y0 = 0;

    vector<double> xtheta_y1;
    vector<double> xtheta_y0;

    for(int i = 0; i < sampleDim; ++i){
        if(xyData[i][0] > 0){
            if(GD::trueIP(xyData[i], wData, factorDim) < 0) n_fail_y1++;
            xtheta_y1.push_back(xyData[i][0] * GD::trueIP(xyData[i] + 1, wData + 1, factorDim - 1));
        } else{
            if(GD::trueIP(xyData[i], wData, factorDim) < 0) n_fail_y0++;
            xtheta_y0.push_back(xyData[i][0] * GD::trueIP(xyData[i] + 1, wData + 1, factorDim - 1));
        }
    }

    double correctness = 100.0 - (100.0 * (n_fail_y0 + n_fail_y1) / sampleDim);
    cout << "Failure rate: (y = 1) " << n_fail_y1 << "/" << xtheta_y1.size() << " + (y = 0) " << n_fail_y0 << "/" ;
    cout << xtheta_y0.size() << " = " <<  (100.0 * (n_fail_y0 + n_fail_y1) / sampleDim) << " %." << endl;
    cout << "Correctness: " << correctness  << " %." << endl;

    if(xtheta_y0.size() == 0 || xtheta_y1.size() ==0){
        cout << "n_test_yi = 0 : cannot compute AUC" << endl;
        return 0.0;
    } else{
        double auc = 0.0;
        for(long i = 0; i < xtheta_y1.size(); ++i){
            for(long j = 0; j < xtheta_y0.size(); ++j){
                if(xtheta_y0[j] <= xtheta_y1[i]) auc++;
            }
        }
        auc /= xtheta_y1.size() * xtheta_y0.size();
        return auc;
        cout << "AUC: " << auc << endl;
    }
}


double GD::calculateMSE(double* wData1, double* wData2, long factorDim) {
    double res= 0.0;

    for(long i = 0; i < factorDim; ++i) {
        res += pow(wData1[i] - wData2[i], 2.0);
    }
    res /= factorDim;

    return res;
}


double GD::calculateNMSE(double* wData1, double* wData2, long factorDim) {
    double res= 0.0;

    for(long i = 0; i < factorDim; ++i) {
        res += pow(wData1[i], 2.0);
    }
    res /= factorDim;

    double mse = GD::calculateMSE(wData1, wData2, factorDim);
    res = mse / res;

    return res;

}
