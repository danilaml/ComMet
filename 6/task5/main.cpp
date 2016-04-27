#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <algorithm>
#include <functional>

using namespace std;

typedef vector<vector<double>> matrix;

vector<double> columnGauss(matrix Ab) {
	for(unsigned i = 0; i < Ab.size(); ++i) {
		unsigned maxind = i;
		for (unsigned k = i + 1; k < Ab.size(); ++k) {
			if (fabs(Ab[k][i]) > fabs(Ab[maxind][i]))
				maxind = k;
		}
		if (i != maxind) swap(Ab[i], Ab[maxind]);

		double head = Ab[i][i];
		for (unsigned j = i + 1; j < Ab[0].size(); ++j) {
			Ab[i][j] /= head;
		}
		for (unsigned j = i + 1; j < Ab.size(); ++j) {
			auto tmp = Ab[j][i];
			for (unsigned k = i + 1; k < Ab[0].size(); ++k) {
				Ab[j][k] -= Ab[i][k]*tmp;
			}
		}
	}

	vector<double> x(Ab.size());
	x[Ab.size() - 1] = Ab.back()[Ab.size()];
	for (int i = Ab.size() - 2; i >= 0; --i) {
		double xk = Ab[i][Ab.size()];
		for (int j = Ab.size() - 1; j > i; --j) {
			xk -= Ab[i][j]*x[j];
		}
		x[i] = xk;
	}
	return x;
}

void printVec(const vector<double> &v) {
	cout << "( ";
	for (auto x : v) {
		cout << x << " ";
	}
	cout << ")" << endl;
}

void printMat(const matrix &A) {
	for (auto &row : A) {
		cout << "( ";
		for(auto x : row) {
			cout << x << " ";
		}
		cout << ")" << endl;
	}
}

double vecNorm(const vector<double> &v) {
	return max(*max_element(v.begin(), v.end()), -(*min_element(v.begin(), v.end())));
}

vector<double> vecSub(const vector<double> &l, const vector<double> &r) {
	vector<double> res(l.size());
	for (unsigned i = 0; i < l.size(); ++i) {
		res[i] = l[i] - r[i];
	}
	return res;
}

void recalcC(double a, double b, vector<double> &c) {
	double k = (b - a)/2;
	for (double &x : c) {
		x *= k;
	}
}

void recalcX(double a, double b, vector<double> &xs) {
	double k1 = (b - a)/2;
	double k2 = (b + a)/2;
	for (double &x : xs) {
		x = x*k1 + k2;
	}
}

vector<double> mechKvadr(const vector<double> &c
						 , const vector<double> &xs
						 , const function<double (double, double)> K
						 , const function<double (double)> f) {
	double N = c.size();
	matrix Gf(N, vector<double>(N + 1));
	for (int i = 0; i < N; ++i) {
		for (int j = 0; j < N; ++j) {
			Gf[i][j] = c[j]*K(xs[i], xs[j]);
			Gf[i][j] += i == j;
		}
		Gf[i][N] = f(xs[i]);
	}
	cout << "Gf:" << endl;
	printMat(Gf);
	cout << endl;
	return columnGauss(Gf);
}

int main()
{
	double a = 0;
	double b = 1;
	auto K = [](double x, double y){ return log(2 + x + y)/3; };
	auto f = [](double x){ return 2 - x - x*x; };
	vector<double> x3 = {-0.7745966692, 0, 0.7745966692};
	vector<double> x5 = {-0.9061798459, -0.5384693101, 0, 0.5384693101, 0.9061798459};
	vector<double> c3 = {5./9., 8./9., 5./9.};
	vector<double> c5 = {0.2369268851, 0.4786286705, 0.5688888889, 0.4786286705, 0.2369268851};
	recalcX(a, b, x3);
	recalcX(a, b, x5);
	recalcC(a, b, c3);
	recalcC(a, b, c5);
	cout << "x3: ";
	printVec(x3);
	cout << "x5: ";
	printVec(x5);
	cout << endl;
	cout << "c3: ";
	printVec(c3);
	cout << "c5: ";
	printVec(c5);
	cout << endl;


	auto res3 = mechKvadr(c3, x3, K, f);
	auto res5 = mechKvadr(c5, x5, K, f);

	auto U3 = [&](double x){
		double res = f(x);
		for (unsigned i = 0; i < x3.size(); ++i) {
			res -= c3[i]*K(x, x3[i])*res3[i];
		}
		return res;
	};

	auto U5 = [&](double x){
		double res = f(x);
		for (unsigned i = 0; i < x5.size(); ++i) {
			res -= c5[i]*K(x, x5[i])*res5[i];
		}
		return res;
	};

	cout.precision(10);
	cout << fixed;
	cout << "u1, u2, u3: ";
	printVec(res3);
	cout << "u1-u5: ";
	printVec(res5);

	double N = 10;
	double maxR = fabs(U3(a) - U5(a));
	for (double z = a; z <= b; z += (b - a)/N) {
		maxR = max(maxR, fabs(U3(z) - U5(z)));
		cout << "U3(" << z << ") = " << U3(z) << endl;
		cout << "U5(" << z << ") = " << U5(z) << endl << endl;
	}

	cout << "||(U5(zk) - U3(zk))|| = " << maxR << endl;

	return 0;
}

