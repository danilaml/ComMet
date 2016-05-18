#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <algorithm>
#include <functional>

using namespace std;

typedef vector<vector<double>> matrix;

double simpsonKF(double a, double b, unsigned int m, function<double (double)> f) {
	double h = (b - a)/m;
	double xk = a;
	double res = 0;
	for (unsigned int i = 0; i < m; ++i) {
		res += f(xk) + 4 * f(xk + h/2) + f(xk + h);
		xk += h;
	}
	res *= h/6;
	return res;
}

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

void printMat(const matrix &A) {
	for (auto &row : A) {
		cout << "( ";
		for(auto x : row) {
			cout << x << " ";
		}
		cout << ")" << endl;
	}
}

void printVec(const vector<double> &v) {
	cout << "( ";
	for (auto x : v) {
		cout << x << " ";
	}
	cout << ")" << endl;
}

function<double (double)> P(int k, int n) {
	if (n == 0) {
		return [](double x){ return 1; };
	} else if (n == 1) {
		return [k](double x){ return (k + 1)*x; };
	} else {
		n -= 2;
		return [k, n](double x){
			return ((n + k + 2)*(2*n + 2*k + 3)*x*P(k, n + 1)(x) - (n + k + 2)*(n + k + 1)*P(k, n)(x))/(n + 2*k + 2)/(n + 2);
		};
	}
}

function<double (double)> W(int i) {
	if (i == 1) {
		return [](double x){ return x*x - 2*x + 1; };
	} else if (i == 2) {
		return [](double x){ return x*x*x - 3*x - 2; };
	} else {
		return [i](double x){ return (1 - x*x)*(1 - x*x) * P(2, i - 3)(x); };
	}
}

function<double (double)> dW(int i) {
	if (i == 1) {
		return [](double x){ return 2*x - 2; };
	} else if (i == 2) {
		return [](double x){ return 3*x*x - 3; };
	} else {
		int n = i - 3;
		return [n](double x){ return -2*(n + 1)*(1-x*x)*P(1, n + 1)(x); };
	}
}

function<double (double)> ddW(int i) {
	if (i == 1) {
		return [](double x){ return 2; };
	} else if (i == 2) {
		return [](double x){ return 6*x; };
	} else {
		int n = i - 3;
		return [n](double x){ return 4*(n + 1)*(n + 2)*P(0, n + 2)(x); };
	}
}

function<double (double)> LW(int i) {
	return [i](double x){
		return -(6. + x)/(7. + 3.*x)*ddW(i)(x) - (1. - x/2.)*dW(i)(x) + (1. + cos(x)/2.)*W(i)(x);
	};
}

int main()
{
	const int N = 7;
	double x0;
	cout << "x0 = ";
	cin >> x0;
	if (fabs(x0) > 1) {
		cout << "Warning: x0 not in [-1, 1]!" << endl;
	}
	auto f = [](double x){ return 1. - x/3.;};
	cout.precision(10);
	cout << fixed;
	cout << string(30, '-') << "MinSqrs" << string(30, '-') << endl;
	cout << " n |    x = -1    |    x = 0     |    x = 1     |  x = " << x0 << endl;
	cout << string(67, '-') << endl;
	for (int n = 3; n <= N; ++n) {
		matrix Ab(n, vector<double>(n + 1));
		for (int i = 1; i <= n; ++i) {
			for (int j = 1; j <= n; ++j) {
				Ab[i - 1][j - 1] = simpsonKF(-1, 1, 20, [i, j](double x){ return LW(j)(x) * LW(i)(x); });
			}
			Ab[i - 1][n] = simpsonKF(-1, 1, 20, [i, &f](double x){ return f(x) * LW(i)(x); });
		}
		auto c = columnGauss(Ab);
		auto y = [n, &c](double x){
			double res = 0;
			for (int i = 0; i < n; ++i) {
				res += c[i] * W(i + 1)(x);
			}
			return res;
		};
		cout << " " << n << " | " << y(-1) << " | " << y(0) << " | " << y(1) << " | " << y(x0) << endl;
		cout << string(67, '-') << endl;
		if (n == N) {
			cout << "Ab: " << endl;
			printMat(Ab);
			cout << "Coefficients: ";
			printVec(c);
		}
	}
	cout << string(30, '-') << "Collocs" << string(30, '-') << endl;
	cout << " n |    x = -1    |    x = 0     |    x = 1     |  x = " << x0 << endl;
	cout << string(67, '-') << endl;
	for (int n = 3; n <= N; ++n) {
		matrix Ab(n, vector<double>(n + 1));
		for (int i = 1; i <= n; ++i) {
			double ti = cos((2*i - 1)*3.1416/2/n);
			for (int j = 1; j <= n; ++j) {
				Ab[i - 1][j - 1] = LW(j)(ti);
			}
			Ab[i - 1][n] = f(ti);
		}
		auto c = columnGauss(Ab);
		auto y = [n, &c](double x){
			double res = 0;
			for (int i = 0; i < n; ++i) {
				res += c[i] * W(i + 1)(x);
			}
			return res;
		};
		cout << " " << n << " | " << y(-1) << " | " << y(0) << " | " << y(1) << " | " << y(x0) << endl;
		cout << string(67, '-') << endl;
		if (n == N) {
			cout << "Ab: " << endl;
			printMat(Ab);
			cout << "Coefficients: ";
			printVec(c);
		}
	}
	return 0;
}

