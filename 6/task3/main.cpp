#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <algorithm>

using namespace std;

typedef vector<vector<double>> matrix;

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

vector<double> columnGauss(matrix Ab, bool printDeterm = false) {
	double determ = 1.0;
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
		determ *= head;
		for (unsigned j = i + 1; j < Ab.size(); ++j) {
			auto tmp = Ab[j][i];
			for (unsigned k = i + 1; k < Ab[0].size(); ++k) {
				Ab[j][k] -= Ab[i][k]*tmp;
			}
		}
		if (printDeterm) {
			printMat(Ab);
			cout << endl;
		}
	}

	if (printDeterm) {
		cout << "determinant: " << fabs(determ) << endl;
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

vector<double> matMultVec(const matrix &A, const vector<double> &v) {
	vector<double> result(A.size());
	for (unsigned i = 0; i < A.size(); ++i) {
		for (unsigned j = 0; j < A.size(); ++j) {

			result[i] += A[i][j]*v[j];
		}
	}
	return result;
}

vector<double> vecSum(const vector<double> &l, const vector<double> &r) {
	vector<double> res(l.size());
	for (unsigned i = 0; i < l.size(); ++i) {
		res[i] = l[i] + r[i];
	}
	return res;
}

vector<double> scMultVec(double s, const vector<double> &v) {
	vector<double> res(v.size());
	for (unsigned i = 0; i < v.size(); ++i) {
		res[i] = s * v[i];
	}
	return res;
}

vector<double> vecSub(const vector<double> &l, const vector<double> &r) {
	vector<double> res(l.size());
	for (unsigned i = 0; i < l.size(); ++i) {
		res[i] = l[i] - r[i];
	}
	return res;
}

double findm(const matrix &Ab) {
	double res = +INFINITY;
	for (unsigned i = 0; i < Ab.size(); ++i) {
		double sum = Ab[i][i];
		for (unsigned j = 0; j < Ab.size(); ++j) {
			sum -= i != j ? fabs(Ab[i][j]) : 0;
		}
		res = min(res, sum);
	}
	return res;
}

double findM(const matrix &Ab) {
	double res = -INFINITY;
	for (unsigned i = 0; i < Ab.size(); ++i) {
		double sum = Ab[i][i];
		for (unsigned j = 0; j < Ab.size(); ++j) {
			sum += i != j ? fabs(Ab[i][j]) : 0;
		}
		res = max(res, sum);
	}
	return res;
}

double matrixNorm(const matrix &B) {
	double res = 0;
	for (unsigned i = 0; i < B[0].size(); ++i) {
		res += fabs(B[0][i]);
	}
	for (unsigned i = 1; i < B.size(); ++i) {
		double tmp = 0;
		for (unsigned j = 0; j < B[0].size(); ++j) {
			tmp += fabs(B[i][j]);
		}
		res = max(res, tmp);
	}
	return res;
}

double vecNorm(const vector<double> &v) {
	return max(*max_element(v.begin(), v.end()), -(*min_element(v.begin(), v.end())));
}

vector<int> genTheta(unsigned int kpow) {
	vector<int> theta = {1};
	while (theta.size() < kpow) {
		vector<int> tmp(theta.size() * 2);
		for (unsigned i = 0; i < theta.size(); ++i) {
			tmp[2 * i] = theta[i];
			tmp[2 * i + 1] = 4 * theta.size() - theta[i];
		}
		theta = tmp;
	}
	return theta;
}

int main()
{
	matrix Ab =  {{2.1790, 0.38540, 0.21058, 2.14603},
				  {0.38540, 3.12918, 0.59380, 6.36415},
				  {0.21058, 0.59380, 4.81900, 6.06336}};
	matrix A =  {{2.1790, 0.38540, 0.21058},
				 {0.38540, 3.12918, 0.59380},
				 {0.21058, 0.59380, 4.81900}};
	vector<double> b = {2.14603, 6.36415, 6.06336};
	cout << "Ab:" << endl;
	cout.precision(10);
	cout << fixed;
	for (auto &row : Ab) {
		cout << "( ";
		for(auto x : row) {
			cout << x << " ";
		}
		cout << ")" << endl;
	}
	auto res = columnGauss(Ab/*, true*/);
	cout << "Gauss x:" << endl;
	for (unsigned i = 0; i < res.size(); ++i) {
		cout << "x[" << i << "] = " << res[i] << endl;
	}
	cout << endl;

	double m = findm(Ab);
	double M = findM(Ab);
	cout << "m: " << m << ", M: " << M << endl;
	double alpha = 2/(m + M);
	cout << "alpha: " << alpha << endl << endl;

	matrix B(Ab.size(), vector<double>(Ab.size()));
	for (unsigned i = 0; i < Ab.size(); ++i) {
		for (unsigned j = 0; j < Ab.size(); ++j) {
			B[i][j] = i == j ? 1 - alpha * Ab[i][j] : - alpha * Ab[i][j];
		}
	}

	cout << "Ba:" << endl;
	printMat(B);
	double normB = matrixNorm(B);
	cout << "||Ba||: " << normB << endl;

	vector<double> c(Ab.size());
	for (unsigned i = 0; i < Ab.size(); ++i) {
		c[i] = alpha * Ab[i][Ab.size()];
	}

	cout << "c: ";
	printVec(c);
	cout << endl;

	double eps = 1e-6;

	vector<double> oldx = {1., 1., 1.};
	auto newx = vecSum(matMultVec(B, oldx), c);
	int k_apri = static_cast<int>(
				ceil(log(eps/vecNorm(vecSub(oldx, newx))*(1 - normB))/log(normB))
				);
	cout << "k_apri: " << k_apri << endl << endl;
	cout << "x0: ";
	printVec(oldx);
	cout << "||r0||: " << vecNorm(vecSub(matMultVec(A, oldx), b)) << endl << endl;
	cout << "||x0 - x*||: " << vecNorm(vecSub(res, oldx)) << endl;
	cout << "x1: ";
	printVec(newx);
	oldx.swap(newx);
	double phi = normB/(1 - normB);
	double apri = phi * vecNorm(vecSub(oldx, newx));
	double apost = phi * vecNorm(vecSub(oldx, newx));
	cout << "||r1||: " << vecNorm(vecSub(matMultVec(A, oldx), b)) << endl;
	cout << "||x1 - x*||: " << vecNorm(vecSub(res, oldx)) << endl;
	cout << "apri_1: " << apri << endl;
	cout << "apost_1: " << apost << endl << endl;
	for (int i = 2;; ++i) {
		newx = vecSum(matMultVec(B, oldx), c);
		if (vecNorm(vecSub(oldx, newx)) < eps)
			break;
		else
			oldx.swap(newx);
		cout << "x" << i << ": ";
		printVec(oldx);
		cout << "||r" << i << "||: " << vecNorm(vecSub(matMultVec(A, oldx), b)) << endl;
		cout << "||x" << i << " - x*||: " << vecNorm(vecSub(res, oldx)) << endl;
		apri *= normB;
		apost = phi * vecNorm(vecSub(oldx, newx));
		cout << "apri_" << i << ": " << apri << endl;
		cout << "apost_" << i << ": " << apost << endl << endl;
	};

	cout << "---------------MPI---------------" << endl << endl;

	unsigned int k = static_cast<int>(pow(2, ceil(log2(sqrt(M/m)*log(2/eps)/2))));

	cout << "k(eps): " << k << endl;

	auto theta = genTheta(k);
	cout << "theta: ( ";
	for (auto x : theta) {
		cout << x << " ";
	}
	cout << " )" << endl << endl;

	vector<double> x = {1, 1, 1};
	cout << "x0: ";
	printVec(x);
	for (unsigned i = 1; i <= k; ++i) {
		double tau = 2/(M + m - (M - m)*cos(theta[i - 1]*3.14159265359/2/k));
		x = vecSum(x, scMultVec(tau, vecSub(b, matMultVec(A, x))));
		cout << "x" << i << ": ";
		printVec(x);
		cout << "||r" << i << "||: " << vecNorm(vecSub(matMultVec(A, x), b)) << endl;
		cout << "||x" << i << " - x*||: " << vecNorm(vecSub(res, x)) << endl << endl;
	};

	return 0;
}

