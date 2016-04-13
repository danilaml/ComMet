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

vector<double> scMultVec(double s, const vector<double> &v) {
	vector<double> res(v.size());
	for (unsigned i = 0; i < v.size(); ++i) {
		res[i] = s * v[i];
	}
	return res;
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

double vecNorm(const vector<double> &v) {
	return max(*max_element(v.begin(), v.end()), -(*min_element(v.begin(), v.end())));
}

void normVec(vector<double> &v) {
	double norm = vecNorm(v);
	for (double &x : v) {
		x /= norm;
	}
}

vector<double> vecSub(const vector<double> &l, const vector<double> &r) {
	vector<double> res(l.size());
	for (unsigned i = 0; i < l.size(); ++i) {
		res[i] = l[i] - r[i];
	}
	return res;
}

void jacob(matrix A, double eps) {
	matrix U = {{1, 0, 0},
				{0, 1, 0},
				{0, 0, 1}};
	matrix Acop = A;
	auto sign = [](double x){ return x > 0 ? 1 : (x < 0 ? -1 : 0); };

	for(int step = 0;; ++step) {
		double maxi = 0;
		double maxj = 1;
		double maxa = A[0][1];
		for (unsigned i = 0; i < A.size() - 1; ++i) {
			for (unsigned j = i + 1; j < A.size(); ++j) {
				if (fabs(A[i][j]) > fabs(maxa)) {
					maxi = i;
					maxj = j;
					maxa = A[i][j];
				}
			}
		}
		cout << "i: " << maxi << ", j: " << maxj << ", max: " << maxa << endl;
		if (fabs(maxa) < eps)
			break;

		double aii = A[maxi][maxi];
		double ajj = A[maxj][maxj];

		double delt = fabs(aii - ajj);
		double sg = sign(maxa*(aii - ajj));
		double d = sqrt(delt*delt + 4*maxa*maxa);
		double c = sqrt(0.5*(1 + delt/d));
		double s = (delt == 0) ? c : sg*sqrt(0.5*(1 - delt/d));

		cout << "c: " << c << ", s: " << s << endl;

		for (unsigned m = 0; m < A.size(); ++m) {
			if (m != maxi && m != maxj) {
				double oldmi = A[m][maxi];
				double oldmj = A[m][maxj];
				A[m][maxi] = c*oldmi + s*oldmj;
				A[maxi][m] = A[m][maxi];
				A[m][maxj] = -s*oldmi + c*oldmj;
				A[maxj][m] = A[m][maxj];
			}
		}
		A[maxi][maxi] = c*c*aii + 2*c*s*maxa + s*s*ajj;
		A[maxj][maxj] = s*s*aii - 2*c*s*maxa + c*c*ajj;
		A[maxi][maxj] = 0;
		A[maxj][maxi] = 0;
		cout << "A" << step << ":" << endl;
		printMat(A);
		cout << endl;
		for (unsigned i = 0; i < U.size(); ++i) {
			double oldi = U[i][maxi];
			double oldj = U[i][maxj];
			U[i][maxi] = c*oldi + s*oldj;
			U[i][maxj] = -s*oldi + c*oldj;
		}
		//printMat(U);
		//cout << endl;
	}

	cout << "Final (diagonalized) A:" << endl;
	printMat(A);
	cout << "U:" << endl;
	printMat(U);
	vector<double> x1, x2, x3;
	x1 = {U[0][0], U[1][0], U[2][0]};
	x2 = {U[0][1], U[1][1], U[2][1]};
	x3 = {U[0][2], U[1][2], U[2][2]};
	normVec(x1);
	normVec(x2);
	normVec(x3);
	cout << "Eigenvectors normalized:" << endl;
	cout << "x1: ";
	printVec(x1);
	cout << "x2: ";
	printVec(x2);
	cout << "x3: ";
	printVec(x3);
	auto nev = [&Acop](vector<double> &x, double lamb){
		return vecNorm(vecSub(matMultVec(Acop, x), scMultVec(lamb, x)));
	};
	cout << "||A*x1 - l1*x1|| : " << nev(x1, A[0][0]) << endl;
	cout << "||A*x2 - l2*x2|| : " << nev(x2, A[1][1]) << endl;
	cout << "||A*x3 - l3*x3|| : " << nev(x3, A[2][2]) << endl;
}

int main()
{
	matrix A2 = {{5, 10, 5},
				 {10, 11, -1},
				 {5, -1, 2}};
	double eps = 1e-5;

	cout.precision(10);
	cout << fixed;

	cout << "Initial A:" << endl;
	printMat(A2);
	cout << "eps: " << eps << endl;
	cout << endl << "Jacobi method:" << endl;

	jacob(A2, eps);

	return 0;
}

