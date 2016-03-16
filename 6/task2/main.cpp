#include <iostream>
#include <vector>
#include <valarray>
#include <iomanip>

using namespace std;

typedef vector<valarray<float>> matrix;

matrix unify(const matrix &A, const valarray<float> &b) {
	matrix result(A.size());
	for (unsigned i = 0; i < A.size(); ++i) {
		valarray<float> tmp(A.size() + 1);
		for (unsigned j = 0; j < A.size(); ++j) {
			tmp[j] = A[i][j];
		}
		tmp[A.size()] = b[i];
		result[i] = tmp;
	}
	return result;
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

valarray<float> simpleGauss(const matrix &Aa, const valarray<float> &b, bool printDeterm = false) {
	matrix A = unify(Aa, b);
	float determ = 1.0;
	for(unsigned i = 0; i < A.size(); ++i) {
		float head = A[i][i];
		A[i] /= head;
		determ *= head;
		for (unsigned j = i + 1; j < A.size(); ++j) {
			A[j] -= valarray<float>(A[i]*A[j][i]);
		}
		printMat(A);
		cout << endl;
	}

	if (printDeterm) {
		cout << "determinant: " << abs(determ) << endl;
	}

	valarray<float> x(A.size());
	x[A.size() - 1] = A.back()[A.size()];
	for (int i = A.size() - 2; i >= 0; --i) {
		float xk = A[i][A.size()];
		for (int j = A.size() - 1; j > i; --j) {
			xk -= A[i][j]*x[j];
		}
		x[i] = xk/A[i][i];
	}
	return x;
}

valarray<float> simpleGaussB(const matrix &Aa, const valarray<float> &b, bool printDeterm = false) {
	matrix A = unify(Aa, b);
	float determ = 1.0;
	for(unsigned i = 0; i < A.size(); ++i) {
		float head = A[i][i];
		for (unsigned j = i + 1; j < A[0].size(); ++j) {
			A[i][j] /= head;
		}
		determ *= head;
		for (unsigned j = i + 1; j < A.size(); ++j) {
			auto tmp = A[j][i];
			for (unsigned k = i + 1; k < A[0].size(); ++k) {
				A[j][k] -= A[i][k]*tmp;
			}
			//A[j] -= valarray<float>(A[i]*A[j][i]);
		}
		if (printDeterm) {
			printMat(A);
			cout << endl;
		}
	}

	if (printDeterm) {
		cout << "determinant: " << abs(determ) << endl;
	}

	valarray<float> x(A.size());
	x[A.size() - 1] = A.back()[A.size()];
	for (int i = A.size() - 2; i >= 0; --i) {
		float xk = A[i][A.size()];
		for (int j = A.size() - 1; j > i; --j) {
			xk -= A[i][j]*x[j];
		}
		x[i] = xk;//A[i][i];
	}
	return x;
}

valarray<float> columnGauss(const matrix &Aa, const valarray<float> &b, bool printDeterm = false) {
	matrix A = unify(Aa, b);
	float determ = 1.0;
	for(unsigned i = 0; i < A.size(); ++i) {
		unsigned maxind = i;
		for (unsigned k = i + 1; k < A.size(); ++k) {
			if (abs(A[k][i]) > abs(A[maxind][i]))
				maxind = k;
		}
		if (i != maxind) swap(A[i], A[maxind]);

		float head = A[i][i];
		for (unsigned j = i + 1; j < A[0].size(); ++j) {
			A[i][j] /= head;
		}
		determ *= head;
		for (unsigned j = i + 1; j < A.size(); ++j) {
			auto tmp = A[j][i];
			for (unsigned k = i + 1; k < A[0].size(); ++k) {
				A[j][k] -= A[i][k]*tmp;
			}
			//A[j] -= valarray<float>(A[i]*A[j][i]);
		}
		if (printDeterm) {
			printMat(A);
			cout << endl;
		}
	}

	if (printDeterm) {
		cout << "determinant: " << abs(determ) << endl;
	}

	valarray<float> x(A.size());
	x[A.size() - 1] = A.back()[A.size()];
	for (int i = A.size() - 2; i >= 0; --i) {
		float xk = A[i][A.size()];
		for (int j = A.size() - 1; j > i; --j) {
			xk -= A[i][j]*x[j];
		}
		x[i] = xk;//A[i][i];
	}
	return x;
}

matrix findInverse(const matrix &A) {
	matrix result(A.size(), valarray<float>(A.size()));
	for (unsigned i = 0; i < A.size(); ++i) {
		valarray<float> col(0.0, A.size());
		col[i] = 1;
		auto tmp = columnGauss(A, col);
		for (unsigned j = 0; j < A.size(); ++j) {
			result[j][i] = tmp[j];
		}
	}
	return result;
}

valarray<float> matrixVector(const matrix &A, const valarray<float> &v) {
	valarray<float> result(A.size());
	for (unsigned i = 0; i < A.size(); ++i) {
		result[i] = valarray<float>(A[i]*v).sum();
	}
	return result;
}

int main()
{
	matrix A =  {{6.7528E-06, -7.5944E-03, 4.30584},
				{6.1528E-03, -0.75944, 1.53584},
				{0.89528, 0.84056, 0.98584}};/*{{2,1,0,-3},{3,-1,0,1},{1,4,-2,-5}};*/
	valarray<float> b = {3.80430, 1.64243, 1.96424};

	valarray<float> res;
	valarray<float> tmp;

	cout << "A:" << endl;
	cout.precision(10);
	int ind = 0;
	for (auto &row : A) {
		cout << "( ";
		for(auto x : row) {
			cout << x << " ";
		}
		cout << "| " << b[ind] << " )" << endl;
		ind++;
	}

//	cout << endl << "simple:" << endl;
//	res = simpleGauss(A, b, true);
//	cout << "x:" << endl;
//	for (auto x : res) {
//		cout << x << endl;
//	}
	cout << endl << "simple:" << endl;
	res = simpleGaussB(A, b, true);
	cout << "x:" << endl;
	ind = 0;
	for (auto x : res) {
		cout << "x[" << ind << "] = " << x << endl;
		ind++;
	}
	tmp = matrixVector(A, res) - b;
	cout << "|| r || = " << max(tmp.max(), - tmp.min()) << endl;
	cout << endl;

	res = columnGauss(A, b, true);
	cout << "x:" << endl;
	ind = 0;
	for (auto x : res) {
		cout << "x[" << ind << "] = " << x << endl;
		ind++;
	}
	tmp = matrixVector(A, res) - b;
	cout << "|| r || = " << max(tmp.max(), - tmp.min()) << endl;
	cout << endl;

	cout << "inverse:" << endl;
	auto invA = findInverse(A);
	printMat(invA);
//	for (auto &row : invA) {
//		cout << "( ";
//		for(auto x : row) {
//			cout << x << " ";
//		}
//		cout << ")" << endl;
//	}
	cout << endl << "x:" << endl;
	res = matrixVector(invA, b);
	ind = 0;
	for (auto x : res) {
		cout << "x[" << ind << "] = " << x << endl;
		ind++;
	}
	tmp = matrixVector(A, res) - b;
	cout << "|| r || = " << max(tmp.max(), - tmp.min()) << endl;

	return 0;
}

