#include <iostream>
#include <functional>
#include <vector>
#include <cmath>
#include <iomanip>

using namespace std;

vector<pair<double, double>> createTable(double a, double b, unsigned int m, function<double (double)> f) {
	vector<pair<double, double>> res;
	double h = (b - a) / m;
	double x = a;
	for (unsigned int i = 0; i <= m; ++i) {
		res.emplace_back(x, f(x));
		x += h;
	}
	return res;
}

vector<vector<double>> createKR(unsigned int n, vector<pair<double, double>> const &xs) {
	vector<vector<double>> res(n + 1);
	size_t m = xs.size();
	res[0] = vector<double>();
	for (auto const &pr : xs) {
		res[0].push_back(pr.second);
	}
	for (unsigned int i = 1; i <= n; ++i) {
		vector<double> kr(m - i);
		for (unsigned int j = 0; j < m - i; ++j) {
			kr[j] = res[i - 1][j + 1] - res[i - 1][j];
		}
		res[i] = kr;
	}
	return res;
}

double newton_start(double x, double a, unsigned int n, double h, vector<vector<double>> const &kr) {
	double t = (x - a) / h;
	double fact = 1;
	double mult = t;
	double res = kr[0][0];
	for (unsigned int i = 1; i <= n; ++i) {
		res += kr[i][0] * mult / fact;
		fact *= i + 1;
		mult *= t - i;
	}
	return res;
}

double newton_end(double x, double b, unsigned int n, double h, vector<vector<double>> const &kr) {
	double t = (x - b) / h;
	double fact = 1;
	double mult = t;
	double res = kr[0].back();
	for (unsigned int i = 1; i <= n; ++i) {
		res += kr[i].back() * mult / fact;
		fact *= i + 1;
		mult *= t + i;
	}
	return res;
}

double gauss_middle(unsigned int i0
					, double t
					, unsigned int n
					, vector<vector<double>> const &kr) {
	double fact = 1;
	double mult = t;
	double res = kr[0][i0];
	for (unsigned int i = 1; i <= n; ++i) {
		res += kr[i][i0 - (i / 2)] * mult / fact;
		fact *= i + 1;
		mult *= i % 2 == 0 ? t + ((i + 1) / 2) : t - ((i + 1) / 2);
	}
	return res;
}

bool enterxloop(function<double (double)> f
				, unsigned int n
				, unsigned int m
				, double a
				, double b
				, vector<pair<double, double>> &xs) {
	double x, h, lb, ub;
	h = (b - a) / m;
	lb = a + ((n + 1)/2) * h;
	ub = b - ((n + 1)/2) * h;
	cout << "Enter x from [" << a << ", " << a + h << "]U[" << b - h << ", " << b << "]" << endl;
	cout << "U[" << lb << ", " << ub << "]" << endl;
	cout << "x = ";
	cin >> x;
	while (x < a || x > b || (x > a + h && x < lb) || (x > ub && x < b - h)) {
		cout << "error: wrong x, enter it again" << endl;
		cout << "x = ";
		cin >> x;
	}
	auto kr = createKR(n, xs);
	double xl;
	if ((x >= a) & (x <= a + h)) {
		cout << endl << "Newton start: " << (xl = newton_start(x, a, n, h, kr)) << endl;
		cout << "|Pn(x) - f(x)| = " << fabs(xl - f(x)) << endl << endl;
	} else if ((x >= b - h) && (x <= b)) {
		cout << endl << "Newton end: " << (xl = newton_end(x, b, n, h, kr)) << endl;
		cout << "|Pn(x) - f(x)| = " << fabs(xl - f(x)) << endl << endl;
	} else {
		cout << endl << "Gauss middle: ";
		unsigned int i0 = (unsigned int) ((x - a) / h);
		cout << (xl = gauss_middle(i0, (x - xs[i0].first)/h, n, kr)) << endl;
		cout << "|Pn(x) - f(x)| = " << fabs(xl - f(x)) << endl << endl;
	}

	return true;
}

int main()
{
	double a, b;
	unsigned int m, n;
	auto f = [](double x){return (x*x / (1 + x*x));};
	//auto f = [](double x){return (x - 2)*(x - 4);};
	cout << "Enter [a, b]" << endl;
	cout << "a = ";
	cin >> a;
	cout << "b = ";
	cin >> b;
	while (b < a) {
		cout << "error: b < a, enter b again" << endl;
		cout << "b = ";
		cin >> b;
	}
	cout << "Enter m" << endl;
	cout << "m = ";
	cin >> m;
	cout << "Enter n (n <= " << m << ")" << endl;
	cout << "n = ";
	cin >> n;
	while (n > m) {
		cout << "error: n > m, enter n again" << endl;
		cout << "n = ";
		cin >> n;
	}
	auto xs = createTable(a, b, m, f);
	cout << "Table:" << endl;
	cout << setprecision(12) << std::fixed;
	cout << setw(15) << "x" << " | " << setw(15) << "f(x)" << endl;
	cout << setfill('-') << setw(35) << '-' << setfill(' ') << endl;
	for (auto const &pr : xs) {
		cout << setw(15) << pr.first << " | " << setw(15) << pr.second << endl;
	}

	while (enterxloop(f, n, m, a, b, xs)) {}

	return 0;
}
