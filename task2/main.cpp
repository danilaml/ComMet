#include <iostream>
#include <functional>
#include <vector>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <cassert>

using namespace std;

vector<pair<double, double>> createTable(double a, double b, unsigned int m, function<double (double)> f) {
	vector<pair<double, double>> res;
	double h = (b - a) / m;
	double x = a;
	for (unsigned int i = 0; i < m; ++i) {
		res.emplace_back(x, f(x));
		x += h;
	}
	return res;
}

void sortTable(double x, vector<pair<double, double>> &xs) {
	sort(xs.begin(), xs.end(),
		 [x](pair<double, double> const &x1, pair<double, double> &x2){return fabs(x - x2.first) - fabs(x - x1.first) > 0;});
}

double newton(double x, unsigned int n, vector<pair<double, double>> const &xs) {
	double px = xs[0].second;
	double xms = 1;
	vector<double> fxs;
	fxs.reserve(n + 1);
	for (unsigned int i = 0; i <= n; ++i) {
		fxs[i] = xs[i].second;
	}
	for (unsigned int i = 0; i < n; ++i) {
		for (unsigned int j = 0; j < n - i; ++j) {
			fxs[j] = (fxs[j + 1] - fxs[j])/(xs[j + 1 + i].first - xs[j].first);
		}
		xms *= x - xs[i].first;
		px += fxs[0] * xms;
	}
	return px;
}

double lagrange(double x, unsigned int n, vector<pair<double, double>> const &xs) {
	double px = 0;
	double li = 1;
	for (unsigned int i = 0; i <= n; ++i) {
		for (unsigned int j = 0; j <= n; ++j) {
			if (i != j)
				li *= (x - xs[j].first)/(xs[i].first - xs[j].first);
		}
		px += xs[i].second * li;
		li = 1;
	}
	return px;
}

int main()
{
	double a, b, x;
	unsigned int m, n;
	auto f = [](double x){return (1 - exp(-x) + x*x);};
	//auto f = [](double x){return (x - 2)*(x - 4);};
	cout << "Enter [a, b]" << endl;
	cout << "a = "; // TODO: check for correctness
	cin >> a;
	cout << "b = ";
	cin >> b;
	if (b < a) {
		cout << "error: b < a, enter b again" << endl;
		cout << "b = ";
		cin >> b;
	}
	cout << "Enter m" << endl;
	cout << "m = ";
	cin >> m; m++;
	cout << "Enter n (n < " << m << ")" << endl;
	cout << "n = ";
	cin >> n;
	if (n > m) {
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
	cout << "Enter x" << endl;
	cout << "x = ";
	cin >> x;
	sortTable(x, xs);
	cout << "Sorted table:" << endl;
	cout << setw(15) << "x" << " | " << setw(15) << "f(x)" << endl;
	cout << setfill('-') << setw(35) << '-' << setfill(' ') << endl;
	for (auto const &pr : xs) {
		cout << setw(15) << pr.first << " | " << setw(15) << pr.second << endl;
	}
	double xl;
	cout << endl << "Lagrange: " << (xl = lagrange(x, n, xs)) << endl;
	cout << "|Pn(x) - f(x)| = " << fabs(xl - f(x)) << endl;
	cout << endl << "Newton: " << (xl = newton(x, n, xs)) << endl;
	cout << "|Pn(x) - f(x)| = " << fabs(xl - f(x)) << endl;
	return 0;
}

