#ifndef TASK1_H
#define TASK1_H

#include <vector>
#include <functional>
#include <cmath>
#include <iostream>

namespace utils {
using namespace std;

bool check_zero(function<double (double)> f, double a, double b)
{
	return f(a) * f(b) < 0;
}

bool check_zero(function<double (double)> f, pair<double, double> const &p)
{
	return check_zero(f, p.first, p.second);
}

vector<pair<double, double>> find_intervals(function<double (double)> f, double start, double end, double n)
{
	//assert(start < end);
	double h = (end - start) / n;
	vector<pair<double, double>> intvals;

	double l = start;
	double r = start + h;
	while (r < end) {
		if (f(l) * f(r) < 0) {
			intvals.emplace_back(l, r);
		}
		l = r;
		r += h;
	}
	return intvals;
}

double aprox_newton(function<double (double)> f,
					function<double (double)> df,
					function<double (double)> ddf,
					pair<double, double> const &ival,
					double e)
{
	double x = ddf(ival.first) * f(ival.first) > 0 ? ival.first : ival.second;
	double x1 = x - f(x) / df(x);
	unsigned int n = 1;
	while (abs(x1 - x) > e && n < 500) {
		x = x1;
		x1 = x - f(x) / df(x);
		n++;
	}
	if (n == 500) {
		cout << "newton didn't converge fast enough" << endl;
	}
	cout << "newton took " << n << " steps" << endl;
	cout << "newton first aprox: " << x1 << endl;
	return x1;
}

double aprox_newton_mod(function<double (double)> f,
						function<double (double)> df,
						function<double (double)> ddf,
						pair<double, double> const &ival,
						double e)
{
	double x = ddf(ival.first) * f(ival.first) > 0 ? ival.first : ival.second;
	double dfx0 = df(x);
	double xk = x - f(x) / dfx0;
	double x1 = xk;
	unsigned int n = 1;
	while (abs(xk - x) > e && n < 500) {
		x = xk;
		xk = xk - f(xk) / dfx0;
		n++;
	}
	if (n == 500) {
		cout << "newton_mod didn't converge fast enough" << endl;
	}
	cout << "newton_mod took " << n << " steps" << endl;
	cout << "newton_mod first aprox: " << x1 << endl;
	return xk;
}

double aprox_binsearch(function<double (double)> f, double start, double end, double e, unsigned int n = 0) {
	n++;
	if (end - start < e) {
		cout << "bisect took " << n << " steps" << endl;
		return start;
	}
	double mid = (start + end) / 2;
	if (f(start) * f(mid) < 0) {
		return aprox_binsearch(f, start, mid, e, n);
	} else {
		return aprox_binsearch(f, mid, end, e, n);
	}
}

double aprox_chord(function<double (double)> f, double start, double end, double e) {
	double x1 = start;
	double x2 = end;
	unsigned int n = 1;
	double xk;
	double first = x1 - f(x1) * (x2 - x1) / (f(x2) - f(x1));
	while (abs(x1 - x2) > e && n < 500) {
		xk = x1 - f(x1) * (x2 - x1) / (f(x2) - f(x1));
		x1 = x2;
		x2 = xk;
		n++;
	}
	if (n == 500) {
		cout << "chord didn't converge fast enough" << endl;
	}
	cout << "chord took " << n << " steps" << endl;
	cout << "chord first aprox: " << first << endl;
	return xk;
}

}

#endif // TASK1_H

