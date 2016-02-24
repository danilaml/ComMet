#include <iostream>
#include <iomanip>
#include <vector>
#include <functional>
#include <cmath>

using namespace std;

vector<double> solveA(const vector<double> &a,
					  const vector<double> &b,
					  const vector<double> &c,
					  const vector<double> &d,
					  bool printCoeffs = true)
{
	auto n = a.size() - 1;
	vector<double> m(a.size());
	vector<double> k(a.size());
	m[0] = - c[0]/b[0];
	k[0] = d[0]/b[0];

	if (printCoeffs)
		cout << "coefficients:" << endl;
	for (unsigned i = 0; i < n; i++) {
		m[i + 1] = - c[i]/(a[i]*m[i] + b[i]);
		if (printCoeffs)
			cout << "m[" << (i + 1) << "] = " << m[i + 1] << endl;
		k[i + 1] = (d[i] - a[i]*k[i])/(a[i]*m[i] + b[i]);
		if (printCoeffs)
			cout << "k[" << (i + 1) << "] = " << k[i + 1] << endl;
	}

	vector<double> x(a.size());
	cout << "x:" << endl;
	x[n] = (d[n] - a[n]*k[n])/(a[n]*m[n] + b[n]);
	cout << "x[" << n << "] = " << x[n] << endl;
	for (int i = n - 1; i >= 0; --i) {
		x[i] = m[i + 1]*x[i + 1] + k[i + 1];
		cout << "x[" << i << "] = " << x[i] << endl;
	}

	return x;
}

vector<double> solveDif(double alpha0, double alpha1, double A,
						double beta0, double beta1, double B,
						double left, double right, int n, bool isFirst,
						function<double (double)> p,
						function<double (double)> q,
						function<double (double)> r,
						function<double (double)> f)
{
	double h = (right - left)/n;
	vector<double> a(n + 1);
	vector<double> b(n + 1);
	vector<double> c(n + 1);
	vector<double> d(n + 1);
	a[0] = 0;
	c[n] = 0;

	double xi = alpha0 + h;
	for (int i = 1; i < n; ++i) {
		a[i] = p(xi) - q(xi)*h/2.0;
		b[i] = -2.0*p(xi) + r(xi)*h*h;
		c[i] = p(xi) + q(xi)*h/2.0;
		d[i] = f(xi)*h*h;
	}
	if (isFirst) {
		b[0] = h*alpha0 - alpha1;
		c[0] = alpha1;
		d[0] = h*A;
		a[n] = - beta1;
		b[n] = h*beta0 + beta1;
		d[n] = h*B;
	} else {
		b[0] = 2.0*h*alpha0+alpha1*(a[1]/c[1] - 3.0);
		c[0] = alpha1*(b[1]/c[1] + 4.0);
		d[0] = 2.0*h*A + alpha1*d[1]/c[1];
		a[n] = - beta1*(4.0 + b[n - 1]/a[n - 1]);
		b[n] = 2.0*h*beta0 + beta1*(3.0 - c[n - 1]/a[n - 1]);
		d[n] = 2.0*h*B - beta1*d[n - 1]/a[n - 1];
	}

	return solveA(a, b, c, d, false);
}

int main()
{
	cout << fixed;
	cout.precision(10);
	unsigned n;
	do {
		cout << "n = ";
		cin >> n;
	} while (n < 1);
	vector<double> a(n + 1);
	vector<double> b(n + 1);
	vector<double> c(n + 1);
	vector<double> d(n + 1);

	a[0] = 0;
	c[n] = 0;
	cout << endl;
	for (unsigned i = 1; i <= n; ++i) {
		cout << "a[" << i << "] = ";
		cin >> a[i];
	}
	cout << endl;
	for (unsigned i = 0; i <= n; ++i) {
		cout << "b[" << i << "] = ";
		cin >> b[i];
	}
	cout << endl;
	for (unsigned i = 0; i < n; ++i) {
		cout << "c[" << i << "] = ";
		cin >> c[i];
	}
	cout << endl;
	for (unsigned i = 0; i <= n; ++i) {
		cout << "d[" << i << "] = ";
		cin >> d[i];
	}
	cout << endl << endl;

	auto x = solveA(a, b, c, d);

	cout << "r[" << 0 << "] = " << (b[0]*x[0] + c[0]*x[1] - d[0]) << endl;
	for (unsigned i = 1; i < n; ++i) {
		cout << "r[" << i << "] = " << (a[i]*x[i - 1] + b[i]*x[i] + c[i]*x[i + 1] - d[i]) << endl;
	}
	cout << "r[" << n << "] = " << (a[n]*x[n - 1] + b[n]*x[n] - d[n]) << endl;

	cout << endl;
	cout << "--------------DiffEq--------------";
	auto p = [](double x){return 1.0;};
	auto q = [](double x){return (1.0 + 2.0*x);};
	auto r = [](double x){return -log(1.0 + x);};
	auto f = [](double x){return (x - 1.0);};
	double alpha0 = 0.5;
	double alpha1 = -1.0;
	double beta0 = -0.7;
	double beta1 = -1.0;
	double A = 0;
	double B = 0;
	cout << endl << "First, n = 100 :" << endl;
	solveDif(alpha0, alpha1, A, beta0, beta1, B, 0, 1, 100, true, p, q, r, f);
	cout << endl << "Second, n = 10 :" << endl;
	solveDif(alpha0, alpha1, A, beta0, beta1, B, 0, 1, 10, false, p, q, r, f);

	return 0;
}

