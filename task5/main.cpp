#include <iostream>
#include <vector>
#include <functional>
#include <cmath>
#include <iomanip>

using namespace std;

double rectangleKF(double a, double b, unsigned int m, function<double (double)> f) {
	double h = (b - a)/m;
	double xk = a;
	double res = 0;
	for (unsigned int i = 0; i < m; ++i) {
		res += f(xk + h/2);
		xk += h;
	}
	res *= h;
	return res;
}

double trapezeKF(double a, double b, unsigned int m, function<double (double)> f) {
	double h = (b - a)/m;
	double xk = a;
	double res = 0;
	for (unsigned int i = 0; i < m; ++i) {
		res += f(xk) + f(xk + h);
		xk += h;
	}
	res *= h/2;
	return res;
}

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

int main()
{
	double a, b, J, Jh;
	unsigned int m;
	//auto f = [](double x){return (x*x*x + x*x);};
	//auto Jf = [](double x){return (x*x*x*x/4 + x*x*x/3);};
	auto f = [](double x){return (1 - exp(-x) + x*x);};
	auto Jf = [](double x){return (x + exp(-x) + x*x*x/3);};
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
	cout.precision(15);
	cout << fixed;
	cout << "J = " << (J = Jf(b) - Jf(a)) << endl << endl;
	cout << "rectangleKF = " << (Jh = rectangleKF(a, b, m, f)) << endl;
	cout << "|J - Jh| = " << fabs(J - Jh) << endl;
	cout << "trapezeKF = " << (Jh = trapezeKF(a, b, m, f)) << endl;
	cout << "|J - Jh| = " << fabs(J - Jh) << endl;
	cout << "simpsonKF = " << (Jh = simpsonKF(a, b, m, f)) << endl;
	cout << "|J - Jh| = " << fabs(J - Jh) << endl;

	return 0;
}

