#include <iostream>
#include <iomanip>
#ifdef _WIN32
#include <clocale>
#endif
#include "task1.h"

void print_fxn(std::function<double (double)> f, double xn)
{
	std::cout << "|f(x) - 0| = " << std::abs(f(xn)) << std::endl;
}

int main()
{
	using namespace utils;
#ifdef _WIN32
	system("chcp 65001");
#endif
	double e = 0.00001;
	double a = -20;
	double b = 20;
	cout << setprecision(15);
	auto f = [](double x){return (x*x - sin(x));};
	auto df = [](double x){return (2*x - cos(x));};
	auto ddf = [](double x){return (2 + sin(x));};
	//auto f = [](double x){return (x*x*x*x + 1/(5 - x) + 1);};
	//auto df = [](double x){return (4*x*x*x + 1/(5 - x)/(5 - x));};
	//auto ddf = [](double x){return (-1);};
	auto res = find_intervals(f, a, b, 100);
	cout << " Численные методы решения нелинейных алгебраических и" << endl;
	cout << " трансцендентных уравнений" << endl;
	cout << "A := " << a << ", B := " << b << ", e := " << e;
	cout << ", h := " << (b-a)/100 << endl << endl;
	double xn;
	for (auto const &pr : res) {
		cout << "interval: " << pr.first << ", " << pr.second << endl;
		cout << "=============" << endl;
		xn = aprox_binsearch(f, pr.first, pr.second, e);
		cout << "bisect first aprox: " << (pr.first + pr.second) / 2 << endl;
		cout << "bisect: " << xn << endl;
		print_fxn(f, xn);
		cout << endl;
		if (check_zero(df, pr) || check_zero(ddf, pr)) {
			cout << "newton is not applicable (f' or f'' is 0)" << endl;
		}
		else {
			xn = aprox_newton(f, df, ddf, pr, e);
			cout << "newton: " << xn << endl;
			print_fxn(f, xn);
		}
		cout << endl;
		xn = aprox_newton_mod(f, df, ddf, pr, e);
		cout << "newton_mod: " << xn << endl;
		print_fxn(f, xn);
		cout << endl;
		xn = aprox_chord(f, pr.first, pr.second, e);
		cout << "chord: " << xn << endl;
		print_fxn(f, xn);
		cout << endl << endl;
	}
	return 0;
}

