#include <iostream>
#include <vector>
#include <functional>
#include <cmath>
#include <iomanip>
#include <tuple>

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

double approx_binsearch(double start, double end, double e, function<double (double)> f) {
	if (end - start < e) {
		return start;
	}
	double mid = (start + end) / 2;
	if (f(start) * f(mid) < 0) {
		return approx_binsearch(start, mid, e, f);
	} else {
		if (f(start) == 0) {
			return start;
		} else if (f(end) == 0) {
			return end;
		}
		return approx_binsearch(mid, end, e, f);
	}
}

typedef tuple<double, double, double> vec3;
double calculateDet3(vec3 row1, vec3 row2, vec3 row3) {
	return (get<0>(row1)*get<1>(row2)*get<2>(row3)
			+ get<1>(row1)*get<2>(row2)*get<0>(row3)
			+ get<2>(row1)*get<0>(row2)*get<1>(row3)
			- get<2>(row1)*get<1>(row2)*get<0>(row3)
			- get<1>(row1)*get<0>(row2)*get<2>(row3)
			- get<0>(row1)*get<2>(row2)*get<1>(row3));
}

vector<double> calculateMus(double a, double b, unsigned int n, function<double (double)> w) {
	vector<double> res(n);
	for (size_t i = 0; i < n; ++i) {
		auto fw = [i, &w](double x){return (pow(x, i)*w(x));};
		res[i] = rectangleKF(a, b, 256, fw);
	}
	return res;
}

function<double (double)> findOmega(const vector<double> &mus) {
	double det = calculateDet3(make_tuple(mus[2], mus[1], mus[0])
							   , make_tuple(mus[3], mus[2], mus[1])
							   , make_tuple(mus[4], mus[3], mus[2]));
	double det1 = calculateDet3(make_tuple(-mus[3], mus[1], mus[0])
							   , make_tuple(-mus[4], mus[2], mus[1])
							   , make_tuple(-mus[5], mus[3], mus[2]));
	double det2 = calculateDet3(make_tuple(mus[2], -mus[3], mus[0])
							   , make_tuple(mus[3], -mus[4], mus[1])
							   , make_tuple(mus[4], -mus[5], mus[2]));
	double det3 = calculateDet3(make_tuple(mus[2], mus[1], -mus[3])
							   , make_tuple(mus[3], mus[2], -mus[4])
							   , make_tuple(mus[4], mus[3], -mus[5]));
	cout << "c = " << det1/det << ", d = " << det2/det << ", p = " << det3/det << endl;
	auto om = [det, det1, det2, det3](double x){
		return (x*x*x + det1/det*x*x + det2/det*x + det3/det);
	};
	return {om};
}

vector<double> findOmegaRoots(double a, double b, function<double (double)> om) {
	vector<double> res;
	unsigned int m = 16 * static_cast<unsigned int>(b - a + 1);
	double h = (b - a)/m;
	double start = a;
	for (size_t i = 0; i < m; ++i) {
		if (om(start) * om(start + h) < 0) {
			res.push_back(approx_binsearch(start, start + h, 0.00000001/*E-8*/, om));
		}
		start += h;
	}
	return res;
}

vector<double> findAs(const vector<double> &xs, const vector<double> &mus) {
	double det = calculateDet3(make_tuple(1, 1, 1)
							   , make_tuple(xs[0], xs[1], xs[2])
							   , make_tuple(xs[0]*xs[0], xs[1]*xs[1], xs[2]*xs[2]));
	double det1 = calculateDet3(make_tuple(mus[0], 1, 1)
								, make_tuple(mus[1], xs[1], xs[2])
								, make_tuple(mus[2], xs[1]*xs[1], xs[2]*xs[2]));
	double det2 = calculateDet3(make_tuple(1, mus[0], 1)
								, make_tuple(xs[0], mus[1], xs[2])
								, make_tuple(xs[0]*xs[0], mus[2], xs[2]*xs[2]));
	double det3 = calculateDet3(make_tuple(1, 1, mus[0])
							   , make_tuple(xs[0], xs[1], mus[1])
							   , make_tuple(xs[0]*xs[0], xs[1]*xs[1], mus[2]));
	return {det1/det, det2/det, det3/det};
}

int main()
{
	double a, b;
	const unsigned int N = 3;
	auto f = [](double x){return sin(x);};
	auto w = [](double x){return (-log(x));};
	while (true) {
		cout << "N = " << N << endl;
		cout << "Enter [a, b]:" << endl;
		cout << "a = ";
		cin >> a;
		cout << "b = ";
		cin >> b;
		while (b < a) {
			cout << "error: b < a, enter b again" << endl;
			cout << "b = ";
			cin >> b;
		}
		cout.precision(10);
		cout << fixed;
		auto mus = calculateMus(a, b, 2*N, w);
		for (size_t i = 0; i < mus.size(); ++i) {
			cout << "mu" << i << " = " << mus[i] << endl;
		}
		cout << endl;
		auto xs = findOmegaRoots(a, b, findOmega(mus));
		for (size_t i = 0; i < xs.size(); ++i) {
			cout << "x" << i << " = " << xs[i] << endl;
		}
		if (xs.size() < 3) {
			cout << "Not enough w roots because mus are not precise enough" << endl << endl;
			continue;
		}
		cout << endl;
		auto as = findAs(xs, mus);
		for (size_t i = 0; i < as.size(); ++i) {
			cout << "A" << i << " = " << as[i] << endl;
		}
		cout << endl << "check: " << endl;
		cout << "mu5 = " << mus[5] << endl;
		cout << "A1*x1^5 + A2*x2^5 + A3*x3^5 = ";
		cout << (as[0]*pow(xs[0],5) + as[1]*pow(xs[1],5) + as[2]*pow(xs[2],5)) << endl;

		cout << endl << "A1*f(x1) + A2*f(x2) + A3*f(x3) = ";
		cout << (as[0]*f(xs[0]) + as[1]*f(xs[1]) + as[2]*f(xs[2])) << endl << endl;
	}

	return 0;
}

