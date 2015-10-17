#include <iostream>
#include <functional>
#include <vector>
#include <cmath>
#include <iomanip>
#include <algorithm>

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

void sortTable(double x, vector<pair<double, double>> &xs) {
	sort(xs.begin(), xs.end(),
		 [x](pair<double, double> const &x1, pair<double, double> &x2){
		return fabs(x - x2.first) - fabs(x - x1.first) > 0;
	});
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

double aprox_binsearch(function<double (double)> f, double start, double end, double e) {
	if (end - start < e) {
		return start;
	}
	double mid = (start + end) / 2;
	if (f(start) * f(mid) < 0) {
		return aprox_binsearch(f, start, mid, e);
	} else {
		if (f(start) == 0) {
			return start;
		} else if (f(end) == 0) {
			return end;
		}
		return aprox_binsearch(f, mid, end, e);
	}
}

bool enterxloop(function<double (double)> f
				, unsigned int m
				, double fu
				, double fl
				, vector<pair<double, double>> &xs
				, vector<pair<double, double>> &xsrev) {
	double n, F;
	cout << "Enter n (n <= " << m << ")" << endl;
	cout << "n = ";
	cin >> n;
	while (n > m) {
		cout << "error: n > m, enter n again" << endl;
		cout << "n = ";
		cin >> n;
	}
	cout << "Enter F" << endl;
	cout << "F = ";
	cin >> F;
	sortTable(F, xsrev);
	double X;
	cout << endl << "Method 1: X = " << (X = lagrange(F, n, xsrev)) << endl;
	cout << "|f(X) - F| = " << fabs(f(X) - F) << endl;
	if (!(F <= fu && F >= fl)) {
		return true;
	}
	size_t j = 0;
	double avg = 0;
	for (size_t i = 0; i < xs.size() - 1; ++i) {
		if ((xs[i].second - F) * (xs[i + 1].second - F) <= 0) {
			++j;
			avg += (xs[i].first + xs[i + 1].first)/2;
		}
	}
	vector<pair<double, double>> xscopy(xs);
	sortTable(avg/j, xscopy);
	auto p = [n, F, &xscopy](double x){
		return (lagrange(x, n, xscopy) - F);
	};
	cout << endl << "Method 2: " << endl;
	j = 0;
	for (size_t i = 0; i < xs.size() - 1; ++i) {
		if (p(xs[i].first) * p(xs[i + 1].first) <= 0) {
			cout << "x" << j << " = ";
			cout << (X = aprox_binsearch(p, xs[i].first, xs[i + 1].first, 0.00000001)) << endl;
			cout << "|f(" << "x" << j++ <<") - F| = " << fabs(f(X) - F) << endl;
		}
	}
	cout << endl;

	return true;
}

void derivatives(double h
				, vector<pair<double, double>> const &xs
				, vector<double> &fds
				, vector<double> &fdds) {
	fds.reserve(xs.size());
	fdds.reserve(xs.size());
	fds[0] = ((-3)*xs[0].second + 4*xs[1].second - xs[2].second)/(2*h);
	fdds[0] = nan("");
	fds[xs.size() - 1] = (3*xs.back().second - 4*xs[xs.size() - 2].second + xs[xs.size() - 3].second)/(2*h);
	fdds[xs.size() - 1] = nan("");
	for (size_t i = 1; i < xs.size() - 1; ++i) {
		fds[i] = (xs[i + 1].second - xs[i - 1].second)/(2*h);
		fdds[i] = (xs[i + 1].second - 2*xs[i].second + xs[i - 1].second)/h/h;
	}
}

int main()
{
	double a, b;
	unsigned int m;
	auto f = [](double x){return (1 - exp(-x) + x*x);};
	auto df = [](double x){return (exp(-x) + 2*x);};
	auto ddf = [](double x){return (-exp(-x) + 2);};
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
	auto xs = createTable(a, b, m, f);
	cout << "Table:" << endl;
	cout << setprecision(12) << std::fixed;
	cout << setw(15) << "x" << " | " << setw(15) << "f(x)" << endl;
	cout << setfill('-') << setw(35) << '-' << setfill(' ') << endl;
	for (auto const &pr : xs) {
		cout << setw(15) << pr.first << " | " << setw(15) << pr.second << endl;
	}
	vector<double> fd, fdd;
	double h = (b - a) / m;
	derivatives(h, xs, fd, fdd);
	cout << "Derivatives:" << endl;
	cout << setprecision(12) << std::fixed;
	cout << setw(15) << "x" << " | " << setw(15) << "f(x)";
	cout << " | " << setw(15) << "fn'(x)" << " | " << setw(16) << "|f'(x) - fn'(x)|";
	cout << " | " << setw(15) << "fn''(x)" << " | " << setw(15) << "|f''(x) - fn'(x)|";
	cout << endl << setfill('-') << setw(36 * 3) << '-' << setfill(' ') << endl;
	for (size_t i = 0; i < xs.size(); ++i) {
		cout << setw(15) << xs[i].first << " | " << setw(15) << xs[i].second << " | ";
		cout << setw(15) << fd[i] << " | " << setw(16) << fabs(fd[i] - df(xs[i].first)) << " | ";
		cout << setw(15) << fdd[i] << " | " << setw(15) << fabs(fdd[i] - ddf(xs[i].first)) << endl;
	}

	vector<pair<double, double>> xsrev;
	double fupper, flower;
	fupper = xs[0].second;
	flower = fupper;
	for (auto const &pr : xs) {
		fupper = fupper < pr.second ? pr.second : fupper;
		flower = flower > pr.second ? pr.second : flower;
		xsrev.emplace_back(pr.second, pr.first);
	}
	while (enterxloop(f, m, fupper, flower, xs, xsrev)) {}

	return 0;
}
