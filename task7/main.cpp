#include <iostream>
#include <vector>
#include <functional>
#include <cmath>
#include <iomanip>
#include <tuple>

using namespace std;

void calculateKR(vector<vector<double>> &xs) {
	for (unsigned int i = 3; i < 7; ++i) {
		for (unsigned int j = 0; j < xs[i - 1].size() - 1; ++j) {
			xs[i].push_back(xs[i - 1][j + 1] - xs[i - 1][j]);
		}
	}
}

int main()
{
	double h;
	unsigned int N;
	double x0 = 0;
	double y0 = 1;
	auto f = [](double x, double y){return (-2*y+y*y);};
	auto y = [](double x){return (2/(exp(2*x) + 1));};
	auto yT = [](double x){return (1 - x + x*x*x/3 - 2/15*pow(x,5) + 17/315*pow(x,7));};
	while (true) {
		cout << "N = ";
		cin >> N;
		cout << "Enter h:" << endl;
		cout << "h = ";
		cin >> h;
		cout << "Solution table:" << endl;
		{
			cout << setprecision(12) << fixed;
			cout << setw(15) << "xi" << " | " << setw(15) << "f(x)" << endl;
			cout << setfill('-') << setw(35) << '-' << setfill(' ') << endl;
			double xi = x0 - 2*h;
			for (unsigned int i = 0; i < N + 3; ++i) {
				cout << setw(15) << xi << " | " << setw(15) << y(xi) << endl;
				xi += h;
			}
		}
		vector<vector<double>> table(N + 2);
		cout << endl << "Taylor table:" << endl;
		{
			cout << setprecision(12) << fixed;
			cout << setw(15) << "xi" << " | " << setw(15) << "f(x)";
			cout << " | " << setw(15) << "|y(x) - yT(x)|" << endl;
			cout << setfill('-') << setw(52) << '-' << setfill(' ') << endl;
			double xi = x0 - 2*h;
			for (unsigned int i = 0; i < 5; ++i) {
				table[0].push_back(xi);
				table[1].push_back(yT(xi));
				table[2].push_back(h*f(xi, yT(xi)));
				cout << setw(15) << xi << " | " << setw(15) << table[1][i];
				cout << " | " << setw(15) << fabs(y(xi) - table[1][i]) << endl;
				xi += h;
			}
		}
		calculateKR(table);
		cout << endl << "Adams table:" << endl;
		{
			cout << setprecision(12) << fixed;
			cout << setw(15) << "xi" << " | " << setw(15) << "f(x)";
			cout << " | " << setw(15) << "|y(x) - yT(x)|" << endl;
			cout << setfill('-') << setw(52) << '-' << setfill(' ') << endl;
			double xi = x0 + 3*h;
			for (unsigned int i = 5; i < N + 3; ++i) {
				double dyi = table[2][i - 1] + 0.5 * table[3][i - 2] + 5.0/12.0*table[4][i - 3]
						+ 3.0/8.0*table[5][i - 4] + 251.0/720.0*table[6][i - 5];
				double yi = table[1][i - 1] + dyi;
				cout << setw(15) << xi << " | " << setw(15) << yi;
				table[0].push_back(xi);
				table[1].push_back(yi);
				table[2].push_back(h*f(xi,yi));

				for (unsigned int j = 0; j < 5; ++j) {
					table[j + 3].push_back(table[j + 2][i - j] - table[j + 2][i - j - 1]);
				}
				cout << " | " << setw(15) << fabs(y(xi) - table[1][i]) << endl;
				xi += h;
			}
		}
		cout << endl << "R-K table:" << endl;
		{
			cout << setprecision(12) << fixed;
			cout << setw(15) << "xi" << " | " << setw(15) << "f(x)";
			cout << " | " << setw(15) << "|y(x) - yT(x)|" << endl;
			cout << setfill('-') << setw(52) << '-' << setfill(' ') << endl;
			double xi = x0;
			double yi = y0;
			double k1, k2, k3, k4;
			for (unsigned int i = 1; i <= N; ++i) {
				k1 = h*f(xi, yi);
				k2 = h*f(xi + h/2, yi + k1/2);
				k3 = h*f(xi + h/2, yi + k2/2);
				k4 = h*f(xi + h, yi + k3);
				yi += (k1 + 2*k2 + 2*k3 + k4)/6;
				xi += h;
				cout << setw(15) << xi << " | " << setw(15) << yi;
				cout << " | " << setw(15) << fabs(y(xi) - yi) << endl;
			}
		}
		cout << endl << "E table:" << endl;
		{
			cout << setprecision(12) << fixed;
			cout << setw(15) << "xi" << " | " << setw(15) << "f(x)";
			cout << " | " << setw(15) << "|y(x) - yT(x)|" << endl;
			cout << setfill('-') << setw(52) << '-' << setfill(' ') << endl;
			double xi = x0;
			double yi = y0;
			for (unsigned int i = 1; i <= N; ++i) {
				yi += h*f(xi,yi);
				xi += h;
				cout << setw(15) << xi << " | " << setw(15) << yi;
				cout << " | " << setw(15) << fabs(y(xi) - yi) << endl;
			}
		}
		cout << endl << "E+ table:" << endl;
		{
			cout << setprecision(12) << fixed;
			cout << setw(15) << "xi" << " | " << setw(15) << "f(x)";
			cout << " | " << setw(15) << "|y(x) - yT(x)|" << endl;
			cout << setfill('-') << setw(52) << '-' << setfill(' ') << endl;
			double xi = x0;
			double yi = y0;
			for (unsigned int i = 1; i <= N; ++i) {
				yi += h*f(xi + h/2,yi + h/2*f(xi,yi));
				xi += h;
				cout << setw(15) << xi << " | " << setw(15) << yi;
				cout << " | " << setw(15) << fabs(y(xi) - yi) << endl;
			}
		}
		cout << endl << "E-K table:" << endl;
		{
			cout << setprecision(12) << fixed;
			cout << setw(15) << "xi" << " | " << setw(15) << "f(x)";
			cout << " | " << setw(15) << "|y(x) - yT(x)|" << endl;
			cout << setfill('-') << setw(52) << '-' << setfill(' ') << endl;
			double xi = x0;
			double yi = y0;
			for (unsigned int i = 1; i <= N; ++i) {
				yi += h/2*(f(xi,yi) + f(xi + h,yi + h*f(xi,yi)));
				xi += h;
				cout << setw(15) << xi << " | " << setw(15) << yi;
				cout << " | " << setw(15) << fabs(y(xi) - yi) << endl;
			}
		}
	}

	return 0;
}

