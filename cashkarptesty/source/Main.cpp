#include <iostream>
#include <cinttypes>
#include <cmath>
#include <iomanip>
#include <vector>
#include <cstring>

#include "Parameters.hpp"
#include "Euler1.hpp"
#include "Euler2.hpp"
#include "RungeKutta4.hpp"
#include "RungeKutta5.hpp"
#include "CashKarp.hpp"

using namespace std;

int main(int argc, const char *argv[])
{
	Real time, eps, h;
	int_fast32_t Nplot;
	Parameters par;
	bool verbose = false;

	for (int i=1;i<argc;++i)
	if (!strcmp(argv[i],"-v"))
		verbose = true;
	else
	{
		cerr << "nieznana opcja `" << argv[i] << "`" << '\n';
		return 1;
	}
	if (verbose)
		cerr << "g = ";
	cin >> par.g;
	if (verbose)
		cerr << "C = ";
	cin >> par.C;
	if (verbose)
		cerr << "mass = ";
	cin >> par.mass;
	if (verbose)
		cerr << "n = ";
	cin >> par.n;
	if (verbose)
		cerr << "time = ";
	cin >> time;
	if (verbose)
		cerr << "epsilon = ";
	cin >> eps;
		if (verbose)
		cerr << "h poczatkowe= ";
	cin >> h;
	if (verbose)
		cerr << "Nplot = ";
	cin >> Nplot;

	State r(4);
	r[0] = 0.0;
	r[1] = 0.0;

	if (verbose)
		cerr << "vx = ";
	cin >> r[2];
	if (verbose)
		cerr << "vy = ";
	cin >> r[3];

	vector<State> solution;
	solution.push_back(r);
	vector<Real> time_vec;
	time_vec.push_back(0);

	Real t = 0;

	while (t<time)
	{
		t += h;
		r =cash_karp(r,h,eps,par);
		solution.push_back(r);
		time_vec.push_back(t);
	}
	int_fast32_t steps = solution.size();
	if (Nplot==0)
		cout << setw(dist) << time << setw(dist) << solution[steps-1][0] << setw(dist) << solution[steps-1][1] << setw(dist) << 
		solution[steps-1][2] << setw(dist) << solution[steps-1][3] << '\n';
	else
	{
		for (int_fast32_t i=0;i<steps;++i)
		{
			if (i%Nplot==0)
				cout << setw(dist) << time_vec[i] << setw(dist) << solution[i][0] << setw(dist) << solution[i][1] << setw(dist) <<
				solution[i][2] << setw(dist) << solution[i][3] << '\n';
		}
	}

	return 0;
}