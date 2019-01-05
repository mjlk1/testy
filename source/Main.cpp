#include <iostream>
#include <cinttypes>
#include <cmath>
#include <iomanip>
#include <vector>
#include <cstring>
#include <functional>

#include "Parameters.hpp"
#include "Euler1.hpp"
#include "Euler2.hpp"
#include "RungeKutta4.hpp"
#include "RungeKutta5.hpp"
#include "CashKarp.hpp"
#include "CashKarp0.hpp"

using namespace std;

int main(int argc, const char *argv[])
{
	Real time, h;
	int_fast32_t Nplot;
	Parameters par;
	string metoda;
	bool verbose = false;

	State r(4);
	r[0] = 6.0;
	r[1] = 1.0;

	for (int i=1;i<argc;++i)
	if (!strcmp(argv[i],"-v"))
		verbose = true;
	else
	{
		cerr << "nieznana opcja `" << argv[i] << "`" << '\n';
		return 1;
	}
	if (verbose)
		cerr << "metoda =";
	cin >> metoda;
	if (verbose)
		cerr << "g = ";
	cin >> par.g;
	if (verbose)
		cerr << "C = ";
	cin >> par.C;
	if (verbose)
		cerr << "n = ";
	cin >> par.n;
	if (verbose)
		cerr << "vx = ";
	cin >> r[2];
	if (verbose)
		cerr << "vy = ";
	cin >> r[3];
	if (verbose)
		cerr << "time = ";
	cin >> time;
	if (verbose)
		cerr << "Nplot = ";
	cin >> Nplot;


	vector<State> solution;
	solution.push_back(r);
	vector<Real> time_vec;
	time_vec.push_back(0);

	Real t = 0;

	if (metoda=="adaptywna")
	{
		Real eps;
		if (verbose)
			cerr << "epsilon = ";
		cin >> eps;
		if (verbose)
			cerr << "krok h = ";
		cin >> h;

		Real h_done = 0;

		r = cashkarp0(r,h,h_done,eps,par);
		t += h_done;

		if (t>time)
		{
			r = rk5(solution.back(),time,par);
			t = time;
		}
		solution.push_back(r);
		time_vec.push_back(t);
		while (t<time)
		{
			r = cashkarp(r,h,h_done,eps,par);
			t += h_done;

			if (t>time)
			{
				r = rk5(solution.back(),time-time_vec.back(),par);
				t = time;
			}
			solution.push_back(r);
			time_vec.push_back(t);
		}
	}
	else
	{
		int_fast32_t steps;
		if (verbose)
			cerr << "liczba krokow = ";
		cin >> steps;
		h = time/steps;

		function<State(const State &, const Real &, const Parameters &)> f;
		if (metoda=="euler1")
			f = euler1;
		else
		if (metoda=="euler2")
			f = euler2;
		else
		if (metoda=="rk4")
			f = rk4;
		else
		if (metoda=="rk5")
			f = rk5;
		else
		{
			cerr << "unknown method `" << metoda << "`\n";
			return 1;
		}

		for (int_fast32_t i=0;i<steps;++i)
		{
			t += h;
			r = f(r,h,par);
			solution.push_back(r);
			time_vec.push_back(t);
		}
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
