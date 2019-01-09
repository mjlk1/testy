#include <iostream>
#include <cinttypes>
#include <cmath>
#include <iomanip>
#include <vector>
#include <cstring>
#include <functional>
#include <immintrin.h>

#include "Parameters.hpp"
#include "Euler1.hpp"
#include "Euler2.hpp"
#include "RungeKutta4.hpp"
#include "RungeKutta5.hpp"
#include "CashKarp.hpp"
#include "CashKarp0.hpp"
#include "Energy.hpp"
#include "Clock.hpp"
#include "RungeKutta5AVX.hpp"

using namespace std;

int main(int argc, const char *argv[])
{
	Real time, h;
	int_fast32_t Nplot;
	Parameters par;
	string metoda;
	bool verbose = false;

	State r(4);
	r[0] = 1.0;
	r[1] = 0.0;

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

	StateAVX y;
	y[0]=r[0];
	y[1]=r[1];
	y[2]=r[2];
	y[3]=r[3];

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

		TimeInterval ti;
		Real h_done = 0;

		beginTimeMeasurement(ti);

		r = cashkarp(r,h,h_done,eps,par);
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

		endTimeMeasurement(ti);

		std::cerr << "done in " << timeIntervalToSeconds(ti) << " seconds\n";
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
		if (metoda!="rk5AVX" && metoda!="rk5avx2")
		{
			cerr << "unknown method `" << metoda << "`\n";
			return 1;
		}

		TimeInterval ti;
		if (metoda=="rk5AVX")
		{
			beginTimeMeasurement(ti);
			for (int_fast32_t i=0;i<steps;++i)
			{
				t += h;
				y = rk5AVX(y,h,par);
				r[0] = y[0];
				r[1] = y[1];
				r[2] = y[2];
				r[3] = y[3];
				solution.push_back(r);
				time_vec.push_back(t);
			}
			endTimeMeasurement(ti);
		}
		else
		if (metoda=="rk5avx2")
		{
			StateAVX *pSolution = (StateAVX*)aligned_alloc(sizeof(StateAVX),sizeof(StateAVX)*(steps+1));
			StateAVX *p = pSolution, *q = pSolution+1, *e = q+steps;
			*p = y;
			beginTimeMeasurement(ti);
			for ( ;q!=e;++p,++q)
			{
				*q = rk5AVXInline(*p,h,par);
			}
			endTimeMeasurement(ti);
			p = pSolution+1;
			for (int_fast32_t i=0;i<steps;++i)
			{
				t += h;
				r[0] = p[i][0];
				r[1] = p[i][1];
				r[2] = p[i][2];
				r[3] = p[i][3];
				solution.push_back(r);
				time_vec.push_back(t);
			}
			free(pSolution);
		}
		else
		{
			beginTimeMeasurement(ti);
			for (int_fast32_t i=0;i<steps;++i)
			{
				t += h;
				r = f(r,h,par);
				solution.push_back(r);
				time_vec.push_back(t);
			}
			endTimeMeasurement(ti);
		}

		std::cerr << "done in " << timeIntervalToSeconds(ti) << " seconds\n";
	}

	int_fast32_t steps = solution.size();
	if (Nplot==0)
		cout << setw(dist) << time << setw(dist) << solution[steps-1][0] << setw(dist) << solution[steps-1][1] << setw(dist) <<
		solution[steps-1][2] << setw(dist) << solution[steps-1][3] << setw(dist) << energy(solution[steps-1],par) << '\n';
	else
	{
		for (int_fast32_t i=0;i<steps;++i)
		{
			if (i%Nplot==0)
				cout << setw(dist) << time_vec[i] << setw(dist) << solution[i][0] << setw(dist) << solution[i][1] << setw(dist) <<
				solution[i][2] << setw(dist) << solution[i][3] << setw(dist) << energy(solution[i],par) << '\n';
		}
	}

	return 0;
}
