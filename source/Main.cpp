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

using namespace std;

int main(int argc, const char *argv[])
{
	Real time;
	int_fast32_t steps, Nplot;
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
		cerr << "Liczba krokow = ";
	cin >> steps;
	Real h = time/steps;
	if (verbose)
	{
		cerr << "Dlugosc kroku wynosi" << '\t' << h << '\n';
		cerr << "Nplot = ";
	}
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

	for (int_fast32_t i=0;i<steps;++i)
	{
		r = rk4(r,h,par);
		solution.push_back(r);
	}

	if (Nplot==0)
		cout << setw(dist) << time << setw(dist) << solution[steps][0] << setw(dist) << solution[steps][1] << setw(dist) << 
		solution[steps][2] << setw(dist) << solution[steps][3] << '\n';
	else
	{
		for (int_fast32_t i=0;i<=steps;++i)
		{
			if (i%Nplot==0)
				cout << setw(dist) << i*h << setw(dist) << solution[i][0] << setw(dist) << solution[i][1] << setw(dist) <<
				solution[i][2] << setw(dist) << solution[i][3] << '\n';
		}
	}

	return 0;
}
