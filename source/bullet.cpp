#include <iostream>
#include <cinttypes>
#include <cmath>
#include <iomanip>
#include <vector>
#include <cstring>

using namespace std;

typedef float Real;
typedef vector<Real> State; //wektor stanu (x,y,vx,vy)
struct Parametry
{
	Real g,C;
};

const Real mass = 1.0;
const Real n = 2.0;
const int_fast32_t  dist = 2*sizeof(Real)+3;

State rhs(const State &r, const Parametry &par)
{
	State deriv;
	Real l = -par.C*pow(r[2]*r[2]+r[3]*r[3], (n-1)/2)/mass;
	deriv.push_back(r[2]);
	deriv.push_back(r[3]);
	deriv.push_back(l*r[2]);
	deriv.push_back(l*r[3]+par.g);
	return deriv;
}

State euler1(const State &r, const Real &h, const Parametry &par)
{
	State sol;
	State s = rhs(r,par);
	for (int_fast32_t i=0;i<4;++i)
		sol.push_back(r[i]+h*s[i]);
	return sol;
}

State euler2(const State &r, const Real &h, const Parametry &par)
{
	State sol;
	State s = rhs(euler1(r,0.5*h,par),par);
	for (int_fast32_t i=0;i<4;++i)
		sol.push_back(r[i]+h*s[i]);
	return sol;
}

State rk4(const State &r, const Real &h, const Parametry &par)
{
	State k1, k2, k3, k4, yk1, yk2, yk3, sol;

	State s = rhs(r,par);
	for (int_fast32_t i=0;i<4;++i)
		k1.push_back(h*s[i]); //obliczamy k1

	for (int_fast32_t i=0;i<4;++i)
		yk1.push_back(r[i]+0.5*k1[i]); //obliczamy yk1

	State sk1 = rhs(yk1,par);
	for (int_fast32_t i=0;i<4;++i)
		k2.push_back(h*sk1[i]); //obliczamy k2

	for (int_fast32_t i=0;i<4;++i)
		yk2.push_back(r[i]+0.5*k2[i]); //obliczamy yk2

	State sk2 = rhs(yk2,par);
	for (int_fast32_t i=0;i<4;++i)
		k3.push_back(h*sk2[i]); //obliczamy k3

	for (int_fast32_t i=0;i<4;++i)
		yk3.push_back(r[i]+k3[i]); //obliczamy yk3

	State sk3 = rhs(yk3,par);
	for (int_fast32_t i=0;i<4;++i)
		k4.push_back(h*sk3[i]); //obliczamy k3

	for (int_fast32_t i=0;i<4;++i)
		sol.push_back(r[i]+(k1[i]+2*k2[i]+2*k3[i]+k4[i])/6); //obliczamy k3

	return sol;
}

int main(int argc, const char *argv[])
{
	Real time;
	int_fast32_t steps, nplot;
	Parametry par;
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
		cerr << "time = ";
	cin >> time;
	if (verbose)
		cerr << "Liczba krokow = ";
	cin >> steps;
	Real h = time/steps;
	if (verbose)
	{
		cerr << "Dlugosc kroku wynosi" << '\t' << h << '\n';
		cerr << "nplot = ";
	}
	cin >> nplot;

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

	if (nplot==0)
		cout << setw(dist) << time << setw(dist) << solution[steps][0] << setw(dist) << solution[steps][1] << setw(dist) << 
		solution[steps][2] << setw(dist) << solution[steps][3] << '\n';
	else
	{
		for (int_fast32_t i=0;i<=steps;++i)
		{
			if (i%nplot==0)
				cout << setw(dist) << i*h << setw(dist) << solution[i][0] << setw(dist) << solution[i][1] << setw(dist) <<
				solution[i][2] << setw(dist) << solution[i][3] << '\n';
		}
	}

	return 0;
}
