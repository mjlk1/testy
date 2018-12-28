#include <iostream>
#include <cinttypes>
#include <cmath>
#include <iomanip>
#include <vector>

using namespace std;

typedef float Real;
typedef vector<Real> State; //wektor stanu (x,y,vx,vy) 

const Real g = 9.81;
const Real C = 0.001;
const Real mass = 1.0;
const Real n = 2.0;
const int_fast32_t  dist = 2*sizeof(Real)+3;

State rhs(const State &r)
{
	State deriv;
	Real l = -C*pow(r[2]*r[2]+r[3]*r[3], (n-1)/2)/mass;
	deriv.push_back(r[2]);
	deriv.push_back(r[3]);
	deriv.push_back(l*r[2]);
	deriv.push_back(l*r[3]-g);
	return deriv;
}

State euler1(const State &r, const Real &h)
{
	State sol;
	for (int_fast32_t i=0;i<4;++i)
		sol.push_back(r[i]+h*rhs(r)[i]);
	return sol;
}

State euler2(const State &r, const Real &h)
{
	State sol;
	State s = rhs(euler1(r,0.5*h));
	for (int_fast32_t i=0;i<4;++i)
		sol.push_back(r[i]+h*s[i]);
	return sol;
}

State rk4(const State &r, const Real &h)
{
	State k1, k2, k3, k4, yk1, yk2, yk3, sol;
	
	for (int_fast32_t i=0;i<4;++i)
		k1.push_back(h*rhs(r)[i]); //obliczamy k1

	for (int_fast32_t i=0;i<4;++i)
		yk1.push_back(r[i]+0.5*k1[i]); //obliczamy yk1

	for (int_fast32_t i=0;i<4;++i)
		k2.push_back(h*rhs(yk1)[i]); //obliczamy k2

	for (int_fast32_t i=0;i<4;++i)
		yk2.push_back(r[i]+0.5*k2[i]); //obliczamy yk2
	
	for (int_fast32_t i=0;i<4;++i)
		k3.push_back(h*rhs(yk2)[i]); //obliczamy k3

	for (int_fast32_t i=0;i<4;++i)
		yk3.push_back(r[i]+k3[i]); //obliczamy yk3

	for (int_fast32_t i=0;i<4;++i)
		k4.push_back(h*rhs(yk3)[i]); //obliczamy k3

	for (int_fast32_t i=0;i<4;++i)
		sol.push_back(r[i]+(k1[i]+2*k2[i]+2*k3[i]+k4[i])/6); //obliczamy k3

	return sol;
}


int main()
{
	Real time;
	int_fast32_t steps, nplot;
	 
	cout << "time = ";
	cin >> time;
	cout << "Liczba krokow = ";
	cin >> steps;
	Real h = time/steps;
	cout << "Dlugosc kroku wynosi" << '\t' << h << '\n';
	cout << "nplot = ";
	cin >> nplot;

	State r(4);
	r[0] = 0.0;
	r[1] = 0.0;
 
 	cout << "vx = ";
	cin >> r[2];
	cout << "vy = ";
	cin >> r[3];

	vector<State> solution;
	solution.push_back(r);

	for (int_fast32_t i=0;i<steps;++i)
	{
		r = rk4(r,h);
		solution.push_back(r);
	}

	for (int_fast32_t i=0;i<=steps;++i)
	{
		if (i%nplot==0)
			cout << setw(dist) << i*h << setw(dist) << solution[i][0] << setw(dist) << solution[i][1] << setw(dist) << solution[i][2] 
			<< setw(dist) << solution[i][3] << '\n';
	}

	return 0;
}