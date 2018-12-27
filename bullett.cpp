#include <iostream>
#include <cinttypes>
#include <cmath>
#include <iomanip>
#include <vector>

using namespace std;

typedef float real;
typedef vector<real> state; //wektor stanu (x,y,vx,vy) 

const real g=9.81;
const real C=0.001;
const real mass=1;
const real n=2;

state RHS(const state &r)
{
	state deriv;
	deriv.push_back(r[2]);
	deriv.push_back(r[3]);
	deriv.push_back(-C*pow(r[2]*r[2]+r[3]*r[3], (n-1)/2)*r[2]/mass);
	deriv.push_back(-C*pow(r[2]*r[2]+r[3]*r[3], (n-1)/2)*r[3]/mass-g);
	return deriv;
}

state Euler1(const state &r, const real &h)
{
	state sol;
	for(int i=0;i<4;++i)
	{
		sol.push_back(r[i]+h*RHS(r)[i]);
	}
	return sol;
}

state Euler2(const state &r, const real &h)
{
	state sol;
	for(int i=0;i<4;++i)
	{
		sol.push_back(r[i]+h*RHS(Euler1(r,0.5*h))[i]);
	}
	return sol;
}

state RK4(const state &r, const real &h)
{

	state k1,k2,k3,k4,yk1,yk2,yk3,sol;
	
	for(int i=0;i<4;++i)
	{
		k1.push_back(h*RHS(r)[i]); //obliczamy k1
	}

	for(int i=0;i<4;++i)
	{
		yk1.push_back(r[i]+0.5*k1[i]); //obliczamy yk1
	}

	for(int i=0;i<4;++i)
	{
		k2.push_back(h*RHS(yk1)[i]); //obliczamy k2
	}

	for(int i=0;i<4;++i)
	{
		yk2.push_back(r[i]+0.5*k2[i]); //obliczamy yk2
	}
	
	for(int i=0;i<4;++i)
	{
		k3.push_back(h*RHS(yk2)[i]); //obliczamy k3
	}

	for(int i=0;i<4;++i)
	{
		yk3.push_back(r[i]+k3[i]); //obliczamy yk3
	}

	for(int i=0;i<4;++i)
	{
		k4.push_back(h*RHS(yk3)[i]); //obliczamy k3
	}


	for(int i=0;i<4;++i)
	{
		sol.push_back(r[i]+(k1[i]+2*k2[i]+2*k3[i]+k4[i])/6); //obliczamy k3
	}

	return sol;
}



int main()
{
	real time,h,hplot;
	 
	cout<<"time="<<'\n';
	cin>>time;
	cout<<"Dlugosc kroku h=?"<<'\n';
	cin>>h;
	cout<<"Dlugosc kroku do rysowania wykresu hplot=?"<<'\n';
	cin>>hplot;

	int steps=100;
	int kplot=h/hp;
	state r(4);
	r[0]=0.0;
	r[1]=0.0;
 

	cout<<"vx=?"<<'\n';
	cin>>r[2];
	cout<<"vy=?"<<'\n';
	cin>>r[3];

	vector<state> solution;
	solution.push_back(r);

	for(int i=1;i<steps;++i)
	{
		r=RK4(r,h);
		solution.push_back(r);
	}

	for(int i=0;i<steps;i+=kplot)
	{
		cout<<i*h<<'\t'<<solution[i][0]<<'\t'<<solution[i][1]<<'\t'<<solution[i][2]<<'\t'<<solution[i][3]<<'\n';
	}

}