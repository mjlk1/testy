#include <cinttypes>
#include <cmath>

#include "Parameters.hpp"
#include "Derivative.hpp"

State cashkarp0(const State &r, Real &h, Real &h_done, const Real &eps, const Parameters &par)
{
	State k1(4), k2(4), k3(4), k4(4), k5(4), k6(4), yk1(4), yk2(4), yk3(4), yk4(4), yk5(4), sol(4), error(4), error_diff(4), delta_0(4);
	State s, sk1, sk2, sk3, sk4, sk5;
	int_fast32_t worst_coord;

	Real b[7][6];
	Real c[7];
	Real c2[7];

	// b i c to wspolczynniki z RK5
	b[2][1] = 0.2;

	b[3][1] = 3.0/40.0;
	b[3][2] = 9.0/40.0;

	b[4][1] = 0.3;
	b[4][2] = -0.9;
	b[4][3] = 6.0/5.0;

	b[5][1] = -11.0/54.0;
	b[5][2] = 2.5;
	b[5][3] = -70.0/27.0;
	b[5][4] = 35.0/27.0;

	b[6][1] = 1631.0/55296.0;
	b[6][2] = 175.0/512.0;
	b[6][3] = 575.0/13824.0;
	b[6][4] = 44275.0/110592.0;
	b[6][5] = 253.0/4096.0;	

	c[1] = 37.0/378.0;
	c[2] = 0.0;
	c[3] = 250.0/621.0;
	c[4] = 125.0/594.0;
	c[5] = 0.0;
	c[6] = 512.0/1771.0;

	//c2 to wspolczynniki oznaczane jako c_i-c*_i. Sluza do obliczania bledu.
	c2[1] = c[1]-2825.0/27648.0;			
	c2[2] = c[2]-0.0;
	c2[3] = c[3]-18575.0/48384.0;  
	c2[4] = c[4]-13525.0/55296.0;
	c2[5] = c[5]-277.0/14336.0;
	c2[6] = c[6]-0.25;

	
	bool h_test = false; 
	
	while (!h_test)
	{
		h_test = true;

		s = derivative(r,par);
		for (int_fast32_t i=0;i<4;++i)
			k1[i] = h*s[i]; //obliczamy k1

		for (int_fast32_t i=0;i<4;++i)
			yk1[i] = r[i]+b[2][1]*k1[i]; //obliczamy yk1

		sk1 = derivative(yk1,par);
		for (int_fast32_t i=0;i<4;++i)
			k2[i] = h*sk1[i]; //obliczamy k2

		for (int_fast32_t i=0;i<4;++i)
			yk2[i] = r[i]+b[3][1]*k1[i]+b[3][2]*k2[i]; //obliczamy yk2

		sk2 = derivative(yk2,par);
		for (int_fast32_t i=0;i<4;++i)
			k3[i] = h*sk2[i]; //obliczamy k3

		for (int_fast32_t i=0;i<4;++i)
			yk3[i] = r[i]+b[4][1]*k1[i]+b[4][2]*k2[i]+b[4][3]*k3[i]; //obliczamy yk3

		sk3 = derivative(yk3,par);
		for (int_fast32_t i=0;i<4;++i)
			k4[i] = h*sk3[i]; //obliczamy k4

		for (int_fast32_t i=0;i<4;++i)
			yk4[i] = r[i]+b[5][1]*k1[i]+b[5][2]*k2[i]+b[5][3]*k3[i]+b[5][4]*k4[i]; //obliczamy yk4

		sk4 = derivative(yk4,par);
		for (int_fast32_t i=0;i<4;++i)
			k5[i] = h*sk4[i]; //obliczamy k5

		for (int_fast32_t i=0;i<4;++i)
			yk5[i] = r[i]+b[6][1]*k1[i]+b[6][2]*k2[i]+b[6][3]*k3[i]+b[6][4]*k4[i]+b[6][5]*k5[i]; //obliczamy yk5

		sk5 = derivative(yk5,par);
		for (int_fast32_t i=0;i<4;++i)
			k6[i] = h*sk5[i]; //obliczamy k6

		for (int_fast32_t i=0;i<4;++i)
			error[i] = c2[1]*k1[i]+c2[3]*k3[i]+c2[4]*k4[i]+c2[5]*k5[i]+c2[6]*k6[i]; //obliczamy blad

		for (int_fast32_t i=0;i<4;++i)
			sol[i] = r[i]+c[1]*k1[i]+c[3]*k3[i]+c[4]*k4[i]+c[6]*k6[i]; //obliczamy wynik

		for (int_fast32_t i=0;i<4;++i)
			delta_0[i]=eps*(abs(sol[i])+abs(k1[i])); //delta_0 zadaje maksymalny dopuszczalny blad dla kazdej zmiennej
		
		for (int_fast32_t i=0;i<4;++i)
		{
			if (abs(delta_0[i])<abs(error[i]))  // tutaj sprawdzamy, czy znajdzie sie jakas zmienna, ktorej blad jest wiekszy niz chcemy
				h_test = false;
		
			error_diff[i] = abs(error[i])-abs(delta_0[i]); // obliczamy roznice miedzy maks dopuszczalnym bledem a tym co jest teraz
		}

		worst_coord= 0;
		for (int_fast32_t i=0;i<4;++i)
		{	
			if (error_diff[i]>error_diff[worst_coord]) // szukamy zmiennej, ktora obliczylismy z najwiekszym bledem
				worst_coord = i;
		}

		if (!h_test)
			h = 0.95*h*pow(abs(delta_0[worst_coord]/error[worst_coord]), 0.25); // zmniejszamy h, 0.95 to wspolczynnik S
		else
		{
			h_done = h; // taki krok zrobilismy
			h = 0.95*h*pow(abs(delta_0[worst_coord]/error[worst_coord]), 0.2); // zwiekszamy h
		}
	}
	return sol;
}
	
