#include <iostream>
#include <cstdlib>  
#include <fstream>
#include <iomanip>
using namespace std;

/* alfy sa rozbite na czesc re i im, alfa_k=a_k+i*b_k. Podstawowa 'zmienna' to tablica tab=[a_{-nmax},...,a_{nmax},b_{-nmax},...,b_{nmax}]. a_k=tab[nmax+k], 
b_k=tab[3*nmax+k] */

template <typename T>
double ped(T* tin, int nmax)
{
	double p=0;
	for(int i =-nmax;i<=nmax;++i)
	{
	p+=i*(tin[nmax+i]*tin[nmax+i]+tin[3*nmax+1+i]*tin[3*nmax+1+i]);
 	}
 	return p;
}

template <typename T>
double nall(T* tin, int nmax)
{
 	double N=0;
 	for(int i =-nmax;i<=nmax;++i)
 	{
 	N+=tin[nmax+i]*tin[nmax+i]+tin[3*nmax+1+i]*tin[3*nmax+1+i];
	}
	return N;
}

template <typename T>
double ekin(T* tin, int nmax)
{
	double ekin=0;
	for(int i =-nmax;i<=nmax;++i)
	{
	ekin+=i*i*(tin[nmax+i]*tin[nmax+i]+tin[3*nmax+1+i]*tin[3*nmax+1+i]);
 	}
 	return ekin;
}

template <typename T>
double epot(T* tin, int nmax, double gamma)
{
	double epot=0;
	for(int i = -nmax ; i <= nmax ; ++i)
	{ 
		for(int j = -nmax ; j <= nmax ; ++j)
    	{
    		for(int k = -nmax ; k <= nmax ; ++k)
        	{
         		for(int l = -nmax ; l <= nmax ; ++l)
            	{
            		if(i+j==k+l)
            		{
					epot+= tin[nmax+i]*tin[nmax+j]*tin[nmax+k]*tin[nmax+l] + tin[nmax+i]*tin[3*nmax+1+j]*tin[3*nmax+1+k]*tin[nmax+l] +
					tin[nmax+i]*tin[3*nmax+1+j]*tin[nmax+k]*tin[3*nmax+1+l] - tin[nmax+i]*tin[nmax+j]*tin[3*nmax+1+k]*tin[3*nmax+1+l] + 
					tin[3*nmax+1+i]*tin[3*nmax+1+j]*tin[3*nmax+1+k]*tin[3*nmax+1+l] + tin[3*nmax+1+i]*tin[nmax+j]*tin[nmax+k]*tin[3*nmax+1+l] + 
					tin[3*nmax+1+i]*tin[nmax+j]*tin[3*nmax+1+k]*tin[nmax+l] - tin[3*nmax+1+i]*tin[3*nmax+1+j]*tin[nmax+k]*tin[nmax+l];
		        	}
            	}
			}
		}
	}
	return gamma*epot;
}


/* To jest cos ala funkcja F zadajaca r nie rozniczkowe. Nie wiedzialem co zrobic  zeby funkcja zwracala tablice, wiec to cos bierze tablice z a i b (t1) i 
jakas druga tablice (t2) i przypisuje wartosci funkcji F do tablicy t2. */
template <typename T>
void F(T* t1,T* t2, int nmax, double gamma, double h)
{
	double sa,sb;
	for(int i=-nmax; i<=nmax;++i)
  	{
  		sa=0;
  		for(int j = -nmax ; j <= nmax ; ++j)
    	{
    		for(int k = -nmax ; k <= nmax ; ++k)
        	{
         		for(int l = -nmax ; l <= nmax ; ++l)
            	{
            		if(k+l-j==i)
            		sa+=t1[3*nmax+1+j]*t1[3*nmax+1+k]*t1[3*nmax+1+l] + t1[nmax+j]*t1[nmax+k]*t1[3*nmax+1+l] + 
            		t1[nmax+j]*t1[3*nmax+1+k]*t1[nmax+l]-t1[3*nmax+1+j]*t1[nmax+k]*t1[nmax+l];
            	}
        	}
        }
    	sb=0;    
    	for(int j = -nmax ; j <= nmax ; ++j)
    	{
    		for(int k = -nmax ; k <= nmax ; ++k)
        	{
         		for(int l = -nmax ; l <= nmax ; ++l)
            	{
            		if(k+l-j==i)
            		sb+=t1[nmax+j]*t1[nmax+k]*t1[nmax+l] + t1[3*nmax+1+j]*t1[3*nmax+1+k]*t1[nmax+l] + 
            		t1[3*nmax+1+j]*t1[nmax+k]*t1[3*nmax+1+l]-t1[nmax+j]*t1[3*nmax+1+k]*t1[3*nmax+1+l];
            	}
        	}
        }
  		t2[nmax+i]=h*(i*i*t1[3*nmax+1+i]+2*gamma*sa);
  		t2[3*nmax+1+i]=-h*(i*i*t1[nmax+i]+2*gamma*sb);
  	}
}
int main()
{
	int nmax = 1;
	long int steps=10000;
	double gamma=0.005;
	double time=100;
	double h=time/((double) steps);
  
	double sol[steps][4*nmax+2]; /*W tej tablicy jest przechowywana zaleznosc a_k i b_k od kroku */
  
	sol[0][0]=-1.225822559260336;
	sol[0][1]=9.769198838765893;
	sol[0][2]=1.2448968731974293;
	sol[0][3]=-0.23699325268771354;
	sol[0][4]=-1.2054840942332938;
	sol[0][5]=0.031419320455167786;
  
  	double tk1[4*nmax+2],tk2[4*nmax+2],tk3[4*nmax+2],tk4[4*nmax+2], tyk1[4*nmax+2], tyk2[4*nmax+2], tyk3[4*nmax+2];
	/*Powyzej zdefiniowalem tablice ktore odpowiadaja  wektorom k1,k2,k3,k4 z RK4, a takze tablice ktore odpowiadaja wektorom y_n+k1/2, y_n+k2/2, y_n+k3.
	Tablice te w kazdym kroku sa nadpisywane */
	for(int i=1;i<steps;++i)
	{
		F(sol[i-1], tk1,nmax,gamma,h);
		for(int j=0;j<=4*nmax+2;++j)
		{
			tyk1[j]=sol[i-1][j]+tk1[j]/2;
	    }
	    
	    F(tyk1, tk2,nmax,gamma,h);
	    for(int j=0;j<=4*nmax+2;++j)
		{
			tyk2[j]=sol[i-1][j]+tk2[j]/2;
	    }
	    
	    F(tyk2, tk3,nmax,gamma,h);
		for(int j=0;j<=4*nmax+2;++j)
		{
			tyk3[j]=sol[i-1][j]+tk3[j];
	    }
	    
		F(tyk3, tk4,nmax,gamma,h);	 
		for(int j=0;j<=4*nmax+2;++j)
		{
			sol[i][j]=sol[i-1][j]+(tk1[j]+2*tk2[j]+2*tk3[j]+tk4[j])/6;
	    }
	
	}
 
 
	double nzero[steps]; 
	for(int i=0;i<steps;++i)
	{
		nzero[i]=sol[i][nmax]*sol[i][nmax]+sol[i][3*nmax+1]*sol[i][3*nmax+1];
	}

	cout<<setprecision(9)<<ekin(sol[0],nmax)+epot(sol[0],nmax,gamma)<<endl;
	cout<<setprecision(9)<<ekin(sol[steps-3],nmax)+epot(sol[steps-3],nmax,gamma)<<endl;

	/*zapis do pliku danych n_0(t)*/
	fstream file;
	file.open("dane.txt", ios::out);
	
	for(int i=0;i<steps;++i)
	{
		file<<time*((double)(i))/((double)steps)<<" "<<nzero[i]<<endl;
	}
	file.close();
	
}
