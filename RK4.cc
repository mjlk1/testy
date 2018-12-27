#include <iostream>
#include <cstdlib>  
#include <fstream>
#include <iomanip>
#include <vector>
#include <boost/array.hpp>
using namespace std;

/* alfy sa rozbite na czesc re i im, alfa_k=a_k+i*b_k. Podstawowa 'zmienna' to tablica tab=[a_{-nmax},...,a_{nmax},b_{-nmax},...,b_{nmax}].
 a_k=tab[nmax+k], b_k=tab[3*nmax+k] */

const int nmax=1;
const double gamma=0.005;

typedef boost::array<double,4*nmax+2> state; //wektor stanu

double ped(const state &x)
{
  double p=0;
  for(int i =-nmax;i<=nmax;++i)
  {
    p+=i*(x[nmax+i]*x[nmax+i]+x[3*nmax+1+i]*x[3*nmax+1+i]);
  }
  return p;
}


double particle_count(const state &x)
{
  double N=0;
  for(int i =-nmax;i<=nmax;++i)
  {
    N+=x[nmax+i]*x[nmax+i]+x[3*nmax+1+i]*x[3*nmax+1+i];
  }
  return N;
}


double kinetic_energy(const state &x)
{
  double ekin=0;
  for(int i =-nmax;i<=nmax;++i)
  {
    ekin+=i*i*(x[nmax+i]*x[nmax+i]+x[3*nmax+1+i]*x[3*nmax+1+i]);
  }
  return ekin;
}


double potential_energy(const state &x )
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
            epot+= x[nmax+i]*x[nmax+j]*x[nmax+k]*x[nmax+l] + x[nmax+i]*x[3*nmax+1+j]*x[3*nmax+1+k]*x[nmax+l] +
            x[nmax+i]*x[3*nmax+1+j]*x[nmax+k]*x[3*nmax+1+l] - x[nmax+i]*x[nmax+j]*x[3*nmax+1+k]*x[3*nmax+1+l] + 
            x[3*nmax+1+i]*x[3*nmax+1+j]*x[3*nmax+1+k]*x[3*nmax+1+l] + x[3*nmax+1+i]*x[nmax+j]*x[nmax+k]*x[3*nmax+1+l] + 
            x[3*nmax+1+i]*x[nmax+j]*x[3*nmax+1+k]*x[nmax+l] - x[3*nmax+1+i]*x[3*nmax+1+j]*x[nmax+k]*x[nmax+l];
          }
        }
      }
    }
  }
  return gamma*epot;
}




void RHS(const state &x1, state  &x2, double h) //funkcja zadajaca rnanie rozniczkowe
{
  double sa,sb;
  for(int i=-nmax ; i <= nmax;++i)
  {
    sa=0;
    sb=0;   
    for(int j = -nmax ; j <= nmax ; ++j)
    {
      for(int k = -nmax ; k <= nmax ; ++k)
      {
        for(int l = -nmax ; l <= nmax ; ++l)
        {
          if( k + l - j == i)
          {
            sa+=x1[3*nmax+1+j]*x1[3*nmax+1+k]*x1[3*nmax+1+l] + x1[nmax+j]*x1[nmax+k]*x1[3*nmax+1+l] + 
            x1[nmax+j]*x1[3*nmax+1+k]*x1[nmax+l]-x1[3*nmax+1+j]*x1[nmax+k]*x1[nmax+l];

            sb+=x1[nmax+j]*x1[nmax+k]*x1[nmax+l] + x1[3*nmax+1+j]*x1[3*nmax+1+k]*x1[nmax+l] + 
            x1[3*nmax+1+j]*x1[nmax+k]*x1[3*nmax+1+l]-x1[nmax+j]*x1[3*nmax+1+k]*x1[3*nmax+1+l];
          }
        }
      }
    }
    x2[nmax+i]=h*(i*i*x1[3*nmax+1+i]+2*gamma*sa);
    x2[3*nmax+1+i]=-h*(i*i*x1[nmax+i]+2*gamma*sb);
  }
}
int main()
{
  int steps=10000;
  double time=1000;
  double h=time/((double) steps);

  vector <state> solution; //tutaj sa przechowywane punkty obliczone w klejnych krokach

  state x; /*definicja wektora stanu */
  
  x[0]=-1.225822559260336;
  x[1]=9.769198838765893;
  x[2]=1.2448968731974293;
  x[3]=-0.23699325268771354;
  x[4]=-1.2054840942332938;
  x[5]=0.031419320455167786;
  
  solution.push_back(x);
  
  state k1,k2,k3,k4,yk1,yk2,yk3;
  /*Powyzej zdefiniowalem tablice ktore odpowiadaja  wektorom k1,k2,k3,k4 z RK4, a takze tablice ktore odpowiadaja wektorom y_n+k1/2,
  y_n+k2/2, y_n+k3. */
  for(int i=1;i<steps;++i)
  {
    RHS(x,k1,h);
    for(int j=0;j<4*nmax+2;++j)
    {
      yk1[j]=x[j]+k1[j]/2;
    }
      
    RHS(yk1,k2,h);
    for(int j=0;j<4*nmax+2;++j)
    {
      yk2[j]=x[j]+k2[j]/2;
    }
      
    RHS(yk2,k3,h);
    for(int j=0;j<4*nmax+2;++j)
    {
      yk3[j]=x[j]+k3[j];
    }
      
    RHS(yk3,k4,h);   
    for(int j=0;j<4*nmax+2;++j)
    {
      x[j]=x[j]+(k1[j]+2*k2[j]+2*k3[j]+k4[j])/6;
    }
    solution.push_back(x);

  
  }
 
 
  double nzero[steps]; 
  for(int i=0;i<steps;++i)
  {
    nzero[i]=solution[i][nmax]*solution[i][nmax]+solution[i][3*nmax+1]*solution[i][3*nmax+1];
  }

  
  cout<<setprecision(9)<<kinetic_energy(solution[0])+potential_energy(solution[0])<<endl;  //sprawdzenie o ile sie energia zmienila
  cout<<setprecision(9)<<kinetic_energy(solution[steps-1])+potential_energy(solution[steps-1])<<endl;



  /*zapis do pliku danych n_0(t)*/
  /*fstream file;
  file.open("dane.txt", ios::out);
  
  for(int i=0;i<steps;++i)
  {
    file<<time*((double)(i))/((double)steps)<<" "<<nzero[i]<<endl;
  }
  file.close(); */
  
}

