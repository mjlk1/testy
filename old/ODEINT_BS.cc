#include <iostream>
#include <vector>
#include <fstream>

#include <boost/numeric/odeint.hpp>

/* alfy sa rozbite na czesc re i im, alfa_k=a_k+i*b_k. Podstawowa 'zmienna' to tablica tab=[a_{-nmax},...,a_{nmax},b_{-nmax},...,b_{nmax}].
 a_k=tab[nmax+k], b_k=tab[3*nmax+k] */
const int nmax=4;
const double gam=0.005;

typedef boost::array<double,4*nmax+2> state_type;

double ped(const state_type &x)
{
  double p=0;
  for(int i =-nmax;i<=nmax;++i)
  {
    p+=i*(x[nmax+i]*x[nmax+i]+x[3*nmax+1+i]*x[3*nmax+1+i]);
  }
  return p;
}


double particle_count(const state_type &x)
{
  double N=0;
  for(int i =-nmax;i<=nmax;++i)
  {
    N+=x[nmax+i]*x[nmax+i]+x[3*nmax+1+i]*x[3*nmax+1+i];
  }
  return N;
}


double kinetic_energy(const state_type &x)
{
  double ekin=0;
  for(int i =-nmax;i<=nmax;++i)
  {
    ekin+=i*i*(x[nmax+i]*x[nmax+i]+x[3*nmax+1+i]*x[3*nmax+1+i]);
  }
  return ekin;
}


double potential_energy(const state_type &x )
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
  return gam*epot;
}


/* The rhs of x' = f(x) */
void RHS( const state_type &x , state_type &dxdt , const double /* t */ )
{
  double sa,sb;
  for(int i=-nmax; i<=nmax;++i)
  {
    sa=0;
    sb=0;    
    for(int j = -nmax ; j <= nmax ; ++j)
    {
      for(int k = -nmax ; k <= nmax ; ++k)
      {
        for(int l = -nmax ; l <= nmax ; ++l)
        {
          if(k+l-j==i)
          {
            sa+=x[3*nmax+1+j]*x[3*nmax+1+k]*x[3*nmax+1+l] + x[nmax+j]*x[nmax+k]*x[3*nmax+1+l] + 
            x[nmax+j]*x[3*nmax+1+k]*x[nmax+l]-x[3*nmax+1+j]*x[nmax+k]*x[nmax+l];

            sb+=x[nmax+j]*x[nmax+k]*x[nmax+l] + x[3*nmax+1+j]*x[3*nmax+1+k]*x[nmax+l] + 
            x[3*nmax+1+j]*x[nmax+k]*x[3*nmax+1+l]-x[nmax+j]*x[3*nmax+1+k]*x[3*nmax+1+l];
          }
        }
      }
    dxdt[nmax+i]=(i*i*x[3*nmax+1+i]+2*gam*sa);
    dxdt[3*nmax+1+i]=-(i*i*x[nmax+i]+2*gam*sb);
    }
  }
}


//[ integrate_observer
struct push_back_state_and_time
{
    std::vector< state_type >& m_states;
    std::vector< double >& m_times;

    push_back_state_and_time( std::vector< state_type > &states , std::vector< double > &times )
    : m_states( states ) , m_times( times ) { }

    void operator()( const state_type &x , double t )
    {
        m_states.push_back( x );
        m_times.push_back( t );
    }
};
//]

struct write_state
{
    void operator()( const state_type &x ) const
    {
        std::cout << x[0] << "\t" << x[1] << "\n";
    }
};


int main()
{
  using namespace std;
  using namespace boost::numeric::odeint;
  //[ state_initialization
  state_type x;
   
  x[0]=0.7780577053810895;
  x[1]=-1.6332089457627652;
  x[2]=-2.1542143774721465;
  x[3]=2.203353598100645;
  x[4]=2.459659063594901;
  x[5]=-1.2893247047018357;
  x[6]=-1.2187063853406443;
  x[7]=0.35132876928612744;
  x[8]=-0.9795038767280292;
  x[9]=1.0695793730674183;
  x[10]=0.05452654682080023;
  x[11]= 1.9462319461908129;
  x[12]=2.803849370527668;
  x[13]= 6.3575917812193214;
  x[14]=3.7355855316576934;
  x[15]= -2.90943171677841;
  x[16]=-0.8621102747837353;
  x[17]= 0.7572974029938235;

  vector<state_type> x_vec;  //Tutaj sa przechowywane punkty rozwiazania
  vector<double> times;
  
  bulirsch_stoer_dense_out< state_type > stepper( 1E-10 , 1E-10 , 1.0 , 1.0 ); //metoda calkowania

  size_t steps=integrate_adaptive(stepper, RHS, x , 0.0 , 100.0 , 0.1, push_back_state_and_time( x_vec , times ));

  cout<<setprecision(9)<<kinetic_energy(x_vec[0])+potential_energy(x_vec[0])<<endl;  //sprawdzenie o ile sie energia zmienila
  cout<<setprecision(9)<<kinetic_energy(x_vec[steps-1])+potential_energy(x_vec[steps-1])<<endl;

}