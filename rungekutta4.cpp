#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;




double F(double x, double y)
{
	double z=2*sin(x)-y;
	return z;
}

int main()
{
	int n;
	cout<<" n"<<endl;
	cin>>n;

double Y[n+1];
double w;
	cout<<"war pocz"<<endl;
	cin>>w;
	double range;
	cout<<"range"<<endl;
	cin>>range;
	double h=range/((double)n);
	Y[0]=w;
	double k1,k2,k3,k4;
	for(int i=1; i<=n;++i)
	{
		k1=h*F(range*((double)(i-1))/((double)n),Y[i-1]);
		k2=h*F(range*((double)(i-1))/((double)n)+h/2,Y[i-1]+k1/2);
		k3=h*F(range*((double)(i-1))/((double)n)+h/2,Y[i-1]+k2/2);
		k4=h*F(range*((double)(i-1))/((double)n)+h,Y[i-1]+k3);
		Y[i]=Y[i-1]+(k1+2*k2+2*k3+k4)/6;
	
	}
	fstream file;
	file.open("dane.txt", ios::out);
	
	for(int i=0;i<n;++i)
	{
		file<<range*((double)(i))/((double)n)<<" "<<Y[i]<<endl;
	}
	file.close();
	


}
