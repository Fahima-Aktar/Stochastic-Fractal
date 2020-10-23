#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <iostream>
#include<fstream>

using namespace std;

double gamm(double x) //Gamma function
{
    double ret = (1.000000000190015 + 
                 76.18009172947146 / (x + 1) +  
                 -86.50532032941677 / (x + 2) + 
                 24.01409824083091 / (x + 3) +  
                 -1.231739572450155 / (x + 4) + 
                 1.208650973866179e-3 / (x + 5) + 
                 -5.395239384953e-6 / (x + 6));
    
    return ret * sqrt(2*M_PI)/x * pow(x + 5.5, x+.5) * exp(-x-5.5);
}

double f(double n, double p, double alpha )//Fractal dimension function, p = probability, n = fractal dimension
{
return ((1+p)*gamm(n+alpha)/gamm(n+(2*alpha))-(gamm(alpha)/gamm(2*alpha)));
}

double BS(double a, double b, int N, double p, double alpha)//function to get fractal dimension by solving the dimension function using bisection method
    {double eps=0.0000001, c=0; 		//a,b are two arbitrary real numbers, N is number of iterations
    if(f(a,p,alpha)*f(b,p,alpha)>0){cout<<"Error!!"<<endl;return 0;}// if f(a)*f(b)>0 then it will show an error, f(a) and f(b) should have opposite sign
    if (f(a,p,alpha)==0){return a;}
    if (f(b,p,alpha)==0){return b;}

 
    for(int i=1; i<=N; i++)
	{if(abs(a-b)<eps){break;}
	c=(a+b)/2;
	if(f(c,p,alpha)==0){return c;}
	else if(f(c,p,alpha)*f(a,p,alpha)<0){b=c;}
	else {a=c;}
 	}
    return c;
    }


/*int main()
{ofstream fout ("df_p_.txt");//produce fractal dimension for different probability at a fixed value of alpha
double p=0.15, alpha=0.5;
for(int i=1; i<7; i++)
{fout<<p<<"    "<<BS(1,0,150,p,alpha)<<endl;
 cout<<p<<"    "<<BS(1,0,150,p,alpha)<<endl;
p+=0.15;}
return 0;
}*/ 

int main()
{ofstream fout ("df_alpha_.txt");//ppproduce fractal dimension as a function of alpha at a fixed probabilty
double p=0.15, alpha=1;
fout<<0.5<<"    "<<BS(1.0,0,150,p,0.5)<<endl;
for(int i=1; i<8; i++)
{fout<<alpha<<"    "<<BS(1.0,0,150,p,alpha)<<endl;
 cout<<alpha<<"    "<<BS(1,0,150,p,alpha)<<endl;
alpha+=1;}
return 0;
}

