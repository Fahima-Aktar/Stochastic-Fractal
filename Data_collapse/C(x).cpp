#include<iostream>
#include<cmath>
#include<ctime>
#include<cstdlib>
#include<vector>
#include<fstream>
#include <boost/random.hpp>
#include <boost/random/beta_distribution.hpp>

using namespace std;
int main(int argc, char *argv[])
{  
ofstream fout ("b13le20k0005.txt");	//final output file

unsigned seed=time(NULL);
srand(time(NULL));

vector<double>c(2000000);		//number of particle per unit length
vector<double>range(2000000);		//width of the histogram

   
int en=20000;		   //number of ensemble
int test_no=300000;	   //number of iteration
double x_interval=.0005;
int ro=0;

for(int f=0;f<en;f++)	  //loop for ensemble average 
    {
    vector<double>max(200000);
    vector<double>min(200000);
    vector<double>particle(200000);
    vector<double>gap(200000);
    vector<double>c1;
    vector<double> occupied_min;
    vector<double> occupied_max;

    max.insert(max.begin()+0,1);
    min.insert(min.begin()+0,0);
    occupied_max.push_back(1);
    occupied_min.push_back(0);
 
    double test=test_no,summation=1,ms=0.0,r=x_interval;
    int so=0;
    int z=1;
    int seed = 2018;
    typedef boost::random::mt19937 RandomNumberGenerator;
    typedef boost::random::beta_distribution<> BetaDistribution;
    RandomNumberGenerator Rng(seed);
    BetaDistribution distribution(2,2);

    for (int idx = 0 ; idx <= test ; ++idx) // generating random number using beta distribution 
        {  
	double rn=distribution(Rng);
        c1.push_back(rn); 
        }
 
    int nth=0;
    for (int i = 0; i <= test; i++)	    //time loop
        { 
        int w=0,x=0;
        double random_num =(rand() / (double) RAND_MAX);
        double y=random_num;
        double e=0.0;

        for (int j = 0; j < occupied_min.size(); j++)
            {
        	if(y > min[j] && y < max[j])
		    { 
		    double b1=c1[nth];
                    nth+=1;
		    double c=(occupied_max[j]-occupied_min[j])*b1;
		    double d1=(occupied_max[j]-occupied_min[j])-c;

                    double p1=rand()/(double)RAND_MAX;		//random number to decide whether the particle is difusing to left or to right
		    
		    if(p1<=0.5){ e=c+occupied_min[j];
                               //cout<<"left"<<endl;
                               }
                    else if(p1>0.5){e=d1+occupied_min[j];
                                   //cout<<"right"<<endl;
                                  }

                   double p=rand()/(double)RAND_MAX;		//randomnumber that carry the probability of keepning the segmented partof the selected length

		   if(p<=0.75){ occupied_max.insert(occupied_max.begin() + j , e);
                              occupied_min.insert(occupied_min.begin() + j + 1, e);
			      w+=1;			//w=1 means the segmented part will stay
		              }
		   else if(p>0.75){if(occupied_max[j]==1){so+=1;}
				  occupied_max[j]=e;
			          x+=1;			//x=1 means the segmented part eill be removed
				  }
		   }

            if(w==1 || x==1){break;}

	    else continue;
            }
	
	double line=0.0;

        for(int check=0;check<occupied_max.size();check++) 	
	    {
	    particle[check]=occupied_max[check]-occupied_min[check];
	    }
	
	for(int s=0;s<occupied_max.size();s++)		//loop to generate a new line after resizing each length
	    {
	    max[s]=min[s+1]=min[s]+pow(particle[s],3);
            line+= pow(particle[s],3);
	    summation=line;
	    }
	}

    double index=0.0;
    for (double d = 0; d <= 1; d += r)
   	{ 
	double q = 0;
        for (int i = 0; i < occupied_max.size(); i++) 
	    {
            double p = (occupied_max[i] - occupied_min[i]);	//particle length
            if (p >= d && p < (d + r))
	       {q++;}		//counting number of particle within specific length
            }
	c[index]+=q/r;
	range[index]=d+r/2;
        index+=1;
        }

    if(f==ro)
	{
	cout<<"ensemble no......................................................."<<f<<endl;
	ro+=1500;
	}
     }
    
for(int i=0;i<(c.size()-1);i++)
    {if(c[i] !=0)
   	{fout<<(range[i])<<'\t'<<(c[i]/en)<<endl;
        cout<<(range[i])<<'\t'<<(c[i]/en)<<endl;
        }
    }

return 0;
}
