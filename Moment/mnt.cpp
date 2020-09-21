#include<iostream>
#include<cmath>
#include<ctime>
#include<cstdlib>
#include<vector>
#include<fstream>
#include <boost/random.hpp>
#include <boost/random/beta_distribution.hpp>


using namespace std;
int main()
{  
ofstream fout ("mnt5.txt"); //final output file

unsigned seed=time(NULL);
srand(time(NULL));

int en;
cout<<"ensemble="<<endl;
cin>>en;		    //input ensemble number

int test_no;
cout<<"text_no="<<endl;
cin>>test_no;		    //input the number of iterations

double time_interval;
cout<<"time_interval="<<endl;
cin>>time_interval;	    //time interval,eachtime after which i want to record data

int starting_time;
cout<<"starting_time="<<endl;
cin>>starting_time;	    //the time from when i want to start recording data


vector<double>N(5000000);   //this vector will keep record of time
vector<double>delta(5000000);//this will have the sum of df th power of each length
int ro=0;

for(int f=0;f<en;f++)		//loop for ensemble average
    {
    vector<double>max(200000);
    vector<double>min(200000);
    vector<double>particle(200000);
    vector<double>gap(200000);
    vector<double> occupied_min;
    vector<double> occupied_max;

    max.insert(max.begin()+0,1);
    min.insert(min.begin()+0,0);
    occupied_max.push_back(1);
    occupied_min.push_back(0);
 
    double test=test_no,summation=1,ms=0.0;
    int so=0,element=0;
    int z=1;
    int seed = 2018;
    typedef boost::random::mt19937 RandomNumberGenerator;
    typedef boost::random::beta_distribution<> BetaDistribution;
    RandomNumberGenerator Rng(seed);
    BetaDistribution distribution(2,2);
    for (int i = 0; i <= test; i++)		//loop for time
        {
        int w=0,x=0;
        double random_num =(rand() / (double) RAND_MAX);
        double y=random_num;
        double e=0.0;
        for (int j = 0; j < occupied_min.size(); j++)	    //loop for selecting the gap that contain the generating random number and devide the according to rule, otherwise generating a new random number
            {
            if(y > min[j] && y < max[j])
		{
		double b1=distribution(Rng); 		    //generating random number using beta distribution
		double c=(occupied_max[j]-occupied_min[j])*b1;
		double d1=(occupied_max[j]-occupied_min[j])-c;
		
		double p1=rand()/(double)RAND_MAX;	    //a random number to decide whether the particle is difusing to left or to right
		if(p1<=0.5){ e=c+occupied_min[j];
                            //cout<<"left"<<endl;
                           }
                else if(p1>0.5){e=d1+occupied_min[j];
                                //cout<<"right"<<endl;
                               }

                double p=rand()/(double)RAND_MAX;	    //randomnumber that carry the probability of keepning the segmented partof the selected length

		if(p<=0.75){occupied_max.insert(occupied_max.begin() + j , e);
                            occupied_min.insert(occupied_min.begin() + j + 1, e);
			    w+=1;			   //w=1 means the segmented part will stay
			   }
		else if(p>0.75){if(occupied_max[j]==1){so+=1;}
				occupied_max[j]=e;
			        x+=1;        		   //x=1 means the segmented part eill be removed
                               }
		}
	    if(w==1 || x==1){break;}

	    else continue;
	    } 
	  
	double line=0.0;


        for(int check=0;check<occupied_max.size();check++)
	    {
	    particle[check]=occupied_max[check]-occupied_min[check];
	    if(check+1<occupied_max.size())
 	    {gap[check]=occupied_min[check+1]-occupied_max[check];}
	    }
    
	for(int s=0;s<occupied_max.size();s++)		//loop to generate a new line after resizing each length
	    {
	    max[s]=min[s+1]=min[s]+pow(particle[s],3);
            line+= pow(particle[s],3);
	    summation=line;
	    }
	double p = 0.0,avg=0.0;
	int num=0;
	if(i==z)
   	    {
    	    for(int k=0; k<occupied_max.size(); k++)
    		{ 
		p+=pow((occupied_max[k] - occupied_min[k]),.7787);
    		}
	    N[element]=i+1;	//record of time
	    delta[element]+=p;  //sum of df th power of each length
	    z+=time_interval;
	    element+=1;
	    }
	}
    if(f==ro)
	{
	cout<<"ensemble no......................................................"<<f<<endl;
	ro+=100;
	}
    }
for(int i=0;i<(N.size()-1);i++)
    {if(delta[i]!=0)
	{fout<<N[i]<<'\t'<<(delta[i]/en)<<endl;}
    }

return 0;
}
