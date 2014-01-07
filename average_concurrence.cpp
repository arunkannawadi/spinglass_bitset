#include<iostream>
#include<iomanip>
#include<bitset>
#include<stdlib.h>
#include<math.h>
#include<Eigen/Dense>
#include<Eigen/Eigenvalues>

using namespace std; 
using namespace Eigen;

#define root2 1.414213562373095
#define Nqubits 4

int main()
{
  struct timeval start;
  struct timeval end;
  double t1,t2,t3 ,t4;

  time_t seconds;
  time(&seconds);
  srand((unsigned int)seconds);
  
  VectorXd random_definite_particle_states(int);
  double* average_concurrence(VectorXd,int);
  double *entanglement;
  
  VectorXd state;
  
  int n,m=1,N=1;
  for(n=0;n<N;n++)
  {
   state = random_definite_particle_states(m);
   entanglement = average_concurrence(state,m);
   cout << "Avg C = " << entanglement[0] << endl;
   cout << "Avg Tangle = " << entanglement[1] << endl;
   cout << "State = " << state << endl;
  }
}

double pr(VectorXd state)
{
 double IPR, PR;
 int k,n;
 
 n = state.size();
 for(IPR=0,k=0;k<n;k++)
  IPR += pow(state(k),4);
 PR = 1.0/IPR;
 return PR;
}

double* average_concurrence(VectorXd state, int m)
{
  int i,j,k,l,L=Nqubits;
  int nchoosek(int,int);
  int LCm; LCm = nchoosek(L,m);
  
  bitset<Nqubits>* definite_particle_basis(int);
  bitset<Nqubits> *basis;
  basis = definite_particle_basis(m);
  
  double v,y,z;
  double C, avgC, avgTangle, sumTangle=0,sumC=0; 
  for(i=1;i<L;i++)
   for(j=0;j<i;j++)
   {
    v=0; y=0; z=0;
    for(k=0;k<LCm;k++)
     {
      if((basis[k][i]==1)&(basis[k][j]==1))
        y += state(k)*state(k);
      else if((basis[k][i]==0)&(basis[k][j]==0))
        v += state(k)*state(k);
      else
      {
        for(l=0;l<k;l++)
         {
          basis[l].flip(i); basis[l].flip(j);
          if(basis[k]==basis[l])
           z += state(k)*state(l);
          basis[l].flip(i); basis[l].flip(j);
         }
      }  
      
     }
     C = ((abs(z)-sqrt(v*y))>0) ? 2.0*(abs(z)-sqrt(v*y)) : 0.0 ;
     cout << "C = " << C << endl;
     sumC += C;
     sumTangle += C*C;
   }
   
    int LC2; LC2 = nchoosek(L,2);
    
    avgC = sumC/LC2;
    avgTangle = sumTangle/LC2;
    
    double *entanglement = new double[2];
    entanglement[0] = avgC;
    entanglement[1] = avgTangle;
    
    return entanglement;
}

VectorXd random_definite_particle_states(int m)
{
  int k,LCm, L = Nqubits;

  double u,v,s,x,y,z;
  int nchoosek(int,int);
  LCm = nchoosek(L,m);
  
  VectorXd state = VectorXd::Zero(LCm);
  
  for(k=0;k<LCm;k++)
  {
    do
    {
     u = 2*((double)rand())/RAND_MAX -1;
     v = 2*((double)rand())/RAND_MAX -1;
     s = u*u + v*v ;
     x = u*sqrt(-2*log(s)/s);
     y = v*sqrt(-2*log(s)/s);
     z = (x+y)/root2;  
     cout << z << endl; 
    } while(s>=1);
    
    state(k) = z;
  }
  
  return state/state.norm();
}

bitset<Nqubits>* definite_particle_basis(int m)
{
 bitset<Nqubits> u;
 int LCm;
 int nchoosek(int,int);
 LCm = nchoosek(Nqubits,m);
 bitset<Nqubits> *basis = new bitset<Nqubits>[(const unsigned int)LCm];
 
 int i,j,k,l;
 for(k=0;k<m;k++)
    u.set(k);
 basis[0] = u;
    
 for(l=1;l<LCm;l++)
 {
   for(i=0;u[i]==0;i++);
   for(j=i;u[j]==1;j++);
   
   u.set(j);
   u.reset(j-1);
   for(k=0;k<j-i-1;k++)
   {
     u.set(j-k-2,0);
   }
   for(k=0;k<j-i-1;k++)
   {  
     u.set(k,1);
   }
   basis[l] = u; 
 }
 
 return basis;
}

int nchoosek(int n, int k)
{
 int r;
 float nck;
 if(2*k>n)
   k = n-k;
   
 for(nck=1.0, r=0;r<k;r++)
 {
  nck = nck*(n-r)/(r+1) ; 
 }
 return (int)nck;
}
