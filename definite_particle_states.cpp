#include<iostream>
#include<iomanip>
#include<bitset>
#include<stdlib.h>

using namespace std; 
#define Nqubits 5

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

bitset<Nqubits>* definite_particle_states(int m)
{
 bitset<Nqubits> u;
 int LCm;
 int nchoosek(int,int);
 LCm = nchoosek(Nqubits,m);
 cout << " LCm = " << LCm << endl;
 bitset<Nqubits> *states = new bitset<Nqubits>[(const unsigned int)LCm];

 
 int i,j,k,l;
 for(k=0;k<m;k++)
    u.set(k);
 states[0] = u;
    
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
   states[l] = u; 
 }
   return states;
}

