#include<iostream>
#include<iomanip>
#include<stdlib.h>
#include<math.h>

using namespace std; 

#define root2 1.414213562373095

int main()
{
  struct timeval start;
  struct timeval end;
  double t1,t2,t3 ,t4;

  time_t seconds;
  time(&seconds);
  srand((unsigned int)seconds);
  double **symm_gauss_rand_mat_gen(int);
  double **hamiltonian(double**, int,int);
  void printmat(double**,int); 
  
  double **J, **H;
  
  J = symm_gauss_rand_mat_gen(4);
  
  H = hamiltonian(J, 2, 4);
  
  printmat(J,4);
  cout << "\n \n" ; 
  printmat(H,6);
  
}  

double **symm_gauss_rand_mat_gen(int dim)
{
  int i,j;
  double** randmat;
  double u,v,s,x,y,z;

  randmat = (double**)malloc(dim*sizeof(double*));
  for(i=0;i<dim;i++)
    randmat[i] = (double*)malloc(dim*sizeof(double));

  for(i=0;i<dim;i++)
    {
      for(j=0;j<i;j++)
	{
	  do
	    {
	      u = 2*((double)rand())/RAND_MAX -1;
	      v = 2*((double)rand())/RAND_MAX -1;
	      s = u*u + v*v ;
	      x = u*sqrt(-2*log(s)/s);
	      y = v*sqrt(-2*log(s)/s);
	      z = (x+y)/root2;   
	    } while(s>=1);      
	  
	  randmat[i][j] = z;
	  randmat[j][i] = z;
	}
      randmat[i][i] = 0;
    }

  return randmat;
  
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

unsigned int *definite_particle_states(const unsigned int L, int m)
{
 int i,k, LCm;
 unsigned int v,w,t;
 unsigned int* states; 
 int nchoosek(int,int);
 
 LCm = nchoosek(L,m);
 states = (unsigned int*)malloc(LCm*sizeof(unsigned int));
 
 
 v = (1 << m) - 1;
 k = 0;
 
 while(k<LCm)
 {
  states[k++] = v;
  t = v | (v-1);
  v = (t+1) | (((~t & -~t) -1) >> (__builtin_ctz(v)+1));
 }
 
 return states;
 
}

double **hamiltonian(double** J, int m, int L)
{
  unsigned int *definite_particle_states(const unsigned int,int);
  int nchoosek(int,int);
  
  const unsigned int constL = L;  
  int LCm; 
  int i,j,k,n1,n2;
  unsigned int u,v,w,x,y;
    
  double tmp, sigmaJ = 0;
  for(i=0;i<L;i++)
     for(j=0;j<i;j++)
       sigmaJ += J[i][j];
    
  LCm = nchoosek(L,m);
  
  unsigned int *states;
  states =  definite_particle_states(constL,m);
  
  double **H;
  H = (double**)malloc(LCm*sizeof(double*));
  for(k=0;k<LCm;k++)
    {
      H[k] = (double*)malloc(LCm*sizeof(double));
      if(H[k]==NULL)
	{
         cout << "\n error allocating H2 " << k ;
	 break;
	}
    }

   for(i=0;i<LCm;i++)
     {
      for(j=0;j<i;j++)
      {
       u = states[i];
       v = states[j];       
       w = u^v; 
       
       x = w&(w-1);
       y = x&(x-1);
       
              
       if(y==0)
       {
       
       k=0; n1=-1;
       while(k==0)
       {
        k = (w&0x01) ? 1 : 0 ;
        w = w >> 1;
        n1++;
       }
       
       k=0; n2=n1;
       while(k==0)
       {
        k = (w&0x01) ? 1 : 0 ;
        w = w >> 1;
        n2++;
       }
         H[i][j] = 2*J[n1][n2];
         H[j][i] = 2*J[n1][n2];
       }
       else
       {
         H[i][j] = 0;
         H[j][i] = 0;
       }     
      }
      
      H[i][i] = 0;
      
      for(tmp=0,k=0;k<LCm;k++)
         tmp -= H[i][k] ;
      H[i][i] = tmp + sigmaJ;         
     }
     
     return H;
}

void printmat(double **mat, int dim)
{
  int i,j;
  for(i=0;i<dim;i++)
    {
      for(j=0;j<dim;j++)
        cout << setprecision(6) << mat[i][j] << "\t";
      cout << "\n" ;
    }
}


