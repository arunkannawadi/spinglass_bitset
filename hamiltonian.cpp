#include<iostream>
#include<iomanip>
#include<bitset>
#include<stdlib.h>
#include<math.h>
#include <Eigen/Dense>
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
  double **symm_gauss_rand_mat_gen(int);
  MatrixXd hamiltonian(double**, int,int);
  MatrixXd diagonalize_symmetric(MatrixXd);
  void printmat(double**,int); 
  void printmat_eigen(MatrixXd,int);
  
  double i,j;
  double **J;
  MatrixXd H, V;
  
  J = symm_gauss_rand_mat_gen(4);
      
  H = hamiltonian(J, 2, 4);
  
  cout << "Random matrix J = \n";
  printmat(J,4);

  cout << " Hamiltonian H = \n" << H << endl;
  
  V = diagonalize_symmetric(H);
  cout << "\n \n The eigenvecotrs are: \n " << V << endl;
  
  
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

/*
long double *definite_particle_numbers(int m)
{
 int i,k, LCm;
 int nchoosek(int,int);
 
 LCm = nchoosek(Nqubits,m);
 long double states[(const unsigned)LCm]; 
 long double v,w,t;
  
 v = (1 << m) - 1;  // will this work?
 k = 0;
 
 while(k<LCm)
 {
  states[k++] = v;
  t = v | (v-1);
  v = (t+1) | (((~t & -~t) -1) >> (__builtin_ctz(v)+1));
 }
 
 return states;
 
}
*/

bitset<Nqubits>* definite_particle_states(int m)
{
 bitset<Nqubits> u;
 int LCm;
 int nchoosek(int,int);
 LCm = nchoosek(Nqubits,m);
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

MatrixXd hamiltonian(double** J, int m, int L)
{
  bitset<Nqubits> *definite_particle_states(int);
  int nchoosek(int,int);
  
  const unsigned int constL = L;  
  int LCm; 
  int i,j,k,n1,n2,n0;
  bitset<Nqubits> u,v,w,x,y;
  bitset<Nqubits> one;
  one.set(0);
    
  double tmp, sigmaJ;
  sigmaJ = -1;
  for(i=1;i<L;i++)
     for(j=0;j<i;j++)
       sigmaJ = sigmaJ + J[i][j];
       
  cout << "sigmaJ = " << sigmaJ << endl;
    
  LCm = nchoosek(L,m);
  
  bitset<Nqubits> *states;
  states =  definite_particle_states(m);
  
  MatrixXd H = MatrixXd::Zero(LCm,LCm);
  
   for(i=0;i<LCm;i++)
     {
      for(j=LCm-1;j>i;j--)
      {
       u = states[i];
       v = states[j];       
       
       for(n1=0;u[n1]==v[n1];n1++);
       for(n2=n1+1;u[n2]==v[n2];n2++);

       for(n0=n2+1;n0<L;n0++)
         if(!(u[n0]==v[n0]))
          n0=L+2;
           
       if(n0==L)
       {       
         H(i,j) = 2*J[n1][n2];
         H(j,i) = 2*J[n1][n2];
       }  
      }
      
      H(i,i) = 0;
      
      for(tmp=0,k=0;k<LCm;k++)
         tmp -= H(i,k) ;
      H(i,i) = tmp + sigmaJ;         
     }
     
     return H;
}

MatrixXd diagonalize_symmetric(MatrixXd A)
{
 if(!(A==A.transpose()))
 {
  cout << "\n The matrix is not symmetric and will be symmetrized" << endl;
  if(!(A.rows()==A.cols()))
  {
    cout << "\n Error! The matrix is rectangular. Returning the original matrix. " << endl;
    return A;
  }
   A = 0.5*(A + A.transpose()); // ensure Matrix is symmetric 
 }

 SelfAdjointEigenSolver<MatrixXd> es(A);
 
 MatrixXd All;
 All.resize(A.rows(),A.cols()+1);
 All << es.eigenvectors(), es.eigenvalues();
 return All;
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


void printmat_eigen(MatrixXd mat, int dim)
{
  int i,j;
  for(i=0;i<dim;i++)
    {
      for(j=0;j<dim;j++)
        cout << setprecision(6) << mat(i,j) << "\t";
      cout << "\n" ;
    }
}

