#include<iostream>
#include<ctime>
#include<fstream>
#include<iomanip>
#include<bitset>
#include<stdlib.h>
#include<math.h>
#include<sstream>
#include<string>
#include<Eigen/Dense>
#include<Eigen/Eigenvalues>

using namespace std; 
using namespace Eigen;

#define root2 1.414213562373095
#define Nqubits 50

int main()
{
  /* For random number generation */
  struct timeval start;
  struct timeval end;
  double t1,t2,t3 ,t4;
  time_t seconds, basis_start, basis_end, conc_start, conc_end;
  time(&seconds);
  srand((unsigned int)seconds);
  
  /* Function declarations */
  double **symm_gauss_rand_mat_gen(int);
  MatrixXd hamiltonian(double**, bitset<Nqubits>*,int,int,int);
  int nchoosek(int,int);
  bitset<Nqubits>* definite_particle_basis(int,int);
  double PR(VectorXd);
  void average_concurrence(double*,VectorXd, bitset<Nqubits>*,int,int);
  double LN3(VectorXd, bitset<Nqubits>*,int,int);
  
  /* Variable declarations */
  int i,j,k,N,LCm,L=Nqubits,m=2,Niter=1;
  LCm = nchoosek(L,m);
  
  double ln3, pr, **J;
  double *entanglement = new double[3];
  MatrixXd H;
  
  cout << "SpinGlass program started with L = " << L << " m = " << m << " and Niter = " << Niter << endl; 
  /* Basis generation */
  bitset<Nqubits> *basis;
  time(&basis_start);
  basis = definite_particle_basis(m,LCm);
  time(&basis_end);
  cout << " Basis calculation done in " << difftime(basis_end,basis_end) << endl;
    
  /* File opening */
  ofstream myfile;
  ostringstream strL, strm;
  string filename = "SG_L";

  strL << L; strm << m;
  filename += strL.str();
  filename += "_m";
  filename += strm.str();
  filename += ".txt";
  
  myfile.open(filename.c_str(),fstream::out | fstream::app);
  cout << "File opened \n";
  
  for(N=0;N<Niter;N++)
  {
   cout << " N = " << Niter << "\t";
   J = symm_gauss_rand_mat_gen(L);
   cout << " J created \t" ;
   H = hamiltonian(J,basis,m,L,LCm);
   cout << " H created \t" ;
     
   SelfAdjointEigenSolver<MatrixXd> ES(H);
   cout << " H diagonalized \n";
   
   time(&conc_start);
   for(k=0;k<LCm;k++)
   {
    average_concurrence(entanglement, ES.eigenvectors().col(k),basis,m,LCm);
    cout << "Concurrence for state " << k << "\t";
    ln3 = LN3(ES.eigenvectors().col(k),basis,m,LCm);
    cout << "LN3 for state " << k << " is calculated \t";
    pr = PR(ES.eigenvectors().col(k));
    cout << "PR for state " << k << "is calculated \n";
    
    myfile << pr << "\t" << entanglement[0] << "\t" << entanglement[1] << "\t" << entanglement[2] << "\t" << ln3 << endl;
    
   }
   time(&conc_end);
   
   cout<< "Time to compute entanglement = " << difftime(conc_end,conc_start);
   
   for(k=0;k<L;k++)
    free(J[k]);
   free(J);
   
  }
  myfile.close();
}

void average_concurrence(double *entanglement, VectorXd state, bitset<Nqubits>* basis, int m, int LCm)
{
  int i,j,k,l,L=Nqubits;
  int nchoosek(int,int);
        
  double v,w,x,y,z;
  double C, avgC, avgTangle, LN2, avgLN2, sumLN2=0,sumTangle=0,sumC=0; 
  for(i=1;i<L;i++)
   for(j=0;j<i;j++)
   {
    v=0; y=0; z=0; w=0; x=0;
    for(k=0;k<LCm;k++)
     {
      if((basis[k][i]==0)&(basis[k][j]==0))
       { v += state(k)*state(k); }
      else if((basis[k][i]==0)&(basis[k][j]==1))
        x += state(k)*state(k);
      else if((basis[k][i]==1)&(basis[k][j]==0))
        w += state(k)*state(k);
      else if((basis[k][i]==1)&(basis[k][j]==1))    
        y += state(k)*state(k);

      basis[k].flip(i); basis[k].flip(j);    
      for(l=0;l<k;l++)
       {
        if(basis[k]==basis[l])
          z += state(k)*state(l);
        } 
       basis[k].flip(i); basis[k].flip(j);
       
     }
     C = ((abs(z)-sqrt(v*y))>0) ? 2.0*(abs(z)-sqrt(v*y)) : 0.0 ;
     sumC += C;
     sumTangle += C*C;
     
     LN2 = log(w+x+0.5*(v+y+sqrt(4*z*z+(v-y)*(v-y)))+0.5*abs(v+y-sqrt(4*z*z+(v-y)*(v-y))));
     sumLN2 += LN2;
   }
   
    int LC2; LC2 = nchoosek(L,2);
    
    avgC = sumC/LC2;
    avgTangle = sumTangle/LC2;
    avgLN2 = sumLN2/LC2;
    
    
    // double *entanglement = new double[3];
    entanglement[0] = avgC;
    entanglement[1] = avgTangle;
    entanglement[2] = avgLN2;
    
    // return entanglement;
}

double LN3(VectorXd state, bitset<Nqubits>* basis, int m, int LCm)
{
     
 SelfAdjointEigenSolver<MatrixXd> es(9);
 bitset<3> B,B1, B2;
 int b,b1, b2, n1,n2,p,i,j,k,l1,l2;
 int nchoosek(int,int);
 int LC3, L=Nqubits; LC3 = nchoosek(L,3);
 
 double Lambda, LN3, avgLN3, sumLN3=0;
 MatrixXd rho = MatrixXd::Zero(9,9);
 VectorXd eigenvals;
 
 /* i,j,k - Qubit indices */
 /* l1,l2  - Basis indices */
 
 for(i=2;i<L;i++)
  for(j=1;j<i;j++)
   for(k=0;k<j;k++)
   {
    for(p=0;p<3;p++)
    {
     for(n1=0;n1<9;n1++)
      for(n2=0;n2<9;n2++) 
       rho(n1,n2) = 0;
     
     for(l1=0;l1<LCm;l1++)
     {
      B1.set(0,basis[l1][k]);
      B1.set(1,basis[l1][j]);
      B1.set(2,basis[l1][i]);
     
      basis[l1].reset(i); basis[l1].reset(j), basis[l1].reset(k);
     
      b1 = (int)B1.to_ulong();
      rho(b1,b1) += state(l1)*state(l1);
     
      for(l2=0;l2<l1;l2++)
      {
       B2.set(0,basis[l2][k]);
       B2.set(1,basis[l2][j]);
       B2.set(2,basis[l2][i]);
      
       basis[l2].reset(i); basis[l2].reset(j), basis[l2].reset(k);
      
       if(basis[l2]==basis[l1])
       {
        if(!(B2[p]==B1[p]))
        { B2[p].flip(); B1[p].flip(); }
         
        b1 = (int)B1.to_ulong();
        b2 = (int)B2.to_ulong();
        rho(b1,b2) += state(l1)*state(l2);
        rho(b2,b1) = rho(b1,b2);
        
        if(!(B2[p]==B1[p]))
         { B2[p].flip(); B1[p].flip(); }
       }
       basis[l2].set(i,B2[2]); basis[l2].set(j,B2[1]); basis[l2].set(k,B2[0]); 
       B2.reset();
      }
      basis[l1].set(i,B1[2]); basis[l1].set(j,B1[1]); basis[l1].set(k,B1[0]);              
      B1.reset();
     }
              
     /* Diagonalize PPT here */
     es.compute(rho,EigenvaluesOnly);
     eigenvals = es.eigenvalues();

     for(Lambda=0,n1=0;n1<9;n1++)
       { Lambda += abs(eigenvals[n1]); }
       
     LN3 = log(Lambda);
     sumLN3 += LN3;
    }     
   }
   avgLN3 = sumLN3/(3*LC3);
   return avgLN3;
}


double PR(VectorXd state)
{
 double ipr, pr;
 int k,n;
 
 n = state.size();
 for(ipr=0,k=0;k<n;k++)
  ipr += pow(state(k),4);
 pr = 1.0/ipr;
 return pr;
}

MatrixXd hamiltonian(double** J, bitset<Nqubits> *basis, int m, int L,int LCm)
{
 
  const unsigned int constL = L;  
  int i,j,k,n1,n2,n0;
  bitset<Nqubits> u,v,w,x,y;
  bitset<Nqubits> one;
  one.set(0);
    
  double tmp, sigmaJ;
  sigmaJ = -1;
  for(i=1;i<L;i++)
     for(j=0;j<i;j++)
       sigmaJ = sigmaJ + J[i][j];
           
  MatrixXd H = MatrixXd::Zero(LCm,LCm);
  
   for(i=0;i<LCm;i++)
     {
      for(j=LCm-1;j>i;j--)
      {
       u = basis[i];
       v = basis[j];       
       
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

/* For generating J */
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

bitset<Nqubits>* definite_particle_basis(int m, int LCm)
{
 bitset<Nqubits> u;
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
