
//============================================================================
// Name        : opinion.cpp
// Author      : Elena Garcia Ramos
// Version     : v1.0
// Copyright   : 
// Description : Computing the phase change in an opinion-dynamics model with 
//               separation of time scale. 
// Original Reference: Physical Review E 83.016111(2011)
// Authors Reference: G. Iniguez, et al.
//============================================================================
#include <unistd.h>
#include <getopt.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <time.h>
using namespace std;

const double PI = 3.1415926535897932384626433832795;

/* periodic boundary conditions */
int pbc(int r, int n){return (r + n)%n;}

/* gaussian random number generator */
void initgauss(double* x, int N);

/* polar form of Box-Muller Transform */
void gauss(double* u);

/* avegare opinion at short distance given by m */
double avg(double*x, int N, int m, int a);

/* avegare opinion at short distance given by m in case x[a] and x[b] are neighbours*/
double avg_neigh(double*x, int N, int m, int a, int b);

int count(double*x, int N);

double fraction(double *x, double *a, int N);

int cluster(double *x, int N);

/* apply one iteration of the euler method for differential equations */
double* euler(double *x, double* a, int n, int m, double dt, double *y);

//void exchg (double *x, int N, int m, int *a, int *b);
int exchg (double *x, double *a, int N, int m);



int main(int argc, char **argv){
  int opt;
  int N = 10, m = 1, T = 10, g = 10;
  double dt = 0.0001;
  int flag = 1, res = 1;
  while ((opt = getopt(argc, argv, "g:N:m:T:t:f")) != -1) {
    switch (opt) {
    case 'g':  g = atoi(optarg);
      break;
    case 'N': N = atoi(optarg);
      break;
    case 'm':  m = atoi(optarg);
      break;
    case 'T':  T = atoi(optarg);
      break;
    case 't':  dt = atof(optarg);
      break;
    case 'f':  flag = 0; res = 0;
      break;
    default: 
      cout <<  "Usage: " << argv[0] << "[-N agents] [-m inter] [-g Euler] [-T steps] [-t dt]"<< endl;
      exit(EXIT_FAILURE);
    }
  }
    
  if (optind > argc) {
    fprintf(stderr, "Expected argument after options\n");
    exit(EXIT_FAILURE);
  }
  cout << "N = " << N << "\n" << "g = " << g << "\n" <<  "m = " << m << "\n" << "T = " << T << "\n" << "dt = " << dt << "\nInterchage: " << flag << endl;
   
  char filename[512];
  char filename2[512];
   
  ofstream outfile;
  ofstream outfile2;
  sprintf(filename,"opinion_N%d_m%d_g%d_T%d.out",N,m,g,T);
  sprintf(filename2,"agents_N%d_m%d_g%d_T%d.out",N,m,g,T);
  outfile.open (filename);
  outfile2.open (filename2);
   
  time_t start,end;
  double* alpha;
  alpha = new double [N];
  double* x;
  x = new double [N];
  double* alpha_ini;
  alpha_ini = new double [N];
  double* x_ini;
  x_ini = new double [N];
  int und;
  time (&start);

  /* initialize rand() */
  //srand (time(NULL));
  srand (123);

  initgauss(x, N);
  initgauss(alpha, N);

  outfile << "steps  n_exchd  n_und " << endl; 
  outfile2 << "  n     x " << endl; 
  
  for (int r = 1; r <= T; r++){
    for (int n = 1 ; n <= g; n++ ){
      euler(x,alpha,N,m,dt,x);
    }
    und=count(x,N);
    if (und != 0){
      if (flag == 1 && res != 0){
	res = exchg(x,alpha,N,m);
      }
      outfile << "  "  <<  r << "     " <<  res << "    " << und << endl;
    }else{
      cout << "All agents are decided."<< endl; exit;
    }
  }
  
  for (int i = 0; i < N; i++){
    outfile2 << i <<  "   "<<  x[i] << " " << alpha[i] << endl;
  }
  
  double frac=fraction(x,alpha,N);
  double clus=cluster(x,N);
  cout << frac << "  " << clus << endl;
  
  time (&end);
  double dif = difftime (end,start);
  printf ("Elasped time is %.2lf seconds. \n", dif );
  
  outfile.close();   
  outfile2.close();   

  delete[] x;
  delete[] alpha;

  return 0;
} 
 
int exchg (double *x, double *a,int N, int m){
  int r = 0;
  int x1, x2;
  double avgx1, avgx2;
  double p, q, tmp;
  for (int i = 0; i < N*N/2; i++){
    avgx1=0; avgx2=0;
    p=0; q=0; tmp=0; 
    x1 = (rand()/(double)RAND_MAX)*N;
    x2 = (rand()/(double)RAND_MAX)*N;
    if( abs(x[x1]) < 1 && abs(x[x2]) < 1 && x1 != x2){
      avgx1=avg(x,N,m,x1);
      avgx2=avg(x,N,m,x2);
      p=1/4.*(abs(x[x1]-avgx1)+abs(x[x2]-avgx2));
      if ( x2 >= pbc(x1-m,N)  && x2 <= pbc(x1+m,N) ){
	avgx1=avg_neigh(x,N,m,x1,x2);
	avgx2=avg_neigh(x,N,m,x2,x1);
	q=1/4.*(abs(x[x1]-avgx2)+abs(x[x2]-avgx1));
      }else
	q=1/4.*(abs(x[x1]-avgx2)+abs(x[x2]-avgx1));
      if( p > q ){
	tmp=x[x1];
	x[x1]=x[x2];
	x[x2]=tmp;
	tmp=a[x1];
	a[x1]=a[x2];
	a[x2]=tmp;
	r++;  
      }
    }
  }
  return r;
}


double* euler(double *x, double* a, int n, int m, double dt, double *y){
  double* fs;
  fs = new double [n];
  double* fl;
  fl = new double [n];
  double sum =0;
  for (int i = 0; i < n; i++){ sum+=x[i]; }
  for (int i = 0; i < n; i++){
    if (abs(x[i]) < 1){
      fs[i]=0;
      fl[i]=0;
      for (int j = 1; j <= m ; j++){
	fs[i] += x[pbc(i+j,n)]+x[pbc(i-j,n)];
      }
      fl[i]=sum-fs[i]-x[i];
      y[i]=x[i]+dt*(abs(x[i])*fs[i]/2/m+a[i]*fl[i]/(n-1-2*m));
      if ( abs(y[i]) > 1 ){y[i]=y[i]/abs(y[i]);}
    }
  } 
  delete[] fs;
  delete[] fl;
  return y;
}

double fraction(double *x, double *a, int N){
  int counter=0;
  double frac;
  for (int i = 0 ; i < N; i++ ){
    if(abs(x[i]) < 1){
      if(a[i] < 0){counter++;}}
  }  
  frac = (double)counter/N;
}

int cluster(double *x, int N){
  int counter=0;
  double a,b;
  a=x[0]/abs(x[0]);
  if(x[N-1]/abs(x[N-1]) != a){counter=1;}
  for (int i = 0 ; i < N; i++ ){
    b=x[i]/abs(x[i]);
    if (a != b){
      counter++;
    }
    a=b;    
  }
  return counter;
}

void initgauss(double *x, int N){
  double r[2];
  for (int i = 0; i < N; i += 2) {    
  repeat:
    gauss(r);
    if ( abs(r[0]) > 1 || abs(r[1]) > 1){
      goto repeat;
    }
    x[i] = r[0];
    x[i+1] = r[1];
  }
}

void gauss(double *u){
  double l = 10, w;
  while( l >= 1.0){
    u[0] = (rand()/(double)RAND_MAX)*2-1;
    u[1] = (rand()/(double)RAND_MAX)*2-1;
    l = u[0]*u[0]+u[1]*u[1];
  }
  w = sqrt( (-2.0 * log( l ) ) / l );
  u[0]=(u[0] * w);
  u[1]=(u[1] * w);
}

double avg(double*x, int N, int m, int a){
  double avg = 0;
  for (int j = 1; j <= m ; j++){
    avg += x[pbc(a+j,N)]+x[pbc(a-j,N)];
  }
  avg=avg/2/m;
   return avg;
}

double avg_neigh(double*x, int N, int m, int a, int b){
  double avg = 0;
  for (int j = 1; j <= m ; j++){
    avg += x[pbc(a+j,N)]+x[pbc(a-j,N)];
  }
  avg=(avg+x[a]-x[b])/2/m;
   return avg;
}

int count(double*x, int N){
  int und = 0;
  for (int i =0; i< N; i++){
    if (abs(x[i]) < 1){und++;}}
  return und;
}
