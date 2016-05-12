#include<stdio.h>
#include<math.h>
#include <stdlib.h>
#include <omp.h>


void calcula_fuerza(double *Fx, double *Fy, double *Fz , double *xi, double *xj, double *yi, double *yj, double *zi, double *zj, double eps );
void  kick(double *p1, double *p2, double *p3, double *v1, double *v2, double *v3, double *a1, double *a2, double *a3, int n, double delta_t);
void  drift(double *p1, double *p2, double *p3, double *v1, double *v2, double *v3, double *a1, double *a2, double *a3, int n, double delta_t);
void suma_fuerza(double *F1, double *F2, double *F3 , double *X, double *Y, double *Z, int n, double epsi );
void Energia_k(double *Ek, double *v1, double *v2, double *v3, int n);
void Energia_P(double *UP, double *x1, double *y1, double *z1, int n);
double calcula_tiempo_total(int n);
double calcula_time_step(int n, double epsilon);
void recibe_input(int argc, char **argv, int *n, double *e);
void escribe_estado(double *xo, double *yo,double *zo,double *x, double *y,double *z,double *vo1, double *vo2,double *vo3,double *v1, double *v2,double *v3, double *Ek1, double *Ek2,double *UP1, double *UP2, int n, int id);

int main(int argc, char **argv)
{
int n;
double epsilon;
  recibe_input(argc, argv, &n, &epsilon);


//int n=1000;
double R;
double pi=3.14159;
double r;
R=(3.0*n/4.0*pi);
r=pow(R,1.0/3.0);

//epsilon=0.001;

//double time_step=1.0;
//int n_steps=5;


double *X, *Y, *Z,*F1, *F2, *F3, *V1, *V2, *V3, *Ek1, *Ek2, *UP1, *UP2; 

X=malloc(sizeof(double) * n); 
Y=malloc(sizeof(double) * n); 
Z=malloc(sizeof(double) * n); 
F1=malloc(sizeof(double) * n); 
F2=malloc(sizeof(double) * n); 
F3=malloc(sizeof(double) * n); 

V1=malloc(sizeof(double) * n); 
V2=malloc(sizeof(double) * n); 
V3=malloc(sizeof(double) * n); 

Ek1=malloc(sizeof(double) * n); 
Ek2=malloc(sizeof(double) * n);

UP1=malloc(sizeof(double) * n); 
UP2=malloc(sizeof(double) * n);


double *Xo,*Yo, *Zo, *Vo1, *Vo2, *Vo3;

Xo=malloc(sizeof(double) * n); 
Yo=malloc(sizeof(double) * n); 
Zo=malloc(sizeof(double) * n); 
Vo1=malloc(sizeof(double) * n); 
Vo2=malloc(sizeof(double) * n); 
Vo3=malloc(sizeof(double) * n); 

int i,j,l,m;
double rad, theta, phi;
double Fx; double Fy; double Fz;
double xi; double xj; double yi; double yj; double zi; double zj;
double TimeT, TotalT;

//============================================================================================================================================


TimeT = calcula_time_step(n, epsilon);
TotalT = calcula_tiempo_total(n);


//-----------------------Inicia posiciones----//
for(i=0; i<n; i++)
{
rad=r*(drand48());
theta=2.0*pi*(drand48());
phi=pi*drand48();
X[i]=(rad)*sin(phi)*cos(theta);
Y[i]=(rad)*sin(phi)*sin(theta);
Z[i]=(rad)*cos(phi);
//printf("%lf \t %lf \t %lf \n",X[i],Y[i],Z[i]);
}
//---------------------------------------------//

//------------------------------------------------------------Calcula fuerzas-----------------------------------------------------------------//

suma_fuerza(F1, F2, F3 , X, Y, Z, n, epsilon );


//for(i=0;i<n;i++)
//{
//printf("%lf \t %lf \t %lf  \n",F1[i],F2[i],F3[i]);
//}

//---------------------------------------------leapfrog------------------//
 
  Energia_P(UP1, X, Y, Z, n);
  kick(X, Y, Z, V1, V2, V3, F1, F2, F3, n, TimeT);

for(i=0; i<n; i++){
	Xo[i]=X[i];
	Yo[i]=Y[i];
	Zo[i]=Z[i];
	Vo1[i]=V1[i];
	Vo2[i]=V2[i];
	Vo3[i]=V3[i];
}



  Energia_k(Ek1, V1, V2, V3,  n);  
  for(i=0;i<TotalT;i++){
	drift(X, Y, Z, V1, V2, V3, F1, F2, F3, n, TimeT);  
	suma_fuerza(F1, F2, F3 , X, Y, Z, n,epsilon );
	kick(X, Y, Z, V1, V2, V3, F1, F2, F3, n, TimeT); 

  }
  Energia_k(Ek2, V1, V2, V3,  n);
  Energia_P(UP2, X, Y, Z, n);
//for(i=0;i<n;i++){
//printf("%lf \t %lf \t %lf \n",X[i],Y[i],Z[i]);
//}





//for(i=0; i<n; i++){
//	printf("%lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \t %lf \n \n",Xo[i],Yo[i],Zo[i],X[i],Y[i],Z[i],Vo1[i],Vo2[i],Vo3[i],V1[i],V2[i],V3[i],Ek1[i],Ek2[i],UP1[i],UP2[i]);
//}

escribe_estado(Xo, Yo, Zo, X, Y, Z, Vo1,  Vo2, Vo3, V1,  V2, V3,  Ek1, Ek2, UP1,  UP2,  n, i);






return 0;
}


//-----------------------------------------------------------------FUNCIONES----------------------------------------------------------------
//..........................................................................................................................................

void escribe_estado(double *xo, double *yo,double *zo,double *x, double *y,double *z,double *vo1, double *vo2,double *vo3,double *v1, double *v2,double *v3, double *Ek1, double *Ek2,double *UP1, double *UP2, int n, int id){
  FILE *out;
  char filename[512];
  int i;
  sprintf(filename, "output_%d.dat", id);
  if(!(out=fopen(filename, "w"))){
    fprintf(stderr, "Problem opening file %s\n", filename);
    exit(0);
  }
  for (i=0; i<n ; i++){
    fprintf(out, "%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n",xo[i],xo[i],zo[i],x[i],y[i],z[i],vo1[i],vo2[i],vo3[i],v1[i],v2[i],v3[i],Ek1[i],Ek2[i],UP1[i],UP2[i]);
	    	    
  }

  fclose(out);
}













//...........................................................Inicio N y E..................................................................
void recibe_input(int argc, char **argv, int *n, double *e){
  if(argc!=3){
    fprintf(stderr, "USAGE: ./nbody N_particles epsilon\n");
    exit(1);
  }
  *n = atoi(argv[1]);
  *e = atof(argv[2]);
}
//........................................................................................................................................
double calcula_time_step(int n, double epsilon){
  double t_dyn; 
  double rho;
  double R;
  double G=4.49E-3;
  double PI=3.14159;
  rho = n / (4.0/3.0) / (PI * pow(epsilon, 3));
  t_dyn = 1.0/sqrt(G * rho);
  return t_dyn;
}



//----------------------------------------------------------Calcula tiempo total----------------------------------------------------------------
double calcula_tiempo_total(int n){
  double t_dyn;
  double rho;
  double R;
  double G=4.49E-3;
  double PI=3.14159;
  R = pow(n, 1.0/3.0);
  rho = n / (4.0/3.0) / (PI * pow(R,3));
  t_dyn = 1.0/sqrt(G * rho);
  return 5 * t_dyn;
}

//----------------------------------------------------------Calcular Energia Cinetica ----------------------------------------------------------
void Energia_k(double *Ek, double *v1, double *v2, double *v3, int n){
int i;
	for(i=0; i<n; i++){
	Ek[i]=(v1[i]*v1[i]+ v2[i]*v2[i]+ v3[i]*v3[i])*0.5;
	}
}


void Energia_P(double *UP, double *x1, double *y1, double *z1, int n){
int i, j;
double dist;
double u;
double g=4.49E-3;
dist=0.0;
	
	for(i=0; i<n; i++){
	
		for(j=0; j<n; j++){
			if(i!=j){
			dist=pow((x1[i]-x1[j] + y1[i]-y1[j]+ z1[i]-z1[j]),2);
			u+=g/sqrt(dist);
			}
			
		}
		UP[i]=u;
	}
}

//-------------------------------------------------------------Calcular Fuerzas -----------------------------------------------------------------




void calcula_fuerza(double *Fx, double *Fy, double *Fz , double *xi, double *xj, double *yi, double *yj, double *zi, double *zj, double eps ){
double G=4.49E-3;
//double Fx, Fy,Fz;
double r;

double xdiff;
double ydiff;
double zdiff;

xdiff=*xi-*xj;
ydiff=*yi-*yj;
zdiff=*zi-*zj;
*Fx=-(G*xdiff)/pow((xdiff*xdiff + eps*eps),1.5);
*Fy=-(G*ydiff)/pow((ydiff*ydiff + eps*eps),1.5);
*Fz=-(G*zdiff)/pow((zdiff*zdiff + eps*eps),1.5);
}
//----------------------------------------------------------Suma fuerzas ---------------------------------------------------------------------


void suma_fuerza(double *F1, double *F2, double *F3 , double *X, double *Y, double *Z, int n, double epsi ){
double xi; double xj; double yi; double yj; double zi; double zj;
double Fx, Fy, Fz;
int j, l;

 #pragma omp parallel for private(l)
for(j=0;j<n;j++){

 for(l=0;l<n;l++){
  if(j!=l){
   xi=X[j];
   xj=X[l];
   yi=Y[j];
   yj=Y[l];
   zi=Z[j];
   zj=Z[l];
   calcula_fuerza(&Fx,&Fy, &Fz, &xi, &xj, &yi, &yj, &zi, &zj,epsi);
   //printf("%lf \t %lf \t %lf  \n \n",Fx,Fy,Fz);
   F1[l]-=Fx;
   F2[l]-=Fy;
   F3[l]-=Fz;

   }
 }
}
}









//---------------------------------------------------------------leapfrog-----------------------------------------------------------------
void  kick(double *p1, double *p2, double *p3, double *v1, double *v2, double *v3, double *a1, double *a2, double *a3, int n, double delta_t){
int i;
#pragma omp parallel for
for(i=0;i<n;i++){
      
      v1[i] += a1[i] * delta_t;
      v2[i] += a2[i] * delta_t;
      v3[i] += a3[i] * delta_t;
  }
}  

void  drift(double *p1, double *p2, double *p3, double *v1, double *v2, double *v3, double *a1, double *a2, double *a3, int n, double delta_t){
int i;
#pragma omp parallel for 
for(i=0;i<n;i++){
      p1[i] += v1[i] * delta_t;
      p2[i] += v2[i] * delta_t;
      p3[i] += v3[i] * delta_t;  
  }

}  




