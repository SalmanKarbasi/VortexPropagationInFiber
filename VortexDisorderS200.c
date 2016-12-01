////////////////////////////////////////////////////////
////*********Transparent boundary condition*********////
////////////////////////////////////////////////////////
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<complex.h>
#include<omp.h>
#include <string.h>
#include <sysexits.h>
/******************************************************/
double randgen()
{
double r, R, t, b;
r=rand();
R=RAND_MAX;
t=r/R;
if (t>=0.5)
	{
	b=1;
	}
	else
	{
	b=0;
	}
return(b);
}
/*****************************************************/
main()
{
srand ((unsigned)time(NULL));
char buffer[1];
char bufferImag[100];
long int   a,jj,ii,gg,ff,mm,Nx,Ny,Nz,Ncw,nc;
double x0, y0, pi, lam, k0, w0, wr, w, wc, dx, dy, dc, n0, dz, na, ng, rnd, Nsites, nn=1;
/*****************************************************/
lam = 0.633;        /*wavelength*/
w0 = 4;    	/* Gaussian pulse width */
wr = pow(2,-0.5)*w0;
w = 200; 	/*width of the simulation window*/
dc = 0.8;
//Nsites = w/dc;  
Nsites =  250; //(int) Nsites; 
dx = 0.2;   /*step size in X direction*/
nc = 4;
//nc = (int)nc;
dy = dx;
//w = Nsites*nc*dx;	/*simulation window width*/
Nx = 1000;//(int)(w/dx);      /*n number of samples in "x" direction*/
Ny = Nx;
//Ncw = (int)(wc/dx);
printf("Nx = %ld, dc = %f, dx = %f, Nsites = %f, w = %f\n", Nx, dc, dx, Nsites, w);
double (*n)[Nx] = malloc(Nx*sizeof(double[Nx]));

    if (n == NULL)
    {
	fputs("Failed to malloc() array n.\n", stderr);
	exit(EX_UNAVAILABLE);
    }

double complex (*k1)[Nx] = malloc(2*Nx*sizeof(double[Nx]));

    if (k1 == NULL)
    {
	fputs("Failed to malloc() array n.\n", stderr);
	exit(EX_UNAVAILABLE);
    }

double complex (*k2)[Nx] = malloc(2*Nx*sizeof(double[Nx]));

    if (k2 == NULL)
    {
	fputs("Failed to malloc() array n.\n", stderr);
	exit(EX_UNAVAILABLE);
    }

double complex (*k3)[Nx] = malloc(2*Nx*sizeof(double[Nx]));

    if (k3 == NULL)
    {
	fputs("Failed to malloc() array n.\n", stderr);
	exit(EX_UNAVAILABLE);
    }

double complex (*k4)[Nx] = malloc(2*Nx*sizeof(double[Nx]));

    if (k4 == NULL)
    {
	fputs("Failed to malloc() array n.\n", stderr);
	exit(EX_UNAVAILABLE);
    }

double complex  kx1[2*Nx];
double complex  ky1[2*Nx];

/*(*kx1)[Nx] = malloc(2 * Nx * sizeof (double));

    if (kx1 == NULL)
    {
	fputs("Failed to malloc() array n.\n", stderr);
	exit(EX_UNAVAILABLE);
    }
/*/
/*
double complex (*ky1)[Nx] = malloc(2 * Nx * sizeof (double));

    if (ky1 == NULL)
    {
	fputs("Failed to malloc() array n.\n", stderr);
	exit(EX_UNAVAILABLE);
    }
/*/
/*double (*E0)[Nx] = malloc(Nx*sizeof(double[Nx]));/*The real definiations*/
//double (*A1)[Nx] = malloc(Nx*sizeof(double[Nx]));


double complex (*E0)[Nx] = malloc(2*Nx*sizeof(double[Nx]));
double complex (*A1)[Nx] = malloc(2*Nx*sizeof(double[Nx]));


    if (E0 == NULL)
    {
	fputs("Failed to malloc() array n.\n", stderr);
	exit(EX_UNAVAILABLE);
    }

double complex (*E1)[Nx] = malloc(2*Nx*sizeof(double[Nx]));

    if (E1 == NULL)
    {
	fputs("Failed to malloc() array n.\n", stderr);
	exit(EX_UNAVAILABLE);
    }

double complex (*E2)[Nx] = malloc(2*Nx*sizeof(double[Nx]));

    if (E2 == NULL)
    {
	fputs("Failed to malloc() array n.\n", stderr);
	exit(EX_UNAVAILABLE);
    }

/*******************************************************/
na = 1.49;        /*refractive index of Air */
ng = 1.59;        /*refractive index of core*/
/*define  the incident boundary condition */
x0 = (Nx/2)*dx;  /*center of Gaussian pulse*/
y0 = (Ny/2)*dx;
//printf("real center=%f\n",x0); 
//defining a Vortex 
for(ii=0;ii<=Nx-1;ii++)
    for (jj=0;jj<=Ny-1;jj++)
               E0[ii][jj]=exp(-(pow((ii+1)*dx-x0,2)+pow((jj+1)*dy-y0,2))/pow(w0,2));

for(ii=0;ii<=Nx-1;ii++)
    for (jj=0;jj<=Ny-1;jj++)
		A1[ii][jj] = sqrt(pow((ii+1)*dx-x0,2)+pow((jj+1)*dy-y0,2))/wr*cexp(I*atan2(((jj+1)*dy-y0),((ii+1)*dx-x0)));
//(jj*dy-y0)/(ii*dx-x0)

for(ii=0;ii<=Nx-1;ii++)
    for (jj=0;jj<=Ny-1;jj++)
	E0[ii][jj] = E0[ii][jj]*A1[ii][jj];	

/********************************************************/
for (gg=0;gg<Nsites;gg++)
	{
		for (ff=0;ff<Nsites;ff++)
			{
		rnd=randgen();
		for (ii=0;ii<=nc-1;ii++)
			     {   
				   for(jj=0;jj<=nc-1;jj++)   
				    	{
			        	        n[ii+nc*gg][jj+nc*ff]=ng+rnd*(na-ng);
				    	}
			     }
			}
	}
/**********************Making a step index refractice index profile*******
for(ii=0;ii<Nx-1;ii++){	
	for (jj=0;jj<Ny-1;jj++){
		if (pow((ii*dx-x0),2)+pow((jj*dy-y0),2)<=pow(30*dx,2))
			n[ii][jj] = 1.5; 
	}}

/***************************Ref Index Profile******************************/
        FILE *f;
        f = fopen("N5E.dat","w");
        for (ii=0;ii<=Nx-1;ii++)
            {
             fprintf(f,"\n");   
            for(jj=0;jj<=Ny-1;jj++)
            {
                fprintf(f,"%f", n[ii][jj]);
               // fprintf(f,"%f", cabs(E2[ii][jj]));
                fprintf(f,"  ");
            }
            }   
                fclose(f);



/*
//#pragma omp parallel for 
     for (ii=0;ii<=Nx-1;ii++){
//#pragma omp parallel for 
      for (jj=0;jj<=Ny-1;jj++){
		if (pow(ii-Nx/2,2)+pow(jj-Nx/2,2)<pow(15,2))
		{  n[ii][jj]=ng;}
		else 
		  {n[ii][jj]=na;} 	
	}
	}*/
/*********************************************************/
pi = 4*atanl(1);		/*pi*/
k0 = 2*pi/lam; 			/*wavevector*/
printf("phase=%.4f\n",atan2(1,-1)*180/pi); //atan2(y,x) [-pi,pi]
n0 = 0.5*(na+ng);	 	/*effective refractive index of media*/
dz = 0.02*k0*pow(dx,2)*n0;
printf("dz=%f\n",dz);	
Nz = (int)(50000/dz);        	/*number of samples in "z" direction */ 
double *SigmaS = malloc(Nz * sizeof (double));

    if (SigmaS == NULL)
    {
	fputs("Failed to malloc() array n.\n", stderr);
	exit(EX_UNAVAILABLE);
    }

double *Sigmax = malloc(Nz * sizeof (double));

    if (Sigmax == NULL)
    {
	fputs("Failed to malloc() array n.\n", stderr);
	exit(EX_UNAVAILABLE);
    }

double *Sigmay = malloc(Nz * sizeof (double));

    if (Sigmay == NULL)
    {
	fputs("Failed to malloc() array n.\n", stderr);
	exit(EX_UNAVAILABLE);
    }

double *power = malloc(Nz * sizeof (double));

    if (power == NULL)
    {
	fputs("Failed to malloc() array n.\n", stderr);
	exit(EX_UNAVAILABLE);
    }

double *Px = malloc(Nz * sizeof (double));

    if (Px == NULL)
    {
	fputs("Failed to malloc() array n.\n", stderr);
	exit(EX_UNAVAILABLE);
    }

double *Py = malloc(Nz * sizeof (double));

    if (Py == NULL)
    {
	fputs("Failed to malloc() array n.\n", stderr);
	exit(EX_UNAVAILABLE);
    }
double *Ps = malloc(Nz * sizeof (double));

    if (Ps == NULL)
    {
	fputs("Failed to malloc() array n.\n", stderr);
	exit(EX_UNAVAILABLE);
    }
double *N = malloc(Nz * sizeof (double));

    if (N == NULL)
    {
	fputs("Failed to malloc() array n.\n", stderr);
	exit(EX_UNAVAILABLE);
    }

/***********************************************************/
for (ii=0;ii<=Nx-1;ii++)
    for(jj=0;jj<=Ny-1;jj++)
     E1[ii][jj]=E0[ii][jj];


//Nz=10;
/****************************************************************/
/********************Runge-Kutta method**************************/
/****************************************************************/
//#pragma omp parallel for 
for (mm=1;mm<=Nz;mm++)
    {
printf("mm=%f\n",mm*dz);	
Px[mm-1]=0;
Py[mm-1]=0; 
Ps[mm-1]=0; 
N[mm-1]=0; 
power[mm-1]=0;
SigmaS[mm-1]=0;
Sigmax[mm-1]=0;
Sigmay[mm-1]=0;

     for (ii=1;ii<=Nx-2;ii++)
#pragma omp parallel for 
      for (jj=1;jj<=Ny-2;jj++)
                  k1[ii][jj] = 
(-I*nn/n0/k0/2)*((E1[ii+1][jj]-2*E1[ii][jj]+E1[ii-1][jj])/pow(dx,2)+(E1[ii][jj+1]-2*E1[ii][jj]+E1[ii][jj-1])/pow(dy,2)+(pow(n[ii][jj],2)-pow(n0,2))*pow(k0,2)*E1[ii][jj]);


    for(ii=0;ii<=Nx-1;ii++)
#pragma omp parallel for 
        for(jj=0;jj<=Ny-1;jj++)
            E2[ii][jj]=E1[ii][jj]+0.5*dz*k1[ii][jj]; 
/***********************************************/


       for (ii=1;ii<=Nx-2;ii++)
#pragma omp parallel for 
        for(jj=1;jj<=Ny-2;jj++)
        k2[ii][jj] = (-I*nn/n0/k0/2)* 
((E2[ii+1][jj]-2*E2[ii][jj]+E2[ii-1][jj])/pow(dx,2)+(E2[ii][jj+1]-2*E2[ii][jj]+E2[ii][jj-1])/pow(dy,2)+(pow(n[ii][jj],2)-pow(n0,2))*pow(k0,2)*E2[ii][jj]);

        for (ii=0;ii<=Nx-1;ii++)
#pragma omp parallel for 
                  for(jj=0;jj<=Ny-1;jj++)      
                        E2[ii][jj]=E1[ii][jj]+0.5*dz*k2[ii][jj];

/***********************************************/


       for (ii=1;ii<=Nx-2;ii++)
#pragma omp parallel for 
        for (jj=1;jj<=Ny-2;jj++)
                k3[ii][jj] = (-I*nn/n0/k0/2)* 
((E2[ii+1][jj]-2*E2[ii][jj]+E2[ii-1][jj])/pow(dx,2)+(E2[ii][jj+1]-2*E2[ii][jj]+E2[ii][jj-1])/pow(dy,2)+(pow(n[ii][jj],2)-pow(n0,2))*pow(k0,2)*E2[ii][jj]);      

                for (ii=1;ii<=Nx-2;ii++)
#pragma omp parallel for 
             for (jj=1;jj<=Ny-2;jj++)
                                E2[ii][jj]=E1[ii][jj]+ dz*k3[ii][jj]; 

/***********************************************/


for (ii=1;ii<=Nx-2;ii++)
#pragma omp parallel for 
    for (jj=1;jj<=Ny-2;jj++)
                k4[ii][jj] = (-I*nn/n0/k0/2)* 
((E2[ii+1][jj]-2*E2[ii][jj]+E2[ii-1][jj])/pow(dx,2)+(E2[ii][jj+1]-2*E2[ii][jj]+E2[ii][jj-1])/pow(dy,2)+(pow(n[ii][jj],2)-pow(n0,2))*pow(k0,2)*E2[ii][jj]);

                for (ii=1;ii<=Nx-2;ii++)
#pragma omp parallel for 
            for (jj=1;jj<=Ny-2;jj++)                                
E2[ii][jj]=E1[ii][jj]+(k1[ii][jj]+2*k2[ii][jj]+2*k3[ii][jj]+k4[ii][jj])*dz/6;

/*************TBC Upper**************************/

for (jj=0;jj<=Ny-1;jj++)
        {   
        if (E2[2][jj]!=0){
        kx1[jj]=I/dx*clog(E2[1][jj]/E2[2][jj]);
        if (creal(kx1[jj])<0)
            { kx1[jj]=0;}
                E2[0][jj]=E2[1][jj]*cexpl(-I*kx1[jj]*dx);
        }}                  

/*************TBC Lower*************************/

for (jj=0;jj<=Ny-1;jj++)
        {   
        if (E2[Nx-3][jj]!=0){
        kx1[jj]=I/dx*clog(E2[Nx-2][jj]/E2[Nx-3][jj]);
                if (creal(kx1[jj])<0)
                        { kx1[jj]=0;}
                E2[Nx-1][jj]=E2[Nx-2][jj]*cexpl(-I*kx1[jj]*dx);
        }}

/*************TBC Left**************************/

for (ii=0;ii<=Nx-1;ii++)
        {
        if (E2[ii][2]!=0){   
        ky1[ii]=I/dy*clog(E2[ii][1]/E2[ii][2]);
                if (creal(ky1[ii])<0)
                        { ky1[ii]=0;}
                E2[ii][0]=E2[ii][1]*cexpl(-I*ky1[ii]*dx);
        }}

/*************TBC Right**************************/

for (ii=0;ii<=Nx-1;ii++)
        {
        if(E2[ii][Ny-3]!=0){   
        ky1[ii]=I/dy*clog(E2[ii][Ny-2]/E2[ii][Ny-3]);
                if (creal(ky1[ii])<0)
                        { ky1[ii]=0;}
                E2[ii][Ny-1]=E2[ii][Ny-2]*cexpl(-I*ky1[ii]*dx);
        }}


/*******************************************************/
 if (mm==1 || mm%1000==0){
	sprintf(buffer, "%ld", mm);
	strcat(buffer,".dat");	 
/****************************************************/
        FILE *ff1;
        ff1 = fopen(buffer,"w");
        for (ii=0;ii<=Nx-1;ii++)
            {
             fprintf(ff1,"\n");   
             for(jj=0;jj<=Ny-1;jj++)
            {
                fprintf(ff1,"%f", cabs(E1[ii][jj]));
                fprintf(ff1,"  ");
            }
            }   
                fclose(ff1);
/*********************************************/
	sprintf(bufferImag, "%ld", mm);
	strcat(bufferImag,"phase.dat");	 

        FILE *ffimag;
        ffimag = fopen(bufferImag,"w");
        for (ii=0;ii<=Nx-1;ii++)
            {
             fprintf(ffimag,"\n");   
             for(jj=0;jj<=Ny-1;jj++)
            {
                fprintf(ffimag,"%f", carg(E1[ii][jj]));
                fprintf(ffimag,"  ");
            }
            }

                fclose(ffimag);
/*********************************************/
	    }
/***********************************************/
        for (ii=0;ii<=Nx-1;ii++)
            for (jj=0;jj<=Ny-1;jj++)
                E1[ii][jj]=E2[ii][jj];

/***********************************************/
FILE *fff;
fff = fopen("counter.dat", "w");
fprintf(fff,"%f", mm*dz);
fprintf(fff,"\n");
fclose(fff); 

}

/*********************************************/

}
