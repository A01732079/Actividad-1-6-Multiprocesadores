#include <stdio.h>
#include <math.h>
#include <time.h>

FILE *fptr1;
FILE *fptr2;
FILE *fptr3;
FILE *fptr4;

int main(){
   long N = 50000;
   fptr1=fopen("Runge1Sec.txt","w");
   fptr2=fopen("Runge2Sec.txt","w");
   fptr3=fopen("Runge3Sec.txt","w");
   fptr4=fopen("Runge4Sec.txt","w");

   //printf("Numero de pasos:%d Atendido por thread:%d\n", N,omp_get_thread_num());
   fprintf(fptr1, "Datos que encuentra el metodo de Runge(variable ind.\t variable dep.\t numero de thread)\n");
   fprintf(fptr2, "Datos que encuentra el metodo de Runge(variable ind.\t variable dep.\t numero de thread)\n");
   fprintf(fptr3, "Datos que encuentra el metodo de Runge(variable ind.\t variable dep.\t numero de thread)\n");
   fprintf(fptr4, "Datos que encuentra el metodo de Runge(variable ind.\t variable dep.\t numero de thread)\n");

   float h,w,ab,t;
   float a=0.0, b=3.1416;
   float w0=b/4;
   float k1,k2,k3,k4;
   int i;
   w = w0;
   h = (b-a)/N;

   fprintf(fptr1, "%f\t %f\n", a, w);
   fprintf(fptr2, "%f\t %f\n", a, w);
   fprintf(fptr3, "%f\t %f\n", a, w);
   fprintf(fptr4, "%f\t %f\n", a, w);

   clock_t start = clock();

   for(i=0; i<N; i++){ //y'=te^(3t)-2y
      t=a+(h*i);
      k1 = h*(w*exp(3*t)-2*w);
      k2 = h*((w+0.5*k1)*exp(3*(t+0.5*h))-2*(w+0.5*k1));
      k3 = h*((w+0.5*k2)*exp(3*(t+0.5*h))-2*(w+0.5*k2));
      k4 = h*((w+k3)*exp(3*(t+h))-2*(w+k3));
      w = w + (1.0/6.0)*(k1 + k2 + k3 + k4);
      fprintf(fptr1, "%f\t %f\n", t, w);
   }

   w = w0;

   for(i=0; i<N; i++){ //y'=1+(t-y)^2
      t=a+(h*i);
      k1 = h*(1+pow((t-w),2));
      k2 = h*(1+pow(((t+0.5*h)-(w+0.5*k1)),2));
      k3 = h*(1+pow(((t+0.5*h)-(w+0.5*k2)),2));
      k4 = h*(1+pow(((t+h)-(w+k3)),2));
      w = w + (1.0/6.0)*(k1 + k2 + k3 + k4);
      fprintf(fptr2, "%f\t %f\n", t, w);
   }

   w = w0;

   for(i=0; i<N; i++){ //y'=1+y/t
      t=a+(h*i);
      k1 = h*(1+w/t);
      k2 = h*(1+(w+0.5*k1)/(t+0.5*h));
      k3 = h*(1+(w+0.5*k2)/(t+0.5*h));
      k4 = h*(1+(w+k3)/(t+h));
      w = w + (1.0/6.0)*(k1 + k2 + k3 + k4);
      fprintf(fptr3, "%f\t %f\n", t, w);
   }

   w = w0;

   for(i=0; i<N; i++){ //y'= cos(2ty)+sen(3ty)
      t=a+(h*i);
      k1 = h*(cos(2.0*t*w)+sin(3.0*t*w));
      k2 = h*(cos(2.0*(t+0.5*h)*(w+0.5*k1))+sin(3.0*(t+0.5*h)*(w+0.5*k1)));
      k3 = h*(cos(2.0*(t+0.5*h)*(w+0.5*k2))+sin(3.0*(t+0.5*h)*(w+0.5*k2)));
      k4 = h*(cos(2.0*(t+h)*(w+k3))+sin(3.0*(t+h)*(w+k3)));
      w = w + (1.0/6.0)*(k1 + k2 + k3 + k4);
      fprintf(fptr4, "%f\t %f\n", t, w);
   }

   clock_t end = clock();
   double elapsed = (double)(end - start)/CLOCKS_PER_SEC;
    
   printf("Time measured: %.3f seconds.\n", elapsed);

   fclose(fptr1);
   fclose(fptr2);
   fclose(fptr3);
   fclose(fptr4);
}