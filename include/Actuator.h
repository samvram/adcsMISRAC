#ifndef ACTUATOR_DEFD
#define ACTUATOR_DEFD

#include"Requisites.h"
#include"Variables.h"

type u(type t)
{
  type ret;
  if( t<0 )
    ret = 0;
  else
    ret = 1;
  return ret;
}

type I_resp(type t,type toff,int t_ind)
{
  type ret = Vo*(u(t)*(1-exp(-T[t_ind].R*t/T[t_ind].L))-
	u(t-toff)*(1-exp(-T[t_ind].R*(t-toff)/T[t_ind].L)))/T[t_ind].R;
  ret = (ret>T[t_ind].Imax)?T[t_ind].Imax:ret;

  return ret;
}

type I_integ(type toff,int t_ind) /* t_ind - torquer index, 0, 1 or 2 */
{
  return Vo*(toff+T[t_ind].L*exp(-T[t_ind].R*ta/T[t_ind].L)*(1-exp(-T[t_ind].R*toff/T[t_ind].L))/T[t_ind].R)/T[t_ind].R;
}
//#define integfunc(toff) Vo*(toff+L*exp(-R*ta/L)*(1-exp(-R*toff/L))/R)/R

type get_toff(type Iavg,int t_ind)
{
  type toff = 0.1,F,dF;
  #define h 0.001
  #define e 0.001
  F = I_integ(toff,t_ind);
  dF = (I_integ(toff+h,t_ind)-I_integ(toff-h,t_ind))/(2*h);
  do
  {
    toff = toff-F/dF;
    
    F = I_integ(toff,t_ind)-Iavg*ta;
    dF = (I_integ(toff+h,t_ind)-I_integ(toff-h,t_ind))/(2*h);
  }while( pow(F,2)>pow(e,2) );
  toff = toff - F/dF;

  return toff;
}

type *Actuation( type *Mreqd,FREE FR )
{
  int i,j;
  type *Mcom = (type*)malloc(3*sizeof(type));
  
  type **Ax = (type**)malloc(3*sizeof(type*));
  for( i=0;i<3;i++ )
    Ax[i] = (type*)malloc(3*sizeof(type));

  for( i=0;i<3;i++ )
    for( j=0;j<3;j++ )
      Ax[i][j] = T[j].Axis[i];
 
  Mcom = matmult(mxinv(Ax),Mreqd,F,K);
	/* Mcom is a vector whose components
		are the magnitudes of moment
		required from each torquer */
  type Iavg,toff;
  int usign[3]; /* unit sign (+1/-1) of signal */
  int n[3]; /* number of pwm cycles that are to be on */
  for( i=0;i<3;i++ )
  {
    Iavg = Mcom[i]/(T[i].N*T[i].A);
    usign[i] = sqrt( Iavg*Iavg )/Iavg;
    Iavg = usign[i]*Iavg;
    toff = get_toff(Iavg,i);
    /* Saturation */
    toff = (toff>ta)?ta:toff;
    n[i] = toff/tp;
  }

  type* Mprod = (type*)calloc(3,sizeof(type));
  type simt; /* time for torquer simulation */
  for( i=0;i<3;i++ )
  {
    for( simt=0;simt<ta;simt+=tp )
      Mprod[i] += T[i].N*I_resp(simt,n[i]*tp,i)*T[i].A*usign[i];
    Mprod[i] = Mprod[i]/ta;
  }
  Mprod = matmult(Ax,Mprod,F,F);
  free(Mcom);
  if( FR==F )
    free(Mreqd);

  return Mprod;
}

#endif
