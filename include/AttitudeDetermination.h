#ifndef ATTITUDEDETERMINATION_H_INCLUDED
#define ATTITUDEDETERMINATION_H_INCLUDED

#include"Requisites.h"
#include"Variables.h"

/*QUEST-----------------------*/
type* quest(type** V_BRF, type** V_ECI, int n)                   /*Arrays with n activated sun-sensors.*/
{                                                               /*The last row, row 6 is reserved for magnetometer*/
    int i;          						   /*OUTPUTS quaternion OF Reference w.r.t BRF*/
    type weight,sigma,del,kap,alpha,beta,gamma,w,NORM,**B,**S,**X,**temp1,*Z,*X1,*M;
    type temp2,temp,*Q_opt;

    for( i=0;i<7;i++ )
    {
      NORM = norm(V_BRF[i],3,K);
/*        if( norm(V_BRF[i],3,K)!=0 )*/
      if( NORM!=0 )
        {
            /*V_BRF[i] = arrscalmult( 1/norm(V_BRF[i],3,K),V_BRF[i],3,F );*/
	    V_BRF[i] = arrscalmult(1.0/NORM,V_BRF[i],3,F);
            V_ECI[i] = arrscalmult( 1/norm(V_ECI[i],3,K),V_ECI[i],3,F );
        }
    }
    /*Sun sensors are roughly 10 times more accurate than magnetometers. So the weight for sun sensor is 10 times that of
    Magnetometer.*/

    weight = 1.0/(10*n+1);
    B = scalar_mult_mat(10*weight,vector_mult( V_BRF[0],V_ECI[0] ),F );
    for( i=1;i<n;i++ )
        B = matrix_add( B,scalar_mult_mat(10*weight,vector_mult( V_BRF[i],V_ECI[i] ),F ),F,F );
    B = matrix_add( B,scalar_mult_mat(weight,vector_mult( V_BRF[6],V_ECI[6]),F),F,F );
    /* B = summation(weights(i) * vector_from_sensor(body_ref_frame)(i) * reference_vector(ECI_frame)(i)*/
    sigma = trace(B,K);
    /* sigma = trace of B*/
    S = matrix_add( B,transpose(B),F,F );
    /* S = B + B_transpose*/
    del = determinant(S);
    /* delta = determinant of S*/
    kap = trace( scalar_mult_mat( determinant(S),mxinv(S),F),F );/*adjoint(S));*/
    /* kappa = trace of adjoint of S*/
    w = 1; /* sum of weights is 1*/
    alpha = (w*w) - (sigma*sigma) + kap ;
    beta = w - sigma;                   /* from Etika's report (matlab file ; quest.m)*/
    gamma = (w+sigma)*alpha - del;

    X = matrix_add(scalar_mult_mat(alpha,eye,K), 
		scalar_mult_mat(beta,S,K),F,F); /* alpha_I + beta_S*/
    temp1 = matrix_multiplication(S,S,3,3,3,3,K,F); /* S_squar*/
    X = matrix_add(X,temp1,F,F); /* alpha_I + beta_S + S_square*/

    Z = arrscalmult( 10*weight,cross(V_ECI[0],V_BRF[0],K,K),3,F );
    for( i=1;i<n;i++ )
        Z = arrsum( Z,arrscalmult( 10*weight,cross(V_ECI[i],V_BRF[i],K,K),3,F ),3,F,F );
    Z = arrsum( Z,arrscalmult(weight,cross(V_ECI[6],V_BRF[6],K,K),3,F),3,F,F);
/* 0.5*(w_sun x v_sun) + 0.5*(w_mag x v_mag)*/

    X1 = matmult(X,Z,F,F); /* required X*/
    M = malloc(4*sizeof(type));      /* M is penultimate optimal quaternion*/
    if(gamma > 0)
    {
        M[0] = gamma;
        M[1] = X1[0];
        M[2] = X1[1];
        M[3] = X1[2];
    }
    else
    {
        M[0] = -gamma;
        M[1] = -X1[0];
        M[2] = -X1[1];
        M[3] = -X1[2];
    }

    temp2 = X1[0]*X1[0] + X1[1]*X1[1] + X1[2]*X1[2] ; /* X1_transpose * X1*/
    free(X1);

    temp = sqrt((gamma*gamma) + (temp2));
    Q_opt = arrscalmult( 1.0/temp,M,4,F );
    Q_opt = qinv(Q_opt,F);/* Required Quaternion; of brf w.r.t eci*/

    if( Q_opt[0]<0 )
        Q_opt = arrscalmult(-1,Q_opt,4,F);

    for( i=0;i<7;i++ )
    {
      free(V_BRF[i]);
      free(V_ECI[i]);
    }
    free(V_BRF);
    free(V_ECI);

    return Q_opt;
}

type *SensorFusion( int ons ) /* Input argument is no. of sunsensors
				which are on */
{
  int i,j;
  type **ECI;
  type **BRF;
  ECI = malloc(7*sizeof(type));
  BRF = malloc(7*sizeof(type));
  for( i=0;i<7;i++ )
  {
    ECI[i] = malloc(3*sizeof(type));
    BRF[i] = malloc(3*sizeof(type));
  }
  for( i=0;i<6;i++ )
  {
    if( i<ons )
    {
      ECI[i][0] = sun_ECI[0];
      ECI[i][1] = sun_ECI[1];
      ECI[i][2] = sun_ECI[2];
      
      BRF[i][0] = sun_brf[0];
      BRF[i][1] = sun_brf[1];
      BRF[i][2] = sun_brf[2];
    }
    else
      for( j=0;j<3;j++ )
      {
        ECI[i][j] = 0;
        BRF[i][j] = 0;
      }
  }
      ECI[i][0] = B_eci[0];
      ECI[i][1] = B_eci[1];
      ECI[i][2] = B_eci[2];
      
      BRF[i][0] = B_brf[0];
      BRF[i][1] = B_brf[1];
      BRF[i][2] = B_brf[2];

  return quest(BRF,ECI,ons);
}

type *RateDetermination()
{ 
  int i; 
  type *dqbdt,*Qw,*W;
  dqbdt = malloc(4*sizeof(type));
    dqbdt[0] = (qd_brf_eci[0]-qd0_brf_eci[0])/step;
    dqbdt[1] = (qd_brf_eci[1]-qd0_brf_eci[1])/step;
    dqbdt[2] = (qd_brf_eci[2]-qd0_brf_eci[2])/step;
    dqbdt[3] = (qd_brf_eci[3]-qd0_brf_eci[3])/step;
  Qw = qmult(qinv(qd_brf_eci,K),arrscalmult(2,dqbdt,4,F),F,F);

  W = malloc(3*sizeof(type));
  for( i=1;i<4;i++ )
    W[i-1] = Qw[i];
  free(Qw);
  free(wd_brf_eci);

  return W;
}

/* KALMAN FILTER */
type *kalman_filter( type* Z0,type* X0,FREE FZ,FREE FX )/* Z0 is calculation; X0 is measurement */
{
    int i,j,k;
    type t_sample = 0.1 ;
    type *X,*Z,**M,*Y,**S_k,**KG,**P,**F_k,**Q,**phi,**R,**H_k;
    X = malloc(7*sizeof(type));
      X[0] = X0[0];
      X[1] = X0[1];
      X[2] = X0[2];
      X[3] = X0[3];
      X[4] = X0[4];
      X[5] = X0[5];
      X[6] = X0[6];
    
    Z = (type*)malloc(7*sizeof(type));

    P = (type**)malloc(7*sizeof(type*))  ;     /* covariance matrix*/
    F_k = (type**)malloc(7*sizeof(type*));
    Q = (type**)malloc(7*sizeof(type*))  ;     /* Q matrix*/
    phi =  (type**)malloc(7*sizeof(type*))  ;     /* phi matrix*/
    R = (type**)malloc(7*sizeof(type*)) ;      /* R matrix*/
    H_k = (type**)malloc(7*sizeof(type*));
    for(i=0;i<7;i++)
    {
        F_k[i] = (type*)malloc(7*sizeof(type));
        phi[i] = (type*)malloc(7*sizeof(type));
        P[i] = (type*)calloc(7,sizeof(type));
        P[i][i] = 1;
        Q[i] = (type*)calloc(7,sizeof(type));
        Q[i][i] = 0.05;
        R[i] = (type*)calloc(7,sizeof(type));
        R[i][i] = 0.0016;
        H_k[i] = (type*)calloc(7,sizeof(type));
        H_k[i][i] = 1;
    }

    F_k[0][0] = 0;                  F_k[0][1] = -I[2][2]/I[0][0];     F_k[0][2] = I[1][1]/I[0][0];      F_k[0][3] = 0;      F_k[0][4] = 0;        F_k[0][5] = 0;      F_k[0][6] = 0;
    F_k[1][0] = I[0][0]/I[1][1];    F_k[1][1] = 0;                    F_k[1][2] = I[2][2]/I[1][1];      F_k[1][3] = 0;      F_k[1][4] = 0;        F_k[1][5] = 0;      F_k[1][6] = 0;
    F_k[2][0] = I[1][1]/I[2][2];    F_k[2][1] = I[0][0]/I[2][2];      F_k[2][2] = 0;                    F_k[2][3] = 0;      F_k[2][4] = 0;        F_k[2][5] = 0;      F_k[2][6] = 0;
    F_k[3][0] = 0;                  F_k[3][1] = 0;                    F_k[3][2] = 0;                    F_k[3][3] = 0;      F_k[3][4] = 0.5;      F_k[3][5] = -0.5;   F_k[3][6] = 0.5;
    F_k[4][0] = 0;                  F_k[4][1] = 0;                    F_k[4][2] = 0;                    F_k[4][3] = -0.5;   F_k[4][4] = 0;        F_k[4][5] = 0.5;    F_k[4][6] = 0.5;
    F_k[5][0] = 0;                  F_k[5][1] = 0;                    F_k[5][2] = 0;                    F_k[5][3] = 0.5;    F_k[5][4] = -0.5;     F_k[5][5] = 0;      F_k[5][6] = 0.5;
    F_k[6][0] = 0;                  F_k[6][1] = 0;                    F_k[6][2] = 0;                    F_k[6][3] = -0.5;   F_k[6][4] = -0.5;     F_k[6][5] = -0.5;   F_k[6][6] = 0;

    for(i=0;i<7;i++)
        for(j=0;j<7;j++)
            phi[i][j] = eye7[i][j] + F_k[i][j]*t_sample;

    for(k=0;k<5;k++)
    {
        M = sum_77( matrix_multiplication( phi, 
			matrix_multiplication(P,transpose_7(phi,K),7,7,7,7,F,F),7,7,7,7,K,F ),Q,F,K);
        X = matmult_7(phi,X,K,F);
        for(j=0;j<7;j++)
            Z[j] = Z0[j];/*rand()%0.002;*/

        Y =  sum_17( Z,scal_17( -1,matmult_7(H_k,X,K,K),F ),K,F );
        S_k = sum_77( matrix_multiplication( H_k,
 	        matrix_multiplication( M,transpose_7(H_k,K),7,7,7,7,K,F ),7,7,7,7,K,F ),R,F,K );
        KG = matrix_multiplication( M,matrix_multiplication( transpose_7(H_k,K),inverse_7(S_k,F),7,7,7,7,F,F ), 7,7,7,7,K,F );
        X = sum_17( X,matmult_7( KG,Y,K,F ),F,F );
        P = matrix_multiplication( sum_77( eye7, scal_77( -1, matrix_multiplication(KG,H_k,7,7,7,7,F,K),F ),K,F ),M,7,7,7,7,F,F );
    }

    for(i=0;i<7;i++)
    {
        free(Q[i]);
        free(P[i]);
        free(F_k[i]);
        free(phi[i]);
        free(R[i]);
        free(H_k[i]);
    }
    free(Q);
    free(P);
    free(H_k);
    free(phi);
    free(R);
    free(F_k);
    free(Z);

    if( FZ==F )
      free(Z0);
    if( FX==F )
      free(X0);

    return X;
}

void Filter()
{
  type *SM,*SC,*SF;
  SM = malloc(7*sizeof(type));
  SC = malloc(7*sizeof(type));
	/* State Measurement, Calculation & Filter */
  SC[0] = w_brf_eci[0];
  SC[1] = w_brf_eci[1];
  SC[2] = w_brf_eci[2];
  SC[3] = q_brf_eci[1];
  SC[4] = q_brf_eci[2];
  SC[5] = q_brf_eci[3];
  SC[6] = q_brf_eci[0];

  SM[0] = wd_brf_eci[0];
  SM[1] = wd_brf_eci[1];
  SM[2] = wd_brf_eci[2];
  SM[3] = qd_brf_eci[1];
  SM[4] = qd_brf_eci[2];
  SM[5] = qd_brf_eci[3];
  SM[6] = qd_brf_eci[0];

  SF = kalman_filter(SC,SM,F,F);
/* UPDATION */
  wd_brf_eci[0] = SF[0];
  wd_brf_eci[1] = SF[1];
  wd_brf_eci[2] = SF[2];

  qd_brf_eci[0] = SF[6];
  qd_brf_eci[1] = SF[3];
  qd_brf_eci[2] = SF[4];
  qd_brf_eci[3] = SF[5];
  
  free(SF);
}

#endif /* ATTITUDEDETERMINATION_H_INCLUDED*/
