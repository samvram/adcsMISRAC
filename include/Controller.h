#ifndef CONTROLLER_H_INCLUDED
#define CONTROLLER_H_INCLUDED

#include"Requisites.h"
#include"Variables.h"

void computeslidingvector()
{
    free(S);
    S = arrsum( matmult( I,wd_brf_drf,K,K ),matmult( LAMBDA_q,qd_brf_drf+1,K,K ),3,F,F );
}

void computeLAMBDA_s()
{
    int i=0,j;
    type **Nmx;

    Nmx = malloc( 3*sizeof(type*) );
    for( i=0;i<3;i++ )
    {
        Nmx[i] = calloc( 3,sizeof(type) );
        for( j=0;j<D;j++ )
        {
            Nmx[i][j] = 0;
            Nmx[i][i] += N[j][i];
        }
    }

    for( i=0;i<3;i++ )
        free(LAMBDA_s[i]);
    free(LAMBDA_s);

    LAMBDA_s = scalar_mult_mat(1.0/norm(S,3,K),matrix_add( eta,Nmx,K,K ),F );

    for( i=0;i<3;i++ )
        free(Nmx[i]);
	free(Nmx);

	/*    LAMBDA_s[0][0] = 0.003; LAMBDA_s[0][1] = 0.000; LAMBDA_s[0][1] = 0.000;
    LAMBDA_s[1][0] = 0.000; LAMBDA_s[1][1] = 0.003; LAMBDA_s[1][1] = 0.000;
    LAMBDA_s[2][0] = 0.000; LAMBDA_s[2][1] = 0.000; LAMBDA_s[2][1] = 0.003;*/
}

type* sliding_controller()
{
    int i,j;
    type *dqbdt,*ucap,*ncaps,*U;

    dqbdt = arrsum( arrscalmult( qd_brf_drf[0],wd_brf_drf,3,K ),cross( qd_brf_drf+1,wd_brf_drf,K,K ),3,F,F );
    dqbdt = matmult( LAMBDA_q,dqbdt,K,F );
    dqbdt = arrscalmult( 0.5,dqbdt,3,F );

    ucap = arrsum( cross( w_brf_eci,matmult( I,w_brf_eci,K,K ),K,F ),matmult( I,alpha_drf_eci,K,K ),3,F,F );
    ucap = arrsum( ucap,arrscalmult(-1,dqbdt,3,F),3,F,F );

    ncaps = malloc( 3*sizeof(type) );         /*Compensation for Disturbances*/
    for( i=0;i<3;i++ )
    {
      ncaps[i] = 0;
      for( j=0;j<D;j++ )
        ncaps[i] += ncap[j][i];
    }

    ucap = arrsum( ucap,arrscalmult(-1,ncaps,3,F),3,F,F );				/*Equivalent Control; */
    computeslidingvector();
    computeLAMBDA_s();

    U = arrsum( ucap,arrscalmult(-1,matmult( LAMBDA_s,S,K,K ),3,F ),3,F,F );	/*Control Law*/
    /*return U;*/
    /*Geometrical Approach------------------------------------------------------------*/
    /*type* ncap = cross( B_brf,cross( B_brf,U ) );
    ncap = arrscalmult( 1.0/norm(ncap,3),ncap,3 );

    if( dot(ncap,U)<0 )
        ncap = arrscalmult( -1,ncap,3 );

    type normTctrl = pow( norm(U,3),2 )/dot(ncap,U);

    type* Bcap = arrscalmult( 1.0/norm(B_brf,3),B_brf,3 );
    type normM = normTctrl/norm(B_brf,3);

    type* Mreqd = arrscalmult(normM,cross(Bcap,ncap),3);
    printf("(%lf,%lf,%lf)\n",Mreqd[0],Mreqd[1],Mreqd[2]);*/
    /*---------------------------------------------------------------------------------------*/

    /*Wisnewski's practice-------------------------------------------------------------------*/

    type *scap,*N_des_par,*Mom;
    scap = arrscalmult( 1/norm(S,3,K),S,3,K );
    N_des_par = arrscalmult( dot( U,scap ),scap,3,F );
    free(U);
    sang = acos(dot( scap,B_brf )/norm(B_brf,3,K))*180./M_PI;
    Mom = arrscalmult( 1.0/pow( norm(B_brf,3,K),2 ),cross( B_brf,N_des_par,K,F ),3,F );
    /*---------------------------------------------------------------------------------------*/
    /* Saturation */
    for( i=0;i<3;i++ )
//	  Mom[i] = (fabs(Mom[i])>0.2)?0.2*Mom[i]/fabs(Mom[i]):Mom[i];
    Mom[i] = (fabs(Mom[i])>10.0)?10.0*Mom[i]/fabs(Mom[i]):Mom[i];
    return Mom;
}

type* DetumblingMoment()
{
    int i;
    type** k;

    k = (type**)malloc( 3*sizeof(type*) );
    for( i=0;i<3;i++ )
        k[i] = (type*)malloc( 3*sizeof(type) );

    k[0][0] = -pow(10,9);    k[0][1] = 0;            k[0][2] = 0;
    k[1][0] = 0;            k[1][1] = -pow(10,9);    k[1][2] = 0;
    k[2][0] = 0;            k[2][1] = 0;            k[2][2] = -pow(10,9);

    /*type* M = matmult( k,B_dot,F,K );

    for( i=0;i<3;i++ )
        if(M[i]>5)
            M = arrscalmult( 5/M[i],M,3,F );

    return M;*/

    return matmult( k,B_dot,F,K );
}

#endif /* CONTROLLER_H_INCLUDED*/
