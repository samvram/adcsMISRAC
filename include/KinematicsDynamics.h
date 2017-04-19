#ifndef KINEMATICSDYNAMICS_H_INCLUDED
#define KINEMATICSDYNAMICS_H_INCLUDED

#include"Requisites.h"
#include"Variables.h"

type** Disturbance_torques()
{
    int i,j;
    type **ret;

    ret = malloc( D*sizeof(type) );
    for( i=0;i<D;i++ )
        ret[i] = malloc( 3*sizeof(type) );

    /*ZERO DISTURBANCE STUB-ZERO DISTURBANCE STUB-*/
    for( i=0;i<D;i++ )
        for( j=0;j<3;j++ )
            ret[i][j] = 0;
    /*ZERO DISTURBANCE STUB-ZERO DISTURBANCE STUB-*/

    return ret;
}

type** state_fn( type** State,type *Tctrl, FREE FR )
{
    int i;
    type **ret,*H,*w,*dHbydt,**n;

    ret = (type**)malloc(2*sizeof(type*));
    ret[0] = arrscalmult(0.5,qmult(State[0],State[1],K,K),4,F);   /*dqbydt*/
    H = matmult(I,State[1]+1,K,K);
    w = (type*)malloc(3*sizeof(type));
    w[0] = State[1][1];
    w[1] = State[1][2];
    w[2] = State[1][3];

    dHbydt =  arrsum(arrscalmult(-1,cross(w,H,F,F),3,F),Tctrl,3,F,K);
    n = Disturbance_torques();
    for( i=0;i<D;i++ )
        dHbydt = arrsum(dHbydt,n[i],3,F,K);
    dHbydt = matmult(Iinv,dHbydt,K,F);

    /*ret[1] = dHbydt-1;*/
    ret[1] = (type*)malloc( 4*sizeof(type) );
    ret[1][0] = 0;
    ret[1][1] = dHbydt[0];
    ret[1][2] = dHbydt[1];
    ret[1][3] = dHbydt[2];

    free(dHbydt);
    for( i=0;i<D;i++ )
        free(n[i]);
    free(n);
    if( FR==F )
    {
        for( i=0;i<2;i++ )
            free(State[i]);
        free(State);
    }
    
    return ret;
}

type** Runge_Kutta( type** Y,type* Tctrl,type h )
{
    type **K1,**K2,**K3,**K4;

    K1 = scal_24( h,state_fn(Y,Tctrl,K),F );
    K2 = scal_24( h,state_fn(sum_24(Y,scal_24(0.5,K1,K),K,F),Tctrl,F),F );
    K3 = scal_24( h,state_fn(sum_24(Y,scal_24(0.5,K2,K),K,F),Tctrl,F),F );
    K4 = scal_24( h,state_fn(sum_24(Y,K3,K,K),Tctrl,F),F );

    return sum_24(Y,scal_24(1.0/6.0,sum_24(sum_24(K1,K4,F,F),scal_24(2,sum_24(K2,K3,F,F),F),F,F),F),F,F);
}

void systemresponse( type* Tctrl,type t_step )
{
    int i;
    type **state;

    state = (type**)malloc( 2*sizeof(type*) );
    for( i=0;i<2;i++ )
        state[i] = (type*)malloc( 4*sizeof(type) );

    for( i=0;i<4;i++ )
        state[0][i] = q_brf_eci[i];

    state[1][0] = 0;
    for( i=1;i<4;i++)
        state[1][i] = w_brf_eci[i-1];
    state = Runge_Kutta(state,Tctrl,t_step);

    for( i=0;i<4;i++ )
        q_brf_eci[i] = state[0][i];

    state[1][0] = 0;
    for( i=1;i<4;i++)
        w_brf_eci[i-1] = state[1][i];

    q_brf_eci = arrscalmult( 1.0/norm(q_brf_eci,4,K),q_brf_eci,4,F );

    for( i=0;i<2;i++ )
        free(state[i]);
    free(state);

    if( q_brf_eci[0]<0 )
      q_brf_eci = arrscalmult(-1.0,q_brf_eci,4,F );
}

#endif /* KINEMATICSDYNAMICS_H_INCLUDED */
