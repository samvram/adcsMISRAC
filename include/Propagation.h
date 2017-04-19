#ifndef PROPAGATION_H_INCLUDED
#define PROPAGATION_H_INCLUDED

#include"Requisites.h"

#define PERIOD 6115.886623075233729 /**bound to be changed*/
int first_time;
void orbit_propagate( type time )
{
    type ddummy;
    int idummy;
    type we,w,t,thta;
    we = 7.292115854918357e-05;/* %rad/sec*/
    w = 2*M_PI/PERIOD;
    static type a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,a19,a20,a21,a22;
    static type b0,b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14,b15,b16,b17,b18,b19,b20,b21,b22;
    static type c0,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,c19,c20,c21,c22;

    if(first_time==0)
    {
        FILE* fp = fopen("fps_cof.dat","r"); /*Coefficient file name*/
        fscanf(fp,"%d%lf%lf%d%d%d%d%d%d%d%d",&idummy,&ddummy,&ddummy,&idummy,&idummy,&idummy,&idummy,&idummy,&idummy,&idummy,&idummy);
        fscanf(fp,"%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf",&a0,&a1,&a2,&a3,&a4,&a5,&a6,&a7,&a8,&a9,&a10,&a11);
/*        printf("a0=%e\n",a0);*/
        fscanf(fp,"%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf",&a12,&a13,&a14,&a15,&a16,&a17,&a18,&a19,&a20,&a21,&a22);
        fscanf(fp,"%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf",&b0,&b1,&b2,&b3,&b4,&b5,&b6,&b7,&b8,&b9,&b10,&b11);
        fscanf(fp,"%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf",&b12,&b13,&b14,&b15,&b16,&b17,&b18,&b19,&b20,&b21,&b22);
        fscanf(fp,"%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf",&c0,&c1,&c2,&c3,&c4,&c5,&c6,&c7,&c8,&c9,&c10,&c11);
        fscanf(fp,"%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf",&c12,&c13,&c14,&c15,&c16,&c17,&c18,&c19,&c20,&c21,&c22);
        fclose(fp);
/*        printf("b0=%e\nc0=%e\n",b0,c0);*/
        first_time=1;
    }
    t = (type)time;

    type swt = sin(w*t),cwt = cos(w*t), s2wet = sin(2.0*we*t), c2wet = cos(2*we*t);
    type swet = sin(we*t), cwet = cos(we*t), s2wt = sin(2.0*w*t), c2wt = cos(2.0*w*t);

    r_sat_eci[0] =  a0*1.0 +
                    a1*t +
                    a2*t*t +
                    a3*1.0*swt  +
                    a4 *t*swt +
                    a5*t*t*swt  +
                    a6*1.0*cwt  +
                    a7 *t*cwt +
                    a8*t*t*cwt  +
                    a9*1.0*swt*swt     +
                    a10 *t*swt*swt   +
                    a11*1.0*swt*cwt    +
                    a12 *t*swt*cwt   +
                    a13*(swt*swt*swt)             +
                    a14 * (swt*swt)*cwt  +
                    a15* s2wet                  +
                    a16 * c2wet                 +
                    a17 * swt*s2wet + a18 * swt*c2wet +
                    a19 * cwt*s2wet + a20 * cwt*c2wet +
                    a21 * swet + a22 * cwet;

    r_sat_eci[1] =  b0*1.0   +
                    b1 *t  +
                    b2*t*t   +
                    b3*1.0*swt  +
                    b4 *t*swt +
                    b5*t*t*swt    +
                    b6*1.0*cwt  +
                    b7 *t*cwt +
                    b8*t*t*cwt  +
                    b9*1.0*swt*swt     +
                    b10 *t*swt*swt   +
                    b11*1.0*swt*cwt    +
                    b12 *t*swt*cwt   +
                    b13*(swt*swt*swt)             +
                    b14 * (swt*swt)*cwt  +
                    b15* s2wet                  +
                    b16 * c2wet                   +
                    b17 * swt*s2wet        +
                    b18 * swt*c2wet        +
                    b19 * cwt*s2wet        +
                    b20 * cwt*c2wet        +
                    b21 * swet                     +
                    b22 * cwet;

    r_sat_eci[2] =  c0*1.0   +
                    c1 *t  +
                    c2*t*t   +
                    c3*1.0*swt  +
                    c4 *t*swt +
                    c5*t*t*swt  +
                    c6*1.0*cwt  +
                    c7 *t*cwt +
                    a8*t*t*cwt  +
                    c9*1.0*swt*swt     +
                    c10 *t*swt*swt   +
                    c11*1.0*swt*cwt    +
                    c12 *t*swt*cwt   +
                    c13*(swt*swt*swt)             +
                    c14 * (swt*swt)*cwt  +
                    c15* s2wet                    +
                    c16 * c2wet                   +
                    c17 * swt* s2wet       +
                    c18 * swt*c2wet        +
                    c19 * cwt* s2wet       +
                    c20 * cwt*c2wet        +
                    c21 * cwt* swet           +
                    c22 * cwt * cwet;

    v_sat_eci[0] =  a1* 1.0 +
                    a2*2.0*t +
                    a3* ( w * 1.0 * cwt) +
                    a4* ( 1.0*1.0*swt + w * t * cwt) +
                    a5* ( 2.0*t*swt + w * t*t * cwt) +
                    a6* ( - w * 1.0 * swt) +
                    a7* ( 1.0*1.0*cwt - w * t * swt) +
                    a8* ( 2.0*t*cwt - w * t*t * swt) +
                    a9* ( w * 1.0 * s2wt) +
                    a10* ( 1.0*1.0*swt*swt + w * t * s2wt) +
                    a11* ( w * 1.0 * c2wt) +
                    a12* ( 1.0*1.0*swt*cwt + w * t * c2wt) +
                    a13*3.0*w*swt*swt *cwt +
                    a14*w *( s2wt*cwt - (swt*swt*swt)) +
                    a15*we* c2wet -
                    a16*we*2.0*s2wet +
                    a17* (w*cwt*s2wet + 2.0*we*swt*c2wet ) +
                    a18* (w*cwt*c2wet - 2.0*we*swt*s2wet ) +
                    a19* (2.0*we*cwt*c2wet - w*swt*s2wet ) -
                    a20* (2.0*we*cwt*c2wet + w*swt*c2wet ) +
                    a21* (we*cwt*cwet- w*swt*swet )-
                    a22*(we*cwt*swet+w*swt*cwet);

    v_sat_eci[1] =  b1* 1.0 +
                    b2*2.0*t +
                    b3* ( w * 1.0 * cwt) +
                    b4* ( 1.0*1.0*swt + w * t * cwt) +
                    b5* ( 2.0*t*swt + w * t*t * cwt) +
                    b6* ( - w * 1.0 * swt) +
                    b7* ( 1.0*1.0*cwt - w * t * swt) +
                    b8* ( 2.0*t*cwt - w * t*t * swt) +
                    b9* ( + w * 1.0 * s2wt) +
                    b10* ( 1.0*1.0*swt*swt + w * t * s2wt) +
                    b11* ( + w * 1.0 * c2wt) +
                    b12* ( 1.0*1.0*swt*cwt + w * t * c2wt) +
                    b13*3.0*w*swt*swt *cwt +
                    b14*w *( s2wt*cwt - (swt*swt*swt)) +
                    b15*we* c2wet -
                    a16*2.0*we* s2wet +
                    b17* (w*cwt*s2wet + 2.0*we*swt*c2wet ) +
                    b18* (w*cwt*c2wet - 2.0*we*swt*s2wet ) +
                    b19* (2.0*we*cwt*c2wet - w*swt*s2wet ) -
                    b20* (we*cwt*c2wet + w*swt*c2wet ) +
                    b21* (we*cwt*cwet- w*swt*swet )-
                    b22*(we*cwt*swet+w*swt*cwet);

    v_sat_eci[2] =  c1* 1.0 +
                    c2*2.0*t +
                    c3* ( w * 1.0 * cwt) +
                    c4* ( 1.0*1.0*swt + w * t * cwt) +
                    c5* ( 2.0*t*swt + w * t*t * cwt) +
                    c6* ( - w * 1.0 * swt) +
                    c7* ( 1.0*1.0*cwt - w * t * swt) +
                    c8* ( 2.0*t*cwt - w * t*t * swt) +
                    c9* ( + w * 1.0 * s2wt) +
                    c10* ( 1.0*1.0*swt*swt + w * t * s2wt) +
                    c11* ( w * 1.0 * c2wt) +
                    c12* ( 1.0*1.0*swt*cwt + w * t * c2wt) +
                    c13*3.0*w*swt*swt *cwt +
                    c14*w *( s2wt*cwt - (swt*swt*swt)) +
                    c15*we* c2wet - 2.0 * a16*we* s2wet +
                    c17* (w*cwt*s2wet + 2.0*we*swt*c2wet ) +
                    c18* (w*cwt*c2wet - 2.0*we*swt*s2wet ) +
                    c19* (2.0*we*cwt*c2wet - w*swt*s2wet ) -
                    c20* (we*cwt*c2wet + w*swt*c2wet ) +
                    c21*(we*cwt*cwet - w*swt* swet) -
                    c22*(we*cwt*swet + w*swt* cwet);
    r_sat_eci = arrscalmult(1e3,r_sat_eci,3,F);
    v_sat_eci = arrscalmult(1e3,v_sat_eci,3,F);

    free(w_orbit);
    w_orbit = arrscalmult(1.0/(norm(r_sat_eci,3,K)*norm(v_sat_eci,3,K)),cross(r_sat_eci,v_sat_eci,K,K),3,F);
    thta = acos(dot(w_orbit,r_sat_eci)/norm(r_sat_eci,3,K));   /*Angle between w_orbit and r_sat_eci*/
    w_orbit = arrscalmult(norm(v_sat_eci,3,K)/(norm(r_sat_eci,3,K)*sin(thta)),w_orbit,3,F);
}

/*type Eanomalyfn( type E,type e,type M )
{
    return E - e*sin(E) - M;
}

type Eanomalyfn_( type E,type e,type M )
{
    return 1 - e*cos(E);
}

type gettrueanomaly( type n,type t,type e)
{
    type E = n*t;
    type allowance = 0.001;
    do
    {
        E = E - Eanomalyfn( E,e,n*t )/Eanomalyfn_( E,e,n*t );
    }while( pow(Eanomalyfn(E,e,n*t),2)>pow( allowance,2) );

    type theta = 2*atan( sqrt( (1+e)/(1-e) )*tan(E/2) );

    return theta;
}

void orbit_propagate( type t )//parameters are string satellite name string and Julian time
Global Parameters modified:
r_sat_eci
w_orbit
v_sat_eci

{
    #define G 6.67*pow(10,-11)
    #define Me 5.972*pow(10,24)

    type mu = G*( Me+Msat );
*/
    /*FROM ENVISAT TLE USED BY ETIKA AGARWAL------------------------------------------------------------------*/
/*    type n = 14.32241750*2*M_PI/(23.59*60*60);   //Mean motion, from revs per day to radians per second
    type a = pow( mu/pow(n,2),1.0/3.0 );         //Semi-Major Axis
    type argp = 074.8897*M_PI/180;      //Argument of Perigee
    type RAAN = 212.1479*M_PI/180;      //Right Ascension of Ascending Node
    type i = 98.5445*M_PI/180;         //Inclination
    type e = 0.0001011;         //Eccentricity
    type theta = gettrueanomaly( n,t,e );     //True anomaly*/
    /*--------------------------------------------------------------------------------------------------------*/
/*
    r_sat_eci[0] = a*(1-e*e)/(1-e*cos(theta));
    r_sat_eci[1] = M_PI/2;
    r_sat_eci[2] = theta;                       //The planar position of the satellite on the xy plane in rtp

    r_sat_eci = rtptoxyz( r_sat_eci,F );

    type *qargp = calloc( 4,sizeof(type) ), *qinc = calloc( 4,sizeof(type) ) ,*qRAAN = calloc( 4,sizeof(type) ), *qtotal;
    qargp[0] = cos(-argp/2);
    qargp[1] = 0;
    qargp[2] = 0;
    qargp[3] = sin(-argp/2);

    qinc[0] = cos(-i/2);
    qinc[1] = sin(-i/2);
    qinc[2] = 0;
    qinc[3] = 0;

    qRAAN[0] = cos(-RAAN/2);
    qRAAN[1] = 0;
    qRAAN[2] = 0;
    qRAAN[3] = sin(-RAAN/2);

    qtotal = qmult( qargp, qmult(qinc, qRAAN,F,F ),F,F );

    type *perpaxis = calloc( 3,sizeof(type) );
    perpaxis[0] = 0;
    perpaxis[1] = 0;
    perpaxis[2] = 1;

    r_sat_eci = qtransform( qtotal,r_sat_eci,K,F );

    perpaxis = qtransform( qtotal,perpaxis,F,F );

    type Vnorm = sqrt( mu*(2/norm(r_sat_eci,3,K)-1/a) );
    type wnorm = Vnorm/norm( r_sat_eci,3,K );

    free(w_orbit);
    w_orbit = arrscalmult( wnorm,perpaxis,3,K );

    free(perpaxis);
}*/

void refgen()
{
    getq_ned_eci( r_sat_eci );
    q_drf_eci = realloc( q_drf_eci,4*sizeof(type) );
    /*    q_drf_eci[0] = q_ned_eci[0];
    q_drf_eci[1] = q_ned_eci[1];
    q_drf_eci[2] = q_ned_eci[2];
    q_drf_eci[3] = q_ned_eci[3];*/
    q_drf_eci[0] = 1;
    q_drf_eci[1] = 0;
    q_drf_eci[2] = 0;
    q_drf_eci[3] = 0;
    /*    type *xdrf_eci,*ydrf_eci,*zdrf_eci;
    xdrf_eci = arrscalmult(1.0/norm(v_sat_eci,3,K),v_sat_eci,3,K);
    ydrf_eci = arrscalmult(-1.0/norm(r_sat_eci,3,K),r_sat_eci,3,K);
    zdrf_eci = cross(xdrf_eci,ydrf_eci,K,K);
    zdrf_eci = arrscalmult(1.0/norm(zdrf_eci,3,K),zdrf_eci,3,F);

    type **Tmx = malloc(3*sizeof(type*));
    Tmx[0] = malloc(3*sizeof(type));
    Tmx[1] = malloc(3*sizeof(type));
    Tmx[2] = malloc(3*sizeof(type));    
    Tmx[0][0] = xdrf_eci[0]; Tmx[0][1] = xdrf_eci[1]; Tmx[0][2] = xdrf_eci[2];
    Tmx[1][0] = ydrf_eci[0]; Tmx[1][1] = ydrf_eci[1]; Tmx[1][2] = ydrf_eci[2];
    Tmx[2][0] = zdrf_eci[0]; Tmx[2][1] = zdrf_eci[1]; Tmx[2][2] = zdrf_eci[2];    

    q_drf_eci = r2q(Tmx);

    free(xdrf_eci);
    free(ydrf_eci);
    free(zdrf_eci);
    free(Tmx[0]);
    free(Tmx[1]);
    free(Tmx[2]);
    free(Tmx);*/

    w_drf_eci = realloc( w_drf_eci,3*sizeof(type) );
    /*    w_drf_eci[0] = -w_orbit[0];
    w_drf_eci[1] = -w_orbit[1];
    w_drf_eci[2] = -w_orbit[2];
    w_drf_eci = qtransform(q_drf_eci,w_drf_eci,K,F);*/
    w_drf_eci[0] = 0;
    w_drf_eci[1] = 0;
    w_drf_eci[2] = 0;
}
#endif /* PROPAGATION_H_INCLUDED */
