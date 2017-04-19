#ifndef VARIABLES_H_INCLUDED
#define VARIABLES_H_INCLUDED

typedef double type;

/**NOTATIONS FOLLOWED-------------------------------------------------------**/
/**A quaternion Q_A_B is defined as quaternion OF CS A w.r.t CS B. This     **/
/**quaternion may be used to transform a vector in CS B to CS A by using the**/
/**function call qrotate( Q_A_B,V_B ).                                      **/
/**COORDINATE SYSTEMS-------------------------------------------------------**/
/**1) ECI - Earth Centered Inertial; Centered at the center of mass of the  **/
/**   Earth, with x axis along the vernal equinox, and z axis along the     **/
/**   axis of earth's rotation in the Geometric North direction.            **/
/**2) ECEF - Earth Centered Earth Fixed; Centered at the center of mass of  **/
/**   the Earth, with x axis along the Greenwich Meridian on the equatorial **/
/**   plane, and z axis along the axis of earth's rotation in the Geometric **/
/**   North direction.                                                      **/
/**3) BRF - Body Reference Frame; Centered at the center of mass of the     **/
/**   spacecraft, with the x axis along the direction of maximum principal  **/
/**   moment of Inertia, and z axis along the direction of minimum principal**/
/**   moment of Inertia.                                                    **/
/**4) NED - North-East-Down coordinate system; Centered at the center of    **/
/**   mass of the spacecraft, with the x axis pointing towards the earth's  **/
/**   North Pole and the z axis along the vector from the satellite's center**/
/**   of gravity to the Earth's center of mass.                             **/
/**5) DRF - Desired Reference Frame; Centered at the center of mass of the  **/
/**   spacecraft, with the x,y and z directions oriented along the desired  **/
/**   directions.                                                           **/
/**6) GBRF - Geometric Body Reference Frame; Centered at the geometric      **/
/**   center of the spacecraft with the x, y and z oriented as in BRF.      **/
/**7) SRFi - Sensor Reference frame of i^{th} sensor; Design coordinate     **/
/**   frame of the sensor, having a constant attitude with respect to BRF.  **/
/**8) ARFi - Actuator Reference frame of i^{th} actuator; Design coordinate **/
/**   frame of the actuator, having a constant attitude with respect to BRF.**/
/**9) ALL - Altitude Latitude Longitude Coordinate system - This is one way **/
/**   of expressing the ECEF coordinate system. Altitude is the radial      **/
/**   distance of the satellite from the center of the Earth. Latitude is   **/
/**   pi/2 - theta and Longitude is phi (theta : angle with z axis,         **/
/**   phi : angle with x axis; i.e, from spherical coordinate form of the   **/
/**   ECEF Reference frame.                                                 **/
/**-------------------------------------------------------------------------**/

/**SPACECRAFT KINEMATIC AND DYNAMIC VARIABLES-------------------------------**/
type *r_sat_eci;        /*Vector to satellite in ECI*/
type *r_sat_ecef;       /*Vector to satellite in ECEF*/
type *v_sat_eci;        /*Velocity vector of satellite w.r.t ECI*/
/*QUATERNIONS*/
type *q_brf_eci;    /*Quaternion of satellite BRF w.r.t ECI*/
type *q_ned_eci;        /*Quaternion of ned w.r.t ECI*/
type *q_drf_eci;        /*Quaternion of desired CS w.r.t ECI*/
type *q_brf_drf;        /*Error quaternion; quaternion of BRF w.r.t DRF; (q_drf_eci)^(-1).(q_brf_eci)  TO/FROM*/
/* DETERMINED STATE VARIABLES */
type *qd_brf_eci;
type *qd0_brf_eci;
type *wd_brf_eci;
type *qd_brf_drf;
type *wd_brf_drf;
/* FILTERED STATE VARIABLES */
type *qf_brf_eci;
type *wf_brf_eci;
/*ALL ANGULAR RATES IN BRF*/
type *w_ned_eci;        /*Angular rates of satellite Orbit Reference Frame, NED [North-East-Down];*/
                        /*These are the angular rates induced due to the satellites motion around the earth in the orbit*/
type *w_brf_eci;	/*Reference angular rates (of body reference frame) in BRF*/
type *w_drf_eci;        /*Reference angular rates (of desired reference frame) in BRF*/
type *w_brf_drf;        /*Angular rate of the satellite w.r.t to the desired coordinate*/
type *w_orbit;          /*Orbital Angular velocity, measured in BRF*/

type *alpha_drf_eci;    /*Angular acceleration of desired reference frame*/

/*-------------------------------------------------------------------------*/

/*SPACECRAFT STRUCTURE VARIABLES-------------------------------------------*/
type Msat = 10;         /*Mass of Satellite in Kg*/
type H = 0.3;           /*Height of satellite*/
type **I, **Iinv;               /*Inertia Tensor (Matrix)*/
type **eye,**eye7;
void setupinertia()
{
    int i=0;
    I = (type**)malloc( 3*sizeof(type*) );
    eye = (type**)malloc(3*sizeof(type*));
    eye7 = (type**)malloc(7*sizeof(type*));
    for( i=0;i<7;i++ )
    {
      eye7[i] = (type*)calloc(7,sizeof(type));
      eye7[i][i] = 1;
    }
    for( i=0; i<3; i++ )
    {
        eye[i] = (type*)malloc(3*sizeof(type));
        I[i] = (type*)malloc( 3*sizeof(type) );
    }
    I[0][0] = 0.018;    I[0][1] = 0;    I[0][2] = 0;   /*x-max I*/
    I[1][0] = 0;    I[1][1] = 0.0018;    I[1][2] = 0;
    I[2][0] = 0;    I[2][1] = 0;    I[2][2] = 0.003;    /*z-min I*/

    eye[0][0] = 1;	eye[0][1] = 0;	eye[0][2] = 0;
    eye[1][0] = 0;	eye[1][1] = 1;	eye[1][2] = 0;
    eye[2][0] = 0;	eye[2][1] = 0;	eye[2][2] = 1;   
}
type *vec_gbrf_brf;     /*Vector from origin of BRF to the origin of gBRF;*/
                        /*Vector from center of mass to geometric center;*/
                        /*Vector measured in BRF*/
/*-------------------------------------------------------------------------*/

/*SUN MODEL VARIABLES------------------------------------------------------*/
type *sun_ECEF;         /* Sun vector from Earth in ECEF frame*/
type *sun_ECI;          /* Sun vector from Satellite in ECI frame*/
type *sun_brf;         /* Sun vector from satellite in BRF frame*/
/*-------------------------------------------------------------------------*/

/*TWO LINE ELEMENT VARIABLES-----------------------------------------*/
struct sat
{
    type epoch_year;
    type epoch_day;
    type frac_year;     /*Current fractional year*/
    type inc;           /* in degrees*/
    type RAAN;          /* in degrees*/
    type e;
    type a;             /*in km*/
    type arg_perigee;   /* in degrees*/
    type mn_anomaly;    /* in degrees*/
    type rev_per_day;
    type rev_epoch;
}sat1;
/*-------------------------------------------------------------------------*/

/*SUN SENSOR VARIABLES-----------------------------------------------------*/
type *q_gbrf_srfi;      /*Quaternion of Geometric Body Reference frame w.r.t Sensor Reference Frame(SRFi)*/
type *vec_srfi_gbrf;    /*vector of SRFi from GBRF (Geometric Center)*/
type *vec_sun_srfi;     /*Sun Vector in SRFi frame*/
type *vec_sun_i_brf;   /*Sun Vector in brf frame coming from i^{th} sun sensor;*/
type sunsensor_i_on;    /* Sensor will give output only if sunsensor_i_on = 1*/
/*-------------------------------------------------------------------------*/

/*MAGNETIC FIELD VARIABLES-------------------------------------------------*/
type *B_eci;             /*Magnetic field in Eci*/
type *B_brf;             /*Magnetic field in brf*/
type *B_prev;
type *B_dot;
/*-------------------------------------------------------------------------*/

/*MAGNETOTORQUER VARIABLES-------------------------------------------------*/
typedef struct{
                type A;
                type C_r;
                type M_c;
                type R,L;
                type rho_d;
                type sigma_T;
                type *position;
                type *Axis;}torquer;

type actuationperiod;   /*Period of time allotted for PWM pulse actuation*/

torquer T[3];           /*3 Torquers*/
type V;                 /*Amplitude of PWM*/
type *T_t_brf;          /*Torquer torque in BRF*/
/*---------------------------------------------------------------------------*/

/*SLIDING CONTROLLER VARIABLES-----------------------------------------------*/
int D;                  /*Number of disturbance torques expected*/
/*BASED ON DESIGN*/
type **N;               /*The peak values of the disturbances; This is a matrix*/
                        /*with 3 columns (for the x,y and z component disturbances)*/
                        /*and one row for each disturbance torque, amounting to D rows.*/
void take_Nvalues()
{
    int i;
    N = (type**)malloc(D*sizeof(type*));
    for( i=0;i<D;i++ )
    {    N[i] = (type*)malloc(3*sizeof(type));
      N[i][0] = 0;	N[i][1] = 0;	N[i][2] = 0;	}

    /*N[0][0] = 0;    N[0][1] = 0;    N[0][2] = 0;
    N[1][0] = 0;    N[1][1] = 0;    N[1][2] = 0;
    N[2][0] = 0;    N[2][1] = 0;    N[2][2] = 0;*/
}
type **ncap;            /*The nominal value of each disturbance torque; This is a matrix*/
                        /*with 3 columns (for the x,y and z component disturbances)*/
                        /*and one row for each disturbance torque, amounting to D rows.*/
void take_ncapvalues()
{
    int i;
    ncap = (type**)malloc(D*sizeof(type*));
    for( i=0;i<D;i++ )
    {    ncap[i] = (type*)malloc(3*sizeof(type));
    	ncap[i][0] = 0;	ncap[i][1] = 0;	ncap[i][2] = 0;	}

    /*ncap[0][0] = 0;    ncap[0][1] = 0;    ncap[0][2] = 0;
    ncap[1][0] = 0;    ncap[1][1] = 0;    ncap[1][2] = 0;
    ncap[2][0] = 0;    ncap[2][1] = 0;    ncap[2][2] = 0;*/
}
type **eta;             /*Error margin matrix; to be fixed based on the required pointing*/
                        /*Accuracy.*/
void take_etavalues()
{
    int i;
    eta = (type**)malloc( 3*sizeof(type*) );
    for( i=0;i<3;i++ )
        eta[i] = (type*)malloc(3*sizeof(type));

    eta[0][0] = 0.00001;    eta[0][1] = 0;    eta[0][2] = 0;
    eta[1][0] = 0;    eta[1][1] = 0.00001;    eta[1][2] = 0;
    eta[2][0] = 0;    eta[2][1] = 0;    eta[2][2] = 0.00001;
}
type **LAMBDA_q;        /*To be fixed according to the required x,y,z order of*/
                        /*convergence in the S plane.*/
void take_LAMBDAqvalues()
{
    int i;
    LAMBDA_q = (type**)malloc( 3*sizeof(type*) );
    for( i=0;i<3;i++ )
        LAMBDA_q[i] = (type*)malloc(3*sizeof(type));

    LAMBDA_q[0][0] = 0.0000002;    LAMBDA_q[0][1] = 0;    LAMBDA_q[0][2] = 0;
    LAMBDA_q[1][0] = 0;    LAMBDA_q[1][1] = 0.0000002;    LAMBDA_q[1][2] = 0;
    LAMBDA_q[2][0] = 0;    LAMBDA_q[2][1] = 0;    LAMBDA_q[2][2] = 0.0000002;
}
/*---------------*/
type **LAMBDA_s;        /*Positive definite matrix of S in control input; fixed once*/
                        /*Design parameters are fixed*/
type *S;                /*Sliding vector*/
void setupcontrollervars()
{
    int i;
    LAMBDA_s = (type**)malloc( 3*sizeof(type*) );
    for( i=0;i<3;i++ )
        LAMBDA_s[i] = (type*)malloc(3*sizeof(type));

    S = (type*)malloc(3*sizeof(type));
}

/*---------------------------------------------------------------------------*/

/*TEMPORARY ARRAYS, MATRICES AND MEMORY MANAGEMENT---------------------------*/

typedef enum{K,F}FREE;    /*K ->Keep the input variables in the function*/
                          /*F ->Free the input variables in the function*/
/*---------------------------------------------------------------------------*/

#endif /* VARIABLES_H_INCLUDED */
