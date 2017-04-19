double step;

#include"Requisites.h"
#include"Variables.h"
#include"IGRF.h"
#include"Controller.h"
#include"KinematicsDynamics.h"
#include"Propagation.h"
#include"Sensors.h"
#include"AttitudeDetermination.h"
#include"Actuator.h"

double PRDS,Runtime;
double SPINRATE;
int SPIN;

int main(int args,char* argv[])
{
  PRDS = 4; Runtime = PRDS*PERIOD;
  SPINRATE = 0; SPIN = 0;
  step = 0.1;
  /* PARSING INPUT */
  int i;

  for( i=1;i<args;i++ )
  {
    argv[i][1] = tolower(argv[i][1]);
    switch(argv[i][1]){
      case 'p': sscanf(argv[i],"-p%lf",&PRDS);
                Runtime = PRDS*PERIOD; break;
      case 'r': sscanf(argv[i],"-r%lf",&Runtime); break;
      case 's': sscanf(argv[i],"-s%lf,%d",&SPINRATE,&SPIN); break;
      case 't': sscanf(argv[i],"-t%lf",&step); break;
      case 'h': fprintf(stderr,"h; b; p(); r(); s(,); t(); NOSPACE\n"); exit(1);
      case 'b': i=1000; break;
    }
  }

  fprintf(stderr,"p(%lf); r(%lf); s(%lf,%d); t(%lf)\n",PRDS,Runtime,SPINRATE,SPIN,step);
  if( i>=1000)  exit(2);

  printf("This starts here.");

    type t;                 /*time vector*/
    t = 0;
    setupinertia();
    Iinv = mxinv(I);
    setupcontrollervars();
    Init_Torquers();

    D = 1;                      /*Number of Disturbance Torques*/

    alpha_drf_eci = (type*)malloc( 3*sizeof(type) );
    alpha_drf_eci[0] = 0;
    alpha_drf_eci[1] = 0;
    alpha_drf_eci[2] = 0;

    q_brf_eci = (type*)malloc( 4*sizeof(type) );
    q_brf_eci[0] = 1/sqrt(3);
    q_brf_eci[1] = 1/sqrt(3);
    q_brf_eci[2] = 0;
    q_brf_eci[3] = 1/sqrt(3);

    qd_brf_eci = (type*)malloc(4*sizeof(type));
    qd0_brf_eci = (type*)malloc(4*sizeof(type));

    for( i=0;i<4;i++ )
    {
      qd_brf_eci[i] = q_brf_eci[i];
      qd0_brf_eci[i] = qd_brf_eci[i];
    }

    r_sat_eci = malloc( 3*sizeof(type) );
    v_sat_eci = malloc( 3*sizeof(type) );
    w_orbit = malloc( 3*sizeof(type) );

w_brf_eci = (type*)calloc( 3,sizeof(type) );

/*    w_brf_eci[0] = SPINRATE;
    w_brf_eci[1] = SPINRATE;
    w_brf_eci[2] = SPINRATE; */
    w_brf_eci[SPIN] = SPINRATE;

    B_prev = (type*)malloc( 3*sizeof(type) );

    /*Sliding Controller*/
    take_etavalues();
    take_LAMBDAqvalues();
    take_ncapvalues();
    take_Nvalues();
    /*------------------*/

    type* Mreq;


    FILE *fst;
    fst = fopen("STATES.dat","w+");
    fprintf( fst,"Time\tq0\tqx\tqy\tqz\tqd0\tqdx\tqdy\tqdz\tqe0\tqex\tqey\tqez\n");

    FILE *fw;
    fw = fopen("RATECHECK.dat","w+");
    fprintf( fw,"Time\twx\twy\twz\twdx\twdy\twdz\n");

    FILE *fc;
    fc = fopen("CONTROL.dat","w+");
    fprintf(fc,"Time\tqe0\tqex\tqey\tqez\twex\twex\twey\twez\tsang\tSx\tSy\tSz\n");

    FILE *fe;
    fe = fopen("EFFORT.dat","w+");
    fprintf(fe,"Time\tM_x\tM_y\tM_z\n");

    FILE* fm;
    fm = fopen("MAG.dat","w+");
    fprintf(fm,"Time\tBeci_x\tBeci_y\tBeci_z\tBbrf_x\tBbrf_y\tBbrf_z\n");

    FILE* fs;
    fs = fopen("SUN.dat","w+");
    fprintf(fs,"Time\tSuneci_x\tSuneci_y\tSuneci_z\tSunbrf_x\tSunbrf_y\tSunbrf_z\n");

    FILE* fo;
    fo = fopen("ORBIT.dat","w+");
    fprintf(fo,"Time\tX\tY\tZ\tXd\tYd\tZd\n");

    type *Eangles;
    type timeinyears;

    /* INITIALIZATION */
    timeinyears = 2015;
    orbit_propagate(t);
    refgen();
    /* MAGNETIC FIELD VECTOR */
    B_brf = igrf12( timeinyears,r_sat_eci );            /*Get Magnetic field in BR*/
    //fprintf(stderr,"(%lf,%lf,%lf)\n",B_brf[0],B_brf[1],B_brf[2]);exit(1);
    /* SUN VECTOR */
    SunSensor(t);
    type* tmp;
    /* DETERMINATION - SENSOR FUSION */
    for( i=0;i<4;i++ )
      qd0_brf_eci[i] = qd_brf_eci[i];
    free(qd_brf_eci);
    qd_brf_eci = SensorFusion(2);
    wd_brf_eci = RateDetermination();

    qd_brf_drf = qmult( qinv(q_drf_eci,K),qd_brf_eci,F,K );
    wd_brf_drf = arrsum( wd_brf_eci,arrscalmult(-1,w_drf_eci,3,K),3,K,F);
    T_t_brf = malloc(3*sizeof(double));
    int count = 0;
    for( t=0;t<Runtime;t+=step )
    {
        count++;
        timeinyears = 2015 + t/(60*60*23.9345*365.25);
        T_t_brf[0] = 0; T_t_brf[1] = 0; T_t_brf[2] = 0;
	systemresponse(T_t_brf,step/2);
	free(T_t_brf);
	/* CONTROL ALGORITHM */
        Mreq = sliding_controller();
	#ifdef M_ON
        Mreq = Actuation(Mreq,F);
	#endif
        T_t_brf = cross( Mreq,B_brf,K,K );
	/* SYSTEM INTEGRATION */
        systemresponse( T_t_brf,step/2 );

	/* STATE UPDATION */
        timeinyears += step/(60*60*23.9345*365.25);
	orbit_propagate(t);                                 /*Propagate satellite*/
        refgen();                                           /*Generate Reference*/
        /* MAGNETIC FIELD VECTOR */
        B_prev[0] = B_brf[0];
        B_prev[1] = B_brf[1];
        B_prev[2] = B_brf[2];
        B_brf = igrf12( timeinyears,r_sat_eci );            /*Get Magnetic field in BR*/
        B_dot = arrscalmult( 1.0/step,arrsum( B_brf,arrscalmult(-1,B_prev,3,K),3,K,F ),3,F );
	/* SUN VECTOR */
	SunSensor(t+step);

	/* DETERMINATION - SENSOR FUSION */
	for( i=0;i<4;i++ )
	  qd0_brf_eci[i] = qd_brf_eci[i];
        free(qd_brf_eci);
	qd_brf_eci = SensorFusion(2);
        wd_brf_eci = RateDetermination();

	free(qd_brf_drf);
	qd_brf_drf = qmult( qinv(q_drf_eci,K),qd_brf_eci,F,K );
	free(wd_brf_drf);
	wd_brf_drf = arrsum( wd_brf_eci,arrscalmult(-1,w_drf_eci,3,K),3,K,F);
	/* KALMAN FILTER */
	#ifdef F_ON
	Filter();
	#endif

	/* FILE OUTPUT */
        if( count%10==0 )
        {
            q_brf_drf = qmult( qinv(q_drf_eci,K),q_brf_eci,F,K );
	    //	    printf("(%lf,%lf,%lf,%lf)\n",q_brf_drf[0],q_brf_drf[1],q_brf_drf[2],q_brf_drf[3]);
	    //	    exit(1);
            Eangles = q2Eangs( q_brf_drf );
	    /*            fprintf(f,"%lf\t%lf\t%lf\t%lf\t%lf\t",t,q_brf_drf[0],q_brf_drf[1],q_brf_drf[2],q_brf_drf[3]);
            fprintf(f,"%lf\t%lf\t%lf\t%lf\t",q_brf_eci[0],q_brf_eci[1],q_brf_eci[2],q_brf_eci[3]);
            fprintf(f,"%lf\t%lf\t%lf\t",Eangles[0],Eangles[1],Eangles[2]);
            fprintf(f,"%lf\t%lf\t%lf\t",w_brf_eci[0],w_brf_eci[1],w_brf_eci[2]);
            fprintf(f,"%lf\t%lf\t%lf\t",T_t_brf[0],T_t_brf[1],T_t_brf[2]);
            fprintf(f,"%lf\t%lf\t%lf\t",B_eci[0],B_eci[1],B_eci[2]);
	    fprintf(f,"%lf\t%lf\t%lf\t",B_brf[0],B_brf[1],B_brf[2]);
            fprintf(f,"%lf\t%lf\t%lf\t",Mreq[0],Mreq[1],Mreq[2]);
            fprintf(f,"%lf\t%lf\t%lf\t",S[0],S[1],S[2]);
            fprintf(f,"%lf\t%lf\t%lf\n",r_sat_eci[0],r_sat_eci[1],r_sat_eci[2]);*/

	    fprintf(fc,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%e\t%e\t%e\n",t,q_brf_drf[0],q_brf_drf[1],q_brf_drf[2],q_brf_drf[3],w_brf_eci[0],w_brf_eci[1],w_brf_eci[2],sang,S[0],S[1],S[2]);

	    fprintf(fe,"%lf\t%lf\t%lf\t%lf\n",t,Mreq[0],Mreq[1],Mreq[2]);
	    /*	    fprintf(fs,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",t,
									sun_ECI[0],sun_ECI[1],sun_ECI[2],
									sun_ECEF[0],sun_ECEF[1],sun_ECEF[2],
									sun_brf[0],sun_brf[1],sun_brf[2]);
	    fprintf(fq,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t",t,
								q_brf_eci[0],q_brf_eci[1],q_brf_eci[2],q_brf_eci[3],
								qd_brf_eci[0],qd_brf_eci[1],qd_brf_eci[2],qd_brf_eci[3]);
	    tmp = qmult(q_brf_eci,qinv(qd_brf_eci,K),K,F);
            fprintf(fq,"%lf\t%lf\t%lf\t%lf\n",tmp[0],tmp[1],tmp[2],tmp[3]);
	    free(tmp);
	    fprintf(fw,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",t,w_brf_eci[0],w_brf_eci[1],w_brf_eci[2],
	    wd_brf_eci[0],wd_brf_eci[1],wd_brf_eci[2]);*/

	    fprintf(fo,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",t,r_sat_eci[0],r_sat_eci[1],r_sat_eci[2],v_sat_eci[0],v_sat_eci[1],v_sat_eci[2]);
	    fprintf(fm,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",t,B_eci[0],B_eci[1],B_eci[2],B_brf[0],B_brf[1],B_brf[2]);
	    fprintf(fs,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",t,sun_ECI[0],sun_ECI[1],sun_ECI[2],sun_brf[0],sun_brf[1],sun_brf[2]);


            free(Eangles);
            free(q_brf_drf);
        }
        printf("%.4f%%",t/(Runtime)*100);
        printf("\tNorm S\t: %e\n",norm(S,3,K) );
        free(Mreq);
	//        free(T_t_brf);
        free(B_dot);
    }
    //    fclose(f);
    fclose(fs);
    fclose(fst);
    fclose(fw);

    /* PLOTTING RESULTS */
    /*system("gnuplot Plot.gp");
    system("evince *pdf");*/
    /*----*/

    return 0;
}

