#ifndef REQUISITES_H_INCLUDED
#define REQUISITES_H_INCLUDED

#include<stdio.h>
#include<stdarg.h>
#include<time.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<ctype.h>

#ifndef M_PI
#define M_PI 3.1416
#endif 

#include"Variables.h"

#define r2d 180.0/M_PI
#define d2r M_PI/180.0

/*Matrix operations-------------------------------------------*/
type determinant(type** a) /* determinant of a 3x3 matrix */
{
    return a[0][0]*((a[1][1]*a[2][2])-(a[1][2]*a[2][1])) -
           a[0][1]*((a[1][0]*a[2][2])-(a[1][2]*a[2][0])) +
           a[0][2]*((a[1][0]*a[2][1])-(a[1][1]*a[2][0]));
}

type** transpose(type** a)  /* transpose of 3x3 matrix*/
{
    int i,j;
    type **c;
    c = (type**)malloc(3*sizeof(type*));
    for(i=0;i<3;i++)
    {
        c[i] = (type*)malloc(3*sizeof(type));
        for(j=0;j<3;j++)
            c[i][j] = a[j][i];
    }

    return c;
}

type** scalar_mult_mat(type a, type** b, FREE FR)   /* scalar*matrix */
{
    int i,j;
    type **c; 
    c = (type**)malloc(3*sizeof(type*));
    for(i=0;i<3;i++)
    {
        c[i] = (type*)malloc(3*sizeof(type));
        for(j=0;j<3;j++)
            c[i][j] = a*b[i][j];
    }
    if( FR==F )
    {
        for( i=0;i<3;i++ )
            free(b[i]);
        free(b);
    }
    return c;
}

type** matrix_multiplication(type** a, type** b, type n, type m, type p, type q, FREE F1, FREE F2)  /* matrix multiplication */
{
    int i,j,k;
    type **c;
    c = (type**)malloc(n*sizeof(type*));
    for(i=0;i<n;i++)
    {
        c[i] = (type*)malloc(q*sizeof(type));
        for(j=0;j<q;j++)
            {
                c[i][j]=0;
                for(k=0;k<p;k++)
                    c[i][j] += a[i][k]*b[k][j];
            }
    }
    if( F1==F )
    {
        for( i=0;i<n;i++ )
            free(a[i]);
        free(a);
    }
    if( F2==F )
    {
        for( i=0;i<p;i++ )
            free(b[i]);
        free(b);
    }
    return c;
}

type** matrix_add(type** a, type** b, FREE F1, FREE F2)  /* addition of 3x3 matrices */
{
    int i,j;
    type **c;
    c = (type**)malloc(3*sizeof(type*));
    for(i=0;i<3;i++)
    {
        c[i] = (type*)malloc(3*sizeof(type));
        for(j=0;j<3;j++)
            c[i][j] = a[i][j] + b[i][j];
    }

    if( F1==F )
    {
        for( i=0;i<3;i++ )
            free(a[i]);
        free(a);
    }
    if( F2==F )
    {
        for( i=0;i<3;i++ )
            free(b[i]);
        free(b);
    }
    return c;
}

type** mxinv( type **M )
{
    int i,j;
    type determinant = 0, **cofactor;
    cofactor = (type**)malloc( 3*sizeof(type*) );
    for( i=0; i<3; i++ )
    {
        cofactor[i] = (type*)malloc( 3*sizeof(type) );
        for( j=0; j<3; j++ )
            cofactor[i][j] = M[(i+1)%3][(j+1)%3]*M[(i+2)%3][(j+2)%3] - M[(i+2)%3][(j+1)%3]*M[(i+1)%3][(j+2)%3];
        determinant += M[0][i]*cofactor[0][i];
    }

    if( determinant==0 )
        return NULL;

    for( i=0; i<3; i++ )
        for( j=0; j<3; j++ )
            cofactor[i][j] = cofactor[j][i]/determinant;

    return cofactor;
}

type trace( type **a,FREE FR )
{
    type t=0;
    int i;
    for( i=0;i<3;i++ )
        t += a[i][i];

    if( FR==F )
    {
      for( i=0;i<3;i++ )
        free(a[i]);
      free(a);
    }
    return t;
}

type** scal_24(type k,type** A, FREE FR)
{
    int i,j;
    type **B;
    B = (type**)malloc( 2*sizeof(type*));
    for(i=0;i<2;i++)
    {
        B[i] = (type*)malloc( 4*sizeof(type));
        for(j=0;j<4;j++)
            B[i][j] = k*A[i][j];
    }

    if( FR==F )
    {
        for( i=0;i<2;i++ )
            free(A[i]);
        free(A);
    }
    return B;
}

type** sum_24(type** A,type** B,FREE F1,FREE F2)
{
    int i,j;
    type **C;
    C = (type**)malloc( 2*sizeof(type*) );
    for(i=0;i<2;i++)
    {
        C[i] = (type*)malloc( 4*sizeof(type) );
        for(j=0;j<4;j++)
            C[i][j] = A[i][j]+B[i][j];
    }

    if( F1==F )
    {
        for( i=0;i<2;i++ )
            free(A[i]);
        free(A);
    }
    if( F2==F )
    {
        for( i=0;i<2;i++ )
            free(B[i]);
        free(B);
    }
    return C;
}

/**** 7X1 and 7X7 matrix operations ***************/

type* scal_17(type k,type* A, FREE FR)
{
    int i;
    type *B;
    B = (type*)malloc( 7*sizeof(type));
    for(i=0;i<7;i++)
      B[i] = k*A[i];

    if( FR==F )
        free(A);
    return B;
}


type** scal_77(type k,type** A, FREE FR)
{
    int i,j;
    type **B;
    B = (type**)malloc( 7*sizeof(type*));
    for(i=0;i<7;i++)
    {
        B[i] = (type*)malloc( 7*sizeof(type));
        for(j=0;j<7;j++)
            B[i][j] = k*A[i][j];
    }

    if( FR==F )
    {
        for( i=0;i<7;i++ )
             free(A[i]);
        free(A);
    }
    return B;
}

type** sum_77(type** A,type** B,FREE F1,FREE F2)
{
    int i,j;
    type **C;
    C = (type**)malloc( 7*sizeof(type*) );
    for(i=0;i<7;i++)
    {
        C[i] = (type*)malloc( 7*sizeof(type) );
        for(j=0;j<7;j++)
            C[i][j] = A[i][j]+B[i][j];
    }

    if( F1==F )
    {
        for( i=0;i<7;i++ )
            free(A[i]);
        free(A);
    }
    if( F2==F )
    {
        for( i=0;i<7;i++ )
            free(B[i]);
        free(B);
    }
    return C;
}

type* sum_17(type* A,type* B,FREE F1,FREE F2)
{
    int i;
    type *C;
    C = (type*)malloc( 7*sizeof(type) );

    for(i=0;i<7;i++)
      C[i] = A[i]+B[i];

    if( F1==F )
        free(A);
    if( F2==F )
        free(B);
 
    return C;
}

type* matmult_7( type **a, type *b, FREE F1, FREE F2 )
{
    int i,j;
    type *c;
    c = (type*)malloc( 7*sizeof(type) );
    for( i=0; i<7; i++ )
    {
      c[i] = 0;
      for( j=0; j<7; j++ )
        c[i] += a[i][j]*b[j];
    }

    if( F1==F )
    {
        for( i=0;i<7;i++ )
            free(a[i]);
        free(a);
    }
    if( F2==F )
        free(b);

    return c;
}


type** transpose_7(type** fac,FREE FR)
{
  int i,j;
  type** b = malloc(7*sizeof(type*));
  for(i=0;i<7;i++)
  {
     b[i] = (type*)malloc(7*sizeof(type));
     for (j=0;j<7;j++)
         b[i][j]=fac[j][i];
  }

    if( FR==F )
    {
        for( i=0;i<7;i++ )
          free(fac[i]);
        free(fac);
    }
    return b;
}

type determinant_7(type** a, type k)
{
  int i,j,m,n,c;
  type s,det;
  s = 1;	det = 0;
  type **b;
  b = malloc(7*sizeof(type*));
  for(i=0;i<7;i++)
    b[i] = (type*)malloc(7*sizeof(type));

  if (k==1)
     det = a[0][0];
  else
    {
     det=0;
     for (c=0;c<k;c++)
       {
        m=0;
        n=0;
        for (i=0;i<k;i++)
            for (j=0;j<k;j++)
              {
                b[i][j]=0;
                if (i != 0 && j != c)
                 {
                   b[m][n]=a[i][j];
                   if (n<(k-2))
                    n++;
                   else
                    {
                     n=0;
                     m++;
                     }
                   }
               }
          det = det + s * (a[0][c] * determinant_7(b,k-1));
          s=-1 * s;
          }
      }

    for( i=0;i<7;i++ )
        free(b[i]);
    free(b);

    return det;
}

type** cofactor_7(type** num)
{
 int i,p,q,m,n,j;
 type **b,**fac;
 b = malloc(7*sizeof(type*));
 fac = malloc(7*sizeof(type*));
  for(i=0;i<7;i++)
  {
        b[i] = (type*)malloc(7*sizeof(type));
        fac[i] = (type*)malloc(7*sizeof(type));
  }

 for (q=0;q<7;q++)
   for (p=0;p<7;p++)
    {
     m=0;
     n=0;
     for (i=0;i<7;i++)
       for (j=0;j<7;j++)
          if (i != q && j != p)
          {
            b[m][n]=num[i][j];
            if (n<(7-2))
             n++;
            else
             {
               n=0;
               m++;
               }
            }
      fac[q][p]=pow(-1,q + p) * determinant_7(b,7-1);
    }
  
  for( i=0;i<7;i++ )
      free(b[i]);
  free(b);

  return fac;
}

type** inverse_7(type** num,FREE FR)
{
int i,j;
  type** b;
  type d;
  type** inverse;
  inverse = malloc(7*sizeof(type*));
  for(i=0;i<7;i++)
        inverse[i] = (type*)malloc(7*sizeof(type));

  d = determinant_7(num,7);
  b = cofactor_7(num);
  b = transpose_7(b,F);
  for (i=0;i<7;i++)
  {
     for (j=0;j<7;j++)
        inverse[i][j]=b[i][j] / d;
     free(b[i]);
  }
  free(b);

  if( FR==F )
  {
    for( i=0;i<7;i++ )
      free(num[i]);
    free(num);
  }

return inverse;
}

/*Array operations--------------------------------------------*/
type norm(type *a, int n, FREE FR)
{
    type NORM = 0;
    int i;
    for( i=0;i<n;i++ )
        NORM += pow( a[i], 2 );/*n^2 = a^2 + b^2 + c^2 */
    NORM = sqrt( NORM );

    if( FR==F )
        free(a);
    return NORM;
}

type* arrsum( type *A, type *B, int n,FREE F1, FREE F2 )
{
    int i;
    type* sum;
    sum = (type*)malloc(n*sizeof(type));
    for( i=0;i<n;i++ )
        sum[i] = A[i] + B[i];

    if( F1==F )
        free(A);
    if( F2==F )
        free(B);
    return sum;
}

type* arrscalmult( type k, type *x, int n, FREE FR )
{
    int i;
    type* xm;
    xm = (type*)malloc( n*sizeof(type) );
    for( i=0;i<n;i++ )
        xm[i] = k*x[i];

    if( FR==F )
        free(x);
    return xm;
}

type* trans(type* a)  /* transpose of 3x1 vector */
{
/*  type* c=malloc(3*sizeof(type));

  c[0] = a[0];
  c[1] = a[1];
  c[2] = a[2];

  return c;*/
  return a;
}

type** vector_mult(type* a, type* b) /* vector * vector_transpose  ///HAS TO BE RENAMED DYAD*/
{
    type** c;
    int i;
    c = (type**)malloc(3*sizeof(type*));   
    for(i=0;i<3;i++)
        c[i]=(type*)malloc(3*sizeof(type));

    c[0][0]=a[0]*b[0]; c[0][1]=a[0]*b[1]; c[0][2]=a[0]*b[2];
    c[1][0]=a[1]*b[0]; c[1][1]=a[1]*b[1]; c[1][2]=a[1]*b[2];
    c[2][0]=a[2]*b[0]; c[2][1]=a[2]*b[1]; c[2][2]=a[2]*b[2];

    return c;
}

/*Multiplication operations-----------------------------------*/
type dot( type *A, type *B )
{
    type Dot=0;
    int i = 0;
    for(i=0;i<3;i++)
        Dot += A[i]*B[i];
    return Dot;
}

type* cross( type *A, type *B, FREE F1, FREE F2 )
{
    int i;
    type* Cp;
    Cp = malloc(3*sizeof(type));
    for( i=0; i<3; i++)
        Cp[i] = A[(i+1)%3]*B[(i+2)%3] - B[(i+1)%3]*A[(i+2)%3];

    if( F1==F )
        free(A);
    if( F2==F )
        free(B);
    return Cp;
}

type* matmult( type **a, type *b, FREE F1, FREE F2 )
{
    int i,j;
    type *c;
    c = (type*)malloc( 3*sizeof(type) );
    for( i=0; i<3; i++ )
    {
        c[i] = 0;
        for( j=0; j<3; j++ )
            c[i] += a[i][j]*b[j];
    }

    if( F1==F )
    {
        for( i=0;i<3;i++ )
            free(a[i]);
        free(a);
    }
    if( F2==F )
        free(b);
    return c;
}

/*Quaternion operations--------------------------------------*/
type* qmult( type *p, type *q,FREE F1,FREE F2 )
{
    type* pq;  /*quaternions p and q are multiplied and  stored in pq.*/
    pq = malloc( 4*sizeof(type) );

    pq[0] = p[0]*q[0] - (p[1]*q[1] + p[2]*q[2] + p[3]*q[3]);
    pq[1] = p[2]*q[3] - p[3]*q[2] + p[0]*q[1] + q[0]*p[1];
    pq[2] = p[3]*q[1] - p[1]*q[3] + p[0]*q[2] + q[0]*p[2];
    pq[3] = p[1]*q[2] - p[2]*q[1] + p[0]*q[3] + q[0]*p[3];

    if( F1==F )
        free(p);
    if( F2==F )
        free(q);
    return pq;
}

type* qconj(type *q)
{
    type* qstar;
    qstar = (type*)malloc(4*sizeof(type));

    qstar[0] =  q[0];
    qstar[1] = -q[1];
    qstar[2] = -q[2];
    qstar[3] = -q[3];

    return qstar;
}

type* qinv( type *p,FREE FR )
{
    type* pinv;
    pinv = (type*)malloc(4*sizeof(type));
    type N = pow( norm( p,4,K ),2 );
    pinv[0] =  p[0]/N;
    pinv[1] = -p[1]/N;
    pinv[2] = -p[2]/N;
    pinv[3] = -p[3]/N;

    if( FR==F )
      free(p);

    return pinv;
}

type* qtransform( type* Q,type* v, FREE F1, FREE F2 )
{
    type *QV,*ret;
    QV = malloc( 4*sizeof(type) );
    ret = malloc(3*sizeof(type));
    QV[0] = 0;
    QV[1] = v[0];
    QV[2] = v[1];
    QV[3] = v[2];

    QV = qmult( qinv(Q,K),QV,F,F ); /**/
    QV = qmult( QV,Q,F,F1 );       /**/

    ret[0] = QV[1];
    ret[1] = QV[2];
    ret[2] = QV[3];

    free(QV);
    if( F2==F )
        free(v);

    return ret;
}

/*Inter notational transformation----------------------------*/
type* r2q( type **a )
{
    type* q;
    q = (type*)malloc(4*sizeof(type));
    q[0] = 0.5*sqrt(1+trace(a,K));
    q[1] = (a[1][2] - a[2][1])/(4*q[0]);
    q[2] = (a[2][0] - a[0][2])/(4*q[0]);
    q[3] = (a[0][1] - a[1][0])/(4*q[0]);  /*q0,q1,q2,q3*/

   return q;
}

type* q2Eangs( type* q )/*Aerospace sequence ZYX - qz.qy.qx */
{
    type m11,m12,m13,m23,m33;
    m11 = 2*pow( q[0],2 ) + 2*pow( q[1],2 )-1;
    m12 = 2*q[1]*q[2] + 2*q[0]*q[3];
    m13 = 2*q[1]*q[3] - 2*q[0]*q[2];
    m23 = 2*q[2]*q[3] + 2*q[0]*q[1];
    m33 = 2*pow( q[0],2 ) + 2*pow( q[3],2 ) - 1;

    type *angles;
    angles = malloc( 3*sizeof(type) );
    angles[0] = atan( m12/m11 )*180/M_PI;    /*Shi */
    angles[1] = asin( -m13 )*180/M_PI;       /*Theta */
    angles[2] = atan( m23/m33 )*180/M_PI;    /*Phi */

    return angles;
}

type* Eangs2q( type* angles )/*Aerospace sequence ZYX - qz.qy.qx */
{
    type* q;
    q = malloc(4*sizeof(type));

    q[0] = cos(angles[0]*d2r/2)*cos(angles[1]*d2r/2)*cos(angles[2]*d2r/2) + sin(angles[0]*d2r/2)*sin(angles[1]*d2r/2)*sin(angles[2]*d2r/2);
    q[1] = cos(angles[0]*d2r/2)*cos(angles[1]*d2r/2)*sin(angles[2]*d2r/2) - sin(angles[0]*d2r/2)*sin(angles[1]*d2r/2)*cos(angles[2]*d2r/2);
    q[2] = cos(angles[0]*d2r/2)*sin(angles[1]*d2r/2)*cos(angles[2]*d2r/2) + sin(angles[0]*d2r/2)*cos(angles[1]*d2r/2)*sin(angles[2]*d2r/2);
    q[3] = sin(angles[0]*d2r/2);

    if( q[0]<0 )
        q = arrscalmult( -1,q,4,F );
    return q;
}

/*Some Transformations----------------------------------------*/
type* xyztortp( type* xyz )
{
    /*z as the axis; r[0] as radial distance; r[1] as theta, angle with z; r[2] as phi, angle with x*/
    type* rtp;
    rtp = (type*)malloc(3*sizeof(type));
    rtp[0] = norm( xyz,3,K );             /*r*/
    rtp[1] = acos( xyz[2]/rtp[0] ); /*acos(z/r)*/
    rtp[2] = atan( xyz[1]/xyz[0] ); /*atan(y/x)*/

    return rtp;
}

type* rtptoxyz( type* polar, FREE FR )
{
    type* xyz;
    xyz = malloc( 3*sizeof(type) );
    xyz[0] = polar[0]*sin(polar[1])*cos(polar[2]);
    xyz[1] = polar[0]*sin(polar[1])*sin(polar[2]);
    xyz[2] = polar[0]*cos(polar[1]);

    if( FR==F )
        free(polar);

    return xyz;
}

type* ECEFtoECI( type* R_ecef, type decyear )    /*Inputs are the vector in ecef and the UT1 decimal years (time of the prime meridian)*/
{
    /*Mean Sidereal Day 23 h 56 m 4.1 s*/
    type MSD = 23 + (56 + 4.1/60)/60;   /*No. of hours per day*/
    type *R_eci;
    type GMST,theta;                    /*Greenwich Mean Sidereal Time, in hours*/
    GMST = 18.697374558 + 24.06570982441908 * (365.25*(decyear - 2000) );
    theta = GMST/MSD*2*M_PI;              /*Right ascension of the Greenwich meridian.*/

    type *q;
    q = malloc( 4*sizeof(type) );

    q[0] = cos(theta/2);
    q[1] = 0;
    q[2] = 0;
    q[3] = sin(theta/2);

    if( q[0]<0 )
        q = arrscalmult( -1,q,4,F );
    R_eci = qtransform( q,R_ecef,F,K );

    return R_eci;
}

type* ECItoECEF( type* R_eci, type decyear )    /*Inputs are the vector in eci and the UT1 decimal Julian years (time of the prime meridian)*/
{
    /*Mean Sidereal Day 23 h 56 m 4.1 s*/
    type MSD = 23 + (56 + 4.1/60)/60;
    type *R_ecef;
    type GMST,theta;  /*Greenwich Mean Sidereal Time, in hours*/
    GMST = 18.697374558 + 24.06570982441908 * (365.25*(decyear - 2000) ); /*Reference: Wikipedia*/
    theta = GMST/MSD*2*M_PI;  /*Right ascension of the Greenwich meridian.*/

    type *q;
    q = malloc( 4*sizeof(type) );
    q[0] = cos(-theta/2);
    q[1] = 0;
    q[2] = 0;
    q[3] = sin(-theta/2);

    if( q[0]<0 )
        q = arrscalmult( -1,q,4,F );
    R_ecef = qtransform( q,R_eci,F,K );

    return R_ecef;
}

void getq_ned_eci( type *xyz )      /*Computes quaternion which takes vector in eci to ned;*/
{                                   /*quaternion OF NED w.r.t ECI*/
                                    /*xyz is the vector to ned origin measured in eci*/
    int i;                          /*Needs to be called at each time instant*/
    type** R;
    R = (type**)malloc( 3*sizeof(type*) );
    for( i=0;i<3;i++ )
        R[i] = (type*)malloc( 3*sizeof(type) );
    type *rtp;
    rtp = xyztortp(xyz);

    R[0][0] = -cos(rtp[1])*cos(rtp[2]);     R[0][1] = -cos(rtp[1])*sin(rtp[2]);     R[0][2] = sin(rtp[1]);
    R[1][0] = -sin(rtp[2]);                 R[1][1] = cos(rtp[2]);                  R[1][2] = 0;
    R[2][0] = -sin(rtp[1])*cos(rtp[2]);     R[2][1] = -sin(rtp[1])*sin(rtp[2]);     R[2][2] = -cos(rtp[1]);

    free(q_ned_eci);
    q_ned_eci = r2q( R );

    if( q_ned_eci[0]<0 )
	  q_ned_eci = arrscalmult( -1.,q_ned_eci,4,F );
    for( i=0;i<3;i++ )
        free(R[i]);
    free(R);
    free(rtp);
}

type** transpose_6(type** fac)
{
  int i,j;

  type** b;
  b = malloc(6*sizeof(type*));
  for(i=0;i<6;i++)
  {
     b[i] = (type*)malloc(6*sizeof(type));
     for (j=0;j<6;j++)
         b[i][j]=fac[j][i];
  }

    return b;
}

type** transpose_3_6(type** fac)
{
  int i,j;
  type** b;
  b = malloc(6*sizeof(type*));

  for(i=0;i<6;i++)
  {
     b[i] = (type*)malloc(3*sizeof(type));
     for (j=0;j<3;j++)
         b[i][j]=fac[j][i];
  }
    return b;
}

type* matmult_4( type **a, type *b, FREE F1, FREE F2 )
{
    int i,j;
    type *c;
    c = (type*)malloc( 4*sizeof(type) );
    for( i=0; i<4; i++ )
    {
      c[i] = 0;
      for( j=0; j<4; j++ )
        c[i] += a[i][j]*b[j];
    }

    if( F1==F )
    {
        for( i=0;i<4;i++ )
            free(a[i]);
        free(a);
    }
    if( F2==F )
        free(b);

    return c;
}

type determinant_6(type** a, type k)
{
  int i,j,m,n,c;
  type s=1,det=0;
  type **b;
  b = malloc(6*sizeof(type*));
  s = 1;	det = 0;
  for(i=0;i<6;i++)
        b[i] = (type*)malloc(6*sizeof(type));

  if (k==1)
     return (a[0][0]);
  else
    {
     det=0;
     for (c=0;c<k;c++)
       {
        m=0;
        n=0;
        for (i=0;i<k;i++)
            for (j=0;j<k;j++)
              {
                b[i][j]=0;
                if (i != 0 && j != c)
                 {
                   b[m][n]=a[i][j];
                   if (n<(k-2))
                    n++;
                   else
                    {
                     n=0;
                     m++;
                    }
                   }
               }
          det=det + s * (a[0][c] * determinant_6(b,k-1));
          s=-1 * s;
          }
    }

    return (det);
}

type** cofactor_6(type** num)
{
  int i,p,q,m,n,j;
  type f;
  type** b;
  type** fac = malloc(6*sizeof(type*));
  b = malloc(6*sizeof(type*)); 
  fac = malloc(6*sizeof(type*));
  f = 6;
  for(i=0;i<6;i++)
  {
        b[i] = (type*)malloc(6*sizeof(type));
        fac[i] = (type*)malloc(6*sizeof(type));
  }

 for (q=0;q<f;q++)
 {
   for (p=0;p<f;p++)
    {
     m=0;
     n=0;
     for (i=0;i<f;i++)
       for (j=0;j<f;j++)
          if (i != q && j != p)
          {
            b[m][n]=num[i][j];
            if (n<(f-2))
             n++;
            else
             {
               n=0;
               m++;
             }
          }
      fac[q][p]=pow(-1,q + p) * determinant_6(b,f-1);
    }
  }

  return fac;
}

type** inverse_6(type** num)
{

  int i,j;
  type **b,d,**inverse;
  inverse = malloc(6*sizeof(type*));
  for(i=0;i<6;i++)
        inverse[i] = (type*)malloc(6*sizeof(type));

  d=determinant_6(num,6);

  b = cofactor_6(num);
  b = transpose_6(b);	/* LEAK */
  for (i=0;i<6;i++)
     for (j=0;j<6;j++)
        inverse[i][j]=b[i][j] / d;
  return inverse;
}

type** matrix_add_6(type** a, type** b, FREE F1, FREE F2)  /* addition of 3x3 matrices */
{
    int i,j;
    type** c;
    c = (type**)malloc(6*sizeof(type*));

    for(i=0;i<6;i++)
    {
        c[i] = (type*)malloc(6*sizeof(type));
        for(j=0;j<6;j++)
            c[i][j] = a[i][j] + b[i][j];
    }

    if( F1==F )
    {
        for( i=0;i<6;i++ )
            free(a[i]);
        free(a);
    }
    if( F2==F )
    {
        for( i=0;i<6;i++ )
            free(b[i]);
        free(b);
    }
    return c;
}

type* matmult_6( type **a, type *b, FREE F1, FREE F2 )
{
    int i,j;
    type *c;
    c = (type*)malloc( 6*sizeof(type) );
    for( i=0; i<6; i++ )
    {
      c[i] = 0;
      for( j=0; j<6; j++ )
        c[i] += a[i][j]*b[j];
    }

    if( F1==F )
    {
        for( i=0;i<6;i++ )
            free(a[i]);
        free(a);
    }
    if( F2==F )
        free(b);
    return c;
}

type** adjoint( type **M )
{
    type determinant, **cofactor;
    int i,j;

    determinant = 0;
    cofactor = (type**)malloc( 3*sizeof(type*) );
    for( i=0; i<3; i++ )
    {
        cofactor[i] = (type*)malloc( 3*sizeof(type) );
        for( j=0; j<3; j++ )
            cofactor[i][j] = M[(i+1)%3][(j+1)%3]*M[(i+2)%3][(j+2)%3] - M[(i+2)%3][(j+1)%3]*M[(i+1)%3][(j+2)%3];
        determinant += M[0][i]*cofactor[0][i];
    }

    return transpose(cofactor);
}

#endif /* REQUISITES_H_INCLUDED */
