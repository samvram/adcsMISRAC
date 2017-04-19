#ifndef FILTER_DEFD
#define FILTER_DEFD

type* kalman_filter()
{
/*    printf("STARTED\n");	*/
    int i;

    type* Z0 = calloc(7,sizeof(type)) ;
    Z0[0] = 0.18 ;
    Z0[1] = 0.18 ;
    Z0[2] = 0.12 ;
    Z0[3] = 0.1 ;
    Z0[4] = -0.2 ;
    Z0[5] = 0.3 ;
    Z0[6] = sqrt(1 - 0.2*0.2 - 0.3*0.3 - 0.1*0.1);

    int j;
    type t_sample = 0.1 ;

    type **Im = (type**)calloc(3,sizeof(type*)) ;   /// Inertia matrix
    for(i=0;i<3;i++)
        Im[i] = (type*)calloc(3,sizeof(type));

    Im[0][0] = 1.2 ; Im[0][1] = 0;   Im[0][2] = 0;
    Im[1][0] = 0 ;   Im[1][1] = 1.4; Im[1][2] = 0;
    Im[2][0] = 0 ;   Im[2][1] = 0;   Im[2][2] = 1.6;

    type* X = calloc(7,sizeof(type)) ;
    X[0] = 0.1 ;
    X[1] = 0.1 ;
    X[2] = 0.2 ;
    X[3] = 0.2 ;
    X[4] = -0.3 ;
    X[5] = 0.1 ;
    X[6] = sqrt(1 - 0.2*0.2 - 0.3*0.3 - 0.1*0.1);

    type **P = (type**)calloc(7,sizeof(type*))  ;     /// covariance matrix
    for(i=0;i<7;i++)
        P[i] = (type*)calloc(7,sizeof(type));

    for(i=0;i<7;i++)
        for(j=0;j<7;j++)
            if(i==j)
                P[i][j] = 1;
            else
                P[i][j] = 0;

    type **F_k = (type**)calloc(7,sizeof(type*));
    for(i=0;i<7;i++)
            F_k[i] = (type*)calloc(7,sizeof(type));

    F_k[0][0] = 0;                  F_k[0][1] = -Im[2][2]/Im[0][0];     F_k[0][2] = Im[1][1]/Im[0][0];      F_k[0][3] = 0;      F_k[0][4] = 0;        F_k[0][5] = 0;      F_k[0][6] = 0;
    F_k[1][0] = Im[0][0]/Im[1][1];    F_k[1][1] = 0;                    F_k[1][2] = Im[2][2]/Im[1][1];      F_k[1][3] = 0;      F_k[1][4] = 0;        F_k[1][5] = 0;      F_k[1][6] = 0;
    F_k[2][0] = Im[1][1]/Im[2][2];    F_k[2][1] = Im[0][0]/Im[2][2];      F_k[2][2] = 0;                    F_k[2][3] = 0;      F_k[2][4] = 0;        F_k[2][5] = 0;      F_k[2][6] = 0;
    F_k[3][0] = 0;                  F_k[3][1] = 0;                    F_k[3][2] = 0;                    F_k[3][3] = 0;      F_k[3][4] = 0.5;      F_k[3][5] = -0.5;   F_k[3][6] = 0.5;
    F_k[4][0] = 0;                  F_k[4][1] = 0;                    F_k[4][2] = 0;                    F_k[4][3] = -0.5;   F_k[4][4] = 0;        F_k[4][5] = 0.5;    F_k[4][6] = 0.5;
    F_k[5][0] = 0;                  F_k[5][1] = 0;                    F_k[5][2] = 0;                    F_k[5][3] = 0.5;    F_k[5][4] = -0.5;     F_k[5][5] = 0;      F_k[5][6] = 0.5;
    F_k[6][0] = 0;                  F_k[6][1] = 0;                    F_k[6][2] = 0;                    F_k[6][3] = -0.5;   F_k[6][4] = -0.5;     F_k[6][5] = -0.5;   F_k[6][6] = 0;


    type **eye = (type**)calloc(7,sizeof(type*));    /// Identity matrix
    for(i=0;i<7;i++)
        eye[i] = (type*)calloc(7,sizeof(type));

    for(i=0;i<7;i++)
        for(j=0;j<7;j++)
        {
            if(i==j)
              eye[i][j]=1;
            else
              eye[i][j]=0;
        }

    type** phi =  (type**)calloc(7,sizeof(type*))  ;     /// phi matrix
    for(i=0;i<7;i++)
        phi[i] = (type*)calloc(7,sizeof(type));

    for(i=0;i<7;i++)
        for(j=0;j<7;j++)
            phi[i][j] = eye[i][j] + F_k[i][j]*t_sample;

    type** Q = (type**)calloc(7,sizeof(type*))  ;     /// Q matrix
    for(i=0;i<7;i++)
        Q[i] = (type*)calloc(7,sizeof(type));

    for(i=0;i<7;i++)
        for(j=0;j<7;j++)
        {
            if(i==j)
              Q[i][j] = 0.05;
            else
              Q[i][j] = 0;
        }

    type** R = (type**)calloc(7,sizeof(type*)) ;      /// R matrix
    for(i=0;i<7;i++)
        R[i] = (type*)calloc(7,sizeof(type));

    for(i=0;i<7;i++)
        for(j=0;j<7;j++)
            if(i==j)
               R[i][j] = 0.04*0.04;
            else
               R[i][j] = 0;

    type** H_k = (type**)calloc(7,sizeof(type*));
    for(i=0;i<7;i++)
        H_k[i] = (type*)calloc(7,sizeof(type));

    for(i=0;i<7;i++)
        for(j=0;j<7;j++)
            H_k[i][j] = eye[i][j];

    type* Z = (type*)calloc(7,sizeof(type));
    type** M;
    type* Y;
    type** S_k;
    type** KG;
    int k;

    for(k=0;k<5;k++)
    {
        M = sum_77( matrix_multiplication( phi, matrix_multiplication(P,transpose_7(phi,K),7,7,7,7,K,F),7,7,7,7,K,F ),Q,F,K);
        X = matmult_7(phi,X,K,F);
        for(j=0;j<7;j++)
            Z[j] = Z0[j] + 0.002;//rand()%0.002;

        Y =  sum_17( Z,scal_17( -1,matmult_7(H_k,X,K,K),F ),K,F );
        S_k = sum_77( matrix_multiplication( H_k,
         matrix_multiplication( M,transpose_7(H_k,K),7,7,7,7,K,F ),7,7,7,7,K,F ),R,F,K );
        KG = matrix_multiplication( M,matrix_multiplication( transpose_7(H_k,K),inverse_7(S_k,F),7,7,7,7,F,F ), 7,7,7,7,K,F );
        X = sum_17( X, matmult_7( KG,Y,K,F ), F, F );
        for( i=0;i<7;i++ )
             free(P[i]);
        free(P);
        P = matrix_multiplication( sum_77( eye, scal_77( -1, matrix_multiplication(KG,H_k,7,7,7,7,F,K),F ),K,F ),M,7,7,7,7,F,F );
    }

    for(i=0;i<7;i++)
    {
     //for(j=0;j<7;j++)
           printf("%lf \t",X[i]);
       printf("\n");
    }

    for(i=0;i<7;i++)
    {
        free(Q[i]);
        free(P[i]);
    }
    free(Q);
    free(P);
    for(i=0;i<3;i++)
        free(Im[i]);
    free(Im);
    for(i=0;i<7;i++)
        free(F_k[i]);
    free(F_k);
    for(i=0;i<7;i++)
        free(eye[i]);
    free(eye);
    for(i=0;i<7;i++)
        free(phi[i]);
    free(phi);
    for(i=0;i<7;i++)
        free(R[i]);
    free(R);
    for(i=0;i<7;i++)
        free(H_k[i]);
    free(H_k);
    free(Z0);
    free(Z);

    return X;
}

#endif
