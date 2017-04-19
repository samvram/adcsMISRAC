#ifndef QUEST_DEFD
#define QUEST_DEFD

type* quest_funct(type* v1, type* w1, type* v2, type* w2)
/*(brf1,eci1,brf2,eci2)*/
				/*OUTPUTS QUATERNION OF BRF w.r.t ECI*/
{
    int i,j;
    type sigma,del,kap;
    type alpha,beta,gamma,w;        /// w = sum of weights

    type **eye = calloc(3,sizeof(type));
    for(i=0; i<3; i++)
        eye[i] = calloc(3,sizeof(type));
    for(i=0; i<3; i++)
        for(j=0; j<3; j++)
            if(i==j)
                eye[i][j]=1;
            else
                eye[i][j]=0;

    type** B = matrix_add(scalar_mult_mat(0.5, vector_mult(w1,trans(v1)),F), scalar_mult_mat(0.5, vector_mult(w2, trans(v2)),F),F,F);  /// CHECKED

    type** S = matrix_add(B, transpose(B),K,F);  /// CHECKED
    sigma = trace(B,F); /// CHECKED
/// sigma = trace of B
    del = determinant(S);  /// CHECKED
/// delta = determinant of S
    kap = trace(adjoint(S),F); /// CHECKED
/// kappa = trace of adjoint of S

    w = 1; /// sum of weights is 1
    alpha = (w*w) - (sigma*sigma) + kap ;
    beta = w - sigma;                   /// from Etica's report (matlab file ; quest.m)
    gamma = (w+sigma)*alpha - del;

    type** temp1;
    type** X = matrix_add(scalar_mult_mat(alpha,eye,F), scalar_mult_mat(beta, S,K),F,F); /// alpha_I + beta_S

    temp1 = matrix_multiplication(S,S,3,3,3,3,K,F); /// S_square
    X = matrix_add(X,temp1,F,F); /// alpha_I + beta_S + S_square

    type* Z = arrsum(arrscalmult(0.5,cross(w1,v1,K,K),3,F), arrscalmult(0.5,cross(w2,v2,K,K),3,F),3,F,F);  /// required Z

    type *X1 = matmult(X,Z,K,F); /// required X

    type* M=calloc(4,sizeof(type));      /// M is penultimate optimal quaternion

    if(gamma > 0)
    {
        M[0] =gamma ;
        M[1] = X1[0] ;
        M[2] = X1[1] ;
        M[3] = X1[2];
    }
    else
    {
        M[0] =-gamma ;
        M[1] = -X1[0] ;
        M[2] = -X1[1] ;
        M[3] = -X1[2];
    }

    type temp2 ;
    temp2 = X1[0]*X1[0] + X1[1]*X1[1] + X1[2]*X1[2] ; /// X1_transpose * X1

    type* Q_opt=calloc(4,sizeof(type));
    type temp;
    temp = sqrt((gamma*gamma) + (temp2));

    for(i=0; i<4; i++)
        Q_opt[i] = M[i]/temp;   /// Required Quaternion

    return Q_opt;
}

#endif
