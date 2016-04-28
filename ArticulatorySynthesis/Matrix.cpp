#include "stdafx.h"




#include "Matrix.h"

#include <iostream>
#include <cmath>
#define  ZERO 1.0e-10

// Mu     short for Multiply
// Mi     short for Minus
// Ad     short for Add
// Inv    short for Inverse

CMatrix::CMatrix(void)
{
}

CMatrix::~CMatrix(void)
{
}

void CMatrix::Matrix_Mu_Vector( vector< vector<double> > M, vector<double> x, vector<double> &v, int dim )
// Matrix_M_Vector is to calculate the result of the production between matrix M and vector v
// dim is the dimenstion of matrix M
{
	int i, j;

	for( i = 0; i < dim; i++ ){
		v[i] = 0.0;
		for( j = 0; j < dim; j++ ){
			v[i] = v[i] + M[i][j] * x[j];
		}
	}

}

void CMatrix::Matrix_Mu_Matrix( vector< vector<double> > M1, vector< vector<double> > M2, vector< vector<double> > &M,int dim_r1, int dim_c1, int &dim_r2, int &dim_c2 )
{
	int i, j, k;

	if( dim_c1 != dim_r2 ){
		cout<<"The number of colume in M1 does not equal to the number of row in M2!";
		exit(0);
	}

	for( i = 0; i < dim_r1; i++ ){
		for( j = 0; j < dim_c2; j++ ){
			M[i][j] = 0;
			for( k = 0; k < dim_c1; k++ ){
				M[i][j] = M[i][j] + M1[i][k] * M2[k][j];
			}
		}
	}
}


void CMatrix::ForwardElimination( double *A, double *b, int dim )
{
	int i, j, k;

	double scale;

	for( i = 0; i < dim; i++ ){
		for( j = i+1; j < dim; j++ ){
			if( fabs( A[j*dim+i] ) > ZERO ){
				// the elementes ahead of the diagonal is set to be zero
				scale = A[i*dim+i] / A[j*dim+i];
				for( k = 0; k <= i; k++ ){
					A[j*dim+k] = 0.0;
				}
				for( k = i+1; k < dim; k++ ){
					A[j*dim+k] = A[j*dim+k] * scale - A[i*dim+k];
				}
				b[j] = b[j] * scale - b[i];
			}
		}
	}
	return;
}


void CMatrix::BackwardSubstitution( double *A, double *b, double *x, int dim )
{
	int i, j;
	
	for( i = dim-1; i>= 0; i-- ){		
		x[i] = b[i];
		for( j = i+1; j < dim; j++ ){
			x[i] = x[i] - A[i*dim+j] * x[j];
		}
		x[i] = x[i] / A[i*dim+i];
	}

	return;
}