#ifndef _MATRIX_H_
#define _MATRIX_H_

#include <vector>
using namespace std;

class CMatrix
{
public:
	CMatrix(void);
	~CMatrix(void);
public:
	void Matrix_Mu_Vector( vector< vector<double> > M, vector<double> x, vector<double> &v, int dim ); 
	void Matrix_Mu_Matrix( vector< vector<double> > M1, vector< vector<double> > M2, vector< vector<double> > &M,int dim_r1, int dim_c1, int &dim_r2, int &dim_c2 );
	void SwapLines( vector< vector<double> > &A );


	void ForwardElimination( double *A, double *b, int dim );
	void BackwardSubstitution( double *A, double *b, double *x, int dim );
	
};

#endif