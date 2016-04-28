#ifndef _DSP_H_
#define _DSP_H_


typedef struct{
	double re;
	double im;
}Complex;

class DSP{
public:
	DSP();
	DSP( double Fs );
	~DSP();
private:
	double fs;

public:
	// set the sampling frequency
	void setfs( double Fs ){ fs = Fs; return; };

	// calculate the coefficients of polynomial
	void poly( Complex *r, double *a, int N );

	// bilinear transform of the
	void bilinear( Complex *z_c, Complex *p_c, double k_c, Complex *z_d, Complex *p_d, Complex *k_d, int nZero, int nPole, double fs );
	
	// lowpass butterworth filter
	void butterworth( int nOrder, double *pWc, int nWc, double fs, int nType, double *a, double *b );

	void filter( double *a, int nOrder_n, double *b, int nOrder_d, double *input, double *output );

	// calculation rule for complex
	void r_mp_c( double r, Complex c, Complex *p );
	void c_mp_c( Complex c1, Complex c2, Complex *p );
	void r_plus_c( double r, Complex c, Complex *s );
	void c_plus_c( Complex c1, Complex c2, Complex *s );
	void r_dv_c( double r, Complex c, Complex *q );
	void c_dv_r( double r, Complex c, Complex *q );
	void c_dv_c( Complex c1, Complex c2, Complex *q );
};

#endif