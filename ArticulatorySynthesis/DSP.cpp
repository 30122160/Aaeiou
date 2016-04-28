
//
#include "stdafx.h"
//
//


#include <stdlib.h>
#include <iostream>
#include <cmath>

#include "DSP.h"

using namespace std;

#define pi (4*atan(1.0))

#define LOWPASS  0
#define HIGHPASS 1
#define BANDPASS 2
#define BANDSTOP 3

DSP::DSP( )
{
}

DSP::DSP( double Fs )
{
	setfs( Fs );
}

DSP::~DSP()
{
}


void DSP::r_mp_c( double r, Complex c, Complex *p )
{
	p->re = c.re * r;
	p->im = c.im * r;

	return;
}

void DSP::c_mp_c( Complex c1, Complex c2, Complex *p )
{
	p->re = c1.re * c2.re - c1.im * c2.im;
	p->im = c1.re * c2.im + c1.im * c2.re;

	return;
}

void DSP::r_dv_c( double r, Complex c, Complex *q )
{

	q->im = -r/( pow(c.re, 2) + pow(c.im,2) ) * c.im;
	q->re = r /( pow(c.re, 2) + pow(c.im,2) ) * c.re;

	return;
}

void DSP::c_dv_r( double r, Complex c, Complex *q )
{
	q->re = c.re / r;
	q->im = c.im / r;

	return;
}

void DSP::c_dv_c( Complex c1, Complex c2, Complex *q )
{
	q->re = (c1.re * c2.re + c1.im * c2.im) / ( pow(c2.re,2) + pow(c2.im,2) );
	q->im = (c2.re * c1.im - c1.re * c2.im) / ( pow(c2.re,2) + pow(c2.im,2) );

	return;
}

void DSP::r_plus_c( double r, Complex c, Complex *s )
{
	s->re = c.re + r;
	s->im = c.im;

	return;
}

void DSP::c_plus_c( Complex c1, Complex c2, Complex *s )
{
	s->re = c1.re + c2.re;
	s->im = c1.im + c2.im;

	return;
}

void DSP::poly( Complex *r, double *a, int N )
// function: calculator the coefficient of a polynomial when it's root is given 计算多项式系数
// r     roots of the polynomial
// a     coefficients of the polynomial arranged in the following manner
//       p = a[0]*x^N + a[1]*x^(N-1) + ... + a[N] 
// N     the order of the polynomial
{
	int i, j;

	Complex ctmp1, ctmp2, ctmp3;

	Complex *pc = new Complex[N+1];

	pc[0].re = 1;
	pc[0].im = 0;
	for( i = 1; i <= N; i++ ){
		pc[i].re = 0.0;
		pc[i].im = 0.0;
	}

	for( i = 0; i< N; i++ ){
		for( j =i+1; j>=1; j-- ){
			// c[j] = c[j] - r[i] * c[j-1];
			c_mp_c( r[i], pc[j-1], &ctmp1 );
			r_mp_c( -1.0, ctmp1, &ctmp2 );
			c_plus_c( pc[j],ctmp2, &ctmp3 );

			pc[j].re = ctmp3.re;
			pc[j].im = ctmp3.im;
		}
	}

	for( i = 0; i <= N; i++ ){
		a[i] = pc[i].re;
	}

	delete [] pc;

	return;
}

void DSP::butterworth( int nOrder, double *pWc, int nWc, double fs, int nType, double *a, double *b )
// function: calculate the coefficients of butterworth filter with N-order    计算n阶巴特沃斯滤波器的系数
// nOrder  the order of butterworth filter
// pWc     cuttoff frequency in Hz
// nWc     # number of cutter of frequncey
//         = 1, for low-pass and high-pass filter
//         = 2, for band-pass and band stop filter
// fs      sampling rate
// a, b    coefficients for butterworth filter with required specification
//         H(z) = ( b[0] + b[1]*z^(-1) + ... + b[m]*z^(-m) )/( a[0] + a[1]*z^(-1) + ... + a[n]*z^(-n) );
// fangqiang 2011/02/18 
{
	int i;      
	int nZero;
    int nPole;
	double *Wc_wp;     // warpped cutoff frequency
	double phi;
	
	double K_c;
	Complex K_d;
	Complex ctmp1, ctmp2;

	Complex *p_c  = new Complex[nOrder];
	Complex *z_c  = new Complex[nOrder];
	Complex *p_dr = new Complex[nOrder];  // reciprocal of poles in z-plane
	Complex *z_dr = new Complex[nOrder];  // reciprocal of zeros in z-plane

     nZero = 0;
	 nPole = nOrder;
	// prewarping of the cutoff frequency
	// calculate the corresponding cuttoff frequency of corresponding analog filter
	Wc_wp = new double[nWc];
	for( i = 0; i < nWc; i++ ){
		Wc_wp[i] = 2 * fs * tan( 1.0/2 * 2*pi * pWc[i]/fs );
	}

	// calculate the poles of normalized butterworth filter
	for( i = 0; i < nOrder; i++ ){
		phi = pi/2 + (2*i+1) * pi / (2*nOrder);
		p_c[i].re = cos( phi );
		p_c[i].im = sin( phi );
	}

	// bilinear transfrom of the poles in S-plane
	// calculate the gain, and poles in Z-plane
	if( nType == LOWPASS ){
		K_c = 1.0;
	}
	else if( nType == HIGHPASS ){
		// calculate the polesd and zeros of normalized high-pass filter
		// based on the poles of the normalized low-pass filter
		ctmp2.re = 1.0;
		ctmp2.im = 0.0;
		for( i = 0; i < nOrder; i++ ){
			// calculate the pole of high-pass filter
			r_dv_c( 1.0, p_c[i], &ctmp1 );
			p_c[i].re = ctmp1.re;
			p_c[i].im = ctmp1.im;

			// calculate the zeros of high-pass filter
			z_c[i].re = 0.0;
			z_c[i].im = 0.0;

			// calculate the gain of the high-pass filter
			c_mp_c( p_c[i], ctmp2, &ctmp2 );
		}
		K_c   = ctmp2.re * pow(-1.0,nPole);
		nZero = nPole; 
	}
	else if( nType == BANDPASS ){

	}
	else if( nType == BANDSTOP ){
	}
	
	bilinear( z_c,  p_c, K_c, z_dr, p_dr, &K_d, nZero, nOrder, fs/Wc_wp[0] );

	// calculate a, b coefficient
	poly( p_dr, a, nOrder );
	poly( z_dr, b, nOrder );

	for( i = 0; i <= nOrder; i++ ){
		b[i] = K_d.re * b[i];
		// a[i] = a[i] / K_d.re;
	}

	delete [] p_c;
	delete [] p_dr;
	delete [] z_dr;

	return;
}

void DSP::bilinear( Complex *z_c, Complex *p_c, double k_c, Complex *z_dr, Complex *p_dr, Complex *pk_d, int nZero, int nPole, double k )
{
	int i = 0;
	Complex ctmp1, ctmp2;

	pk_d->re = k_c;
	pk_d->im = 0.0;
	// bilinear transform of poles in S-plane
	for( i = 0; i < nPole; i++ ){
		ctmp1.re = 2 * k - p_c[i].re;
		ctmp1.im = -1.0 * p_c[i].im;

		ctmp2.re = 2 * k + p_c[i].re;
		ctmp2.im = p_c[i].im;

		c_dv_c( ctmp2, ctmp1, &p_dr[i] );
		c_dv_c( *pk_d, ctmp1, pk_d );
	}

	// bilinear transform of zeros in S-plane
	for( i = 0; i < nZero; i++ ){
		ctmp1.re = 2 * k - z_c[i].re;
		ctmp1.im = -1.0 * z_c[i].im;

		ctmp2.re = 2 * k + z_c[i].re;
		ctmp2.im = z_c[i].im;

		c_dv_c( ctmp2, ctmp1, &z_dr[i] );
		c_mp_c( *pk_d, ctmp1, pk_d );
	}

	// and extra zeros 
	for( i = 0; i < nPole-nZero; i++ ){
		z_dr[i].re = -1.0;
		z_dr[i].im = 0.0;
	}

	return;
}


void DSP::filter(double *a, int nOrder_d, double *b, int nOrder_n, double *input, double *output)
// a               point to the array of coefficents of the denorminator分母
// b               point to the array of coefficents of the numerator分子
// nOrder_d        the order of the denorminate
// nOrder_n        the order of the numerator
{
	int i;
	double r;
	double tmp = 0.0;       // numerator of the filter


	for( i = 0; i < nOrder_n; i++ ){
		tmp = tmp + b[i]*input[i];
	}

	r = tmp;
	for( i = 1; i < nOrder_d; i++ ){
		r = r - a[i]*output[i];
	}

	// upadate the output vector
	for( i = nOrder_d-1; i >= 1; i-- ){
		output[i] = output[i-1];
	}
	output[0] = r; 

	return;
}