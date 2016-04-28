// #pragma once

#ifndef _vt_sim_h_
#define _vt_sim_h_

#include "DSP.h"
#include "Matrix.h"

#define TR 0
#define VT 1
#define NT 2
#define PFL 3
#define PFR 4

#define AREA_MIN 1.0e-10  // the threshold of minimal cross-sectional area, 1.0e-7 m^2 

typedef struct{
	// geometric parameter
	double a[2];  // cross sectional area
	double x[2];  // length
	double s;     // perimeter

	bool bRigid;  // = true, for rigid wall;

	// acoustic elements
	double Ud;    // flow source introduced by the variation of the tube volumn
	double Ls;    // half of the impedance in serier
	double Rs;    // half of the resistance in serier
	double Rs_l;  // resistance on the left arm of the T-type circuit
	double Rs_r;  // resistance on the right arm of the T-type circuit
	double Cp;    // Capacity in paraller
	double Lw;    // impedance introduced by wall vibration
	double Rw;    // resistance introduced by wall vibration
	double Cw;    // capacity introdued  by wall vibration
	double Yw;    // total induction introdued by wall vibration

	// residue introduced by acoustic elements
	double Qls_l; // residue introduced by inductor on the left arm of T-type circuit
	double Qls_r; // residue introduced by inductor on the right arm of T-type circuit
	double Qlw;   // residue introduced by inductor of wall vibration
	double Qcp;   // residue introduced by capacity in parallel
	double Vcw;   // residue introduced by capacity of wall vibration

	// air flow enter and leave the tube, and air pressure at the middle of the tube
	double Ui;    // air stream flow into a tube
	double Uo;    // air stream flow out of a tube
	double u3;    // air stream introduced by wall vibration
	double P;     // air pressure at the middle point of the tube

	// intermidiate variables for computing coefficients of linear equations and residues state 
	double constP;
	double constU;

	double Ns[3]; // dipole noise pressure sources for Maeda model

	// varaibles for dipole and monopole sources in Peter's model
	double Un[3];     // monoploe flow source
	double Pn[3];     // dipole pressure source
	double Amp_u[3];  // amplitude of the noise flow source
	double Amp_p[3];  // amplitude of the noise pressure source
	double rp[3];     // Gaussian noise for pressure source 
	double ru[3];     // Gaussian noise for flow source  

	double Re;        // Reynold # within the tube

}Tube;

typedef struct{
	int nPos;  // position of the constriction along the vocal tract
	double a;  // the cross-sectional area of the constriction
	double x;  // the length of the constriction

	double dist;  // this distance from the outlet of constrion to obstacle
    double eta_n; // coefficient for dipole pressure source
	              // if air stream imping with teeth, eta_n = 1.0; else, eta_n = 0.5;

	int nPos_u; // position of monopole flow source
	int nPos_p; // position of dipole pressure source
}Constriction;


class CVT_sim
{
public:
	CVT_sim( double Fs );
	~CVT_sim(void);
public:

	// variables concerned with the vocal system
	Tube* tr;  // short tubes comprise the trachea
	Tube* gl;  // short tubes for glottis
	Tube* vt;  // short tubes comprise the vocal tract         
	Tube* nt;  // short tubes comprise the nasal tract
	Tube* pfl; // short tubes comprise the left piriform fossa
	Tube* pfr; // short tubes comprise the right piriform fossa

	Constriction Ac;

	int ntr;    // # of short tubes for trachea
	int ngl;    // # of short tubes for glottis
	int nvt;    // # of short tubes comprise the vocal tract from glottis to lips
	int nnt;    // # of short tubes from naso-pharyngeal port to nostril    
	int nvc;    // # of short tubes change their cross-sectional area in the nasal tract
	            // when coupling with vocal tract
	int npfl;   // # of short tubes for left piriform fossa
	int npfr;   // # of short tubes for right piriform fossa

	int nPos_th;  // the position of teeth along the vocal tract
	int nPos_nt;  // the position of naso-pharyngeal along the vocal tract
	              // nPos_nt = 9, if nvt = 17;
	              // # of coupling length between vt and nt is 3cm 
	int nPos_sin; // the position of side cavity along the nasal tract
	              // nPos_sin = 7 in Maeda's case
	int nPos_pf;  // the position of piriform fossa along the vocal tract

	int nS_tr;   // the position of the 1st equation for trachea
	int nS_gl;   // the position of the 1st equation for glottis
	int nS_vt;   // the position of the 1st equation for vocal tract
	int nS_nt;   // the position of the 1st equation for nasal tract
	int nS_pfl;  // the position of the 1st equation for left piriform fossa
	int nS_pfr;  // the position of the 1st equation for right piriform fossa

	double dist_nt;  // the distance from nasal-pharyngeal port to the upper incisor

	// variables related with lips and nostrils
	double Grl, Srl;       // parallel impedance at the lip end
	double Grn, Srn;       // parallel impedance at the nostril end
	double Pl, Pn, Pw;     // pressure radiated at lips and nostrils
	double Ul[2];          // volume velocity at the lip end
	double Un[2];          // volume velocity at the nostril
	double Uw[2];;         // sum of volume velocity introduced by wall vibration
	double Vrl, Vrn;       // residue introduced by radiation impdance at lips and nostril
	double Qls_g;          // residue introduced by inductance of glottis

	double constPl, constUl;
	double constPn, constUn;

	// acoustic elements for nasal sinus
	double R_sin;
	double L_sin;
	double C_sin;
	double Y_sin;     

	double Udc;      // dc airflow
	double Uac;      // flow at the supraglottal constriction

	double Rg;
	double Lg;

	// configurations for connections with subglottal system and side branches
	bool nt_on;
	bool tr_on;
	bool pfr_on;
	bool pfl_on;
	bool wall_rad;    // wall radiation
	bool sin_on;
	bool bSubstitute; // if == true, replace the Impedance of wall vibration with Impedance of nasal Sinus
	                  // if == false, connect the Impedance of Nasal Sinus to Impedace of wall vibration in parallel

	// constants for computing Rs, Ls, and Cp
	double Krs;
	double Kls;
	double Kcp;

	double Klw;
	double Krw;
	double Kcw;

	// constants for computing Rg, Lg
	double Krg;
	double Kdg;
	double Klg;

	// constants for computing Gr, Grn, Sr, Srn
	double Kgr;
	double Ksr;

	// constants for Bernoulli effects
	double Kb;

	// radiated pressure at the distance of 1m
	double Kr;

	// subglottal pressure
	double Ps;
public:
	double fs;   // sampling rate 

	// matrixes for linear equation
	double *A; // coefficient matrix on the right side of linear equation
	double *b; // constants on the right side of linear equation
	double *X; // solved result


	int nEq;   // the # of linear equations

public:
	CMatrix matrix;  // for matrix computation
	DSP dsp;        // for filtering signals

	double Au[3], Bu[3];  // filter coefficient for flow-source
    double Ap[3], Bp[3];  // filter coefficient for pressure-source
public:
	// involve trachea, left and right piriform fossa into the system
	void turnOntract( bool bOn, int nType );

	void Init( double *a, double *x, int nTube );
	void InitResidue( double *ag, int nAg, double Ps );

	void loadconfig( char *fn, int nType ); // load the configuration of the tract
	void acoust_elm_t( Tube *pTube, int nTube, double m, double b, double k ); // calculate the acoustic elements
	void evalNoiseSource( Tube *pTube, int nTube );
	void glottalImpedance( double *ag, int nTube );
	void radiationImpedance( double a, double &gr, double &sr );
	void sinusImpedance( const double vol, const double len, const double apt, double &R_sin, double &L_sin, double &C_sin, double &Y_sin );
	void modifyNT_sin( int nPos_sin, double r_sin, double l_sin, double c_sin,  double y_sin );
	void addBernoulli( Tube *pt, int nTube );
	void evalFlow ( double *x );

	double glottalArea_dy( double T0, double t, double OQ, double SQ, double A0 );
	double glottalArea_dc( );
	double nonZero( double x );
	//double frictionNoise( int Entry, double mem1[4][3], double mem2[4][3], double *mem );
	//void filter( double *x, double *y, double *a, double *b );

	void evalUdc( const double ag, const double Ps );
	void findConstriction( );
	void refresh_ax( double *a, double *x, double *a_nv );
	void refresh_residue( Tube *pt, const int nTube ); 
	void refresh_noise( Tube *pt, const int nTube );
	void refresh( );

	void clearUP( Tube *pt, const int nTube );
	void clearResidue( Tube *pt, const int nTube );
	void clearNoiseSource( Tube *pt, const int nTube );

	double td_sim( double *ag, double *a, double *x, double *a_nv, double Ps );
	void makeMatrix( double Ps );
	double randn( );


	// transform sparse matrix to intel MKL format
	void matrixTran( double *a, double *a_tran );
};

#endif