


#include "stdafx.h"


#include "VT_sim.h"
#include <stdio.h>
#include <iostream>
#include <math.h>
#include <time.h>
#include "PhysicalConstants.h"

using namespace std;

//#define CHECK
//#define CHECK_F
//#define CHECK_B
//#define CHECK_R
//#define CHECK_RE

CVT_sim::CVT_sim( double Fs )
{
	fs = Fs;
	dsp.setfs( fs );

	gl = NULL;
	vt = NULL;
	nt = NULL;
	tr = NULL;
	pfl = NULL;
	pfr = NULL;

	// turn off the branches when initilaize the vocal system初始化声道系统时关闭分支
	nt_on = false;
	tr_on  = false;
	pfl_on = false;
	pfr_on = false;
	sin_on = true;

	wall_rad = true;

	// set the configuration of the system 设置系统配置
	turnOntract( false, NT );
	turnOntract( false, TR );
	turnOntract( false, PFL );
	turnOntract( false, PFR );
}


CVT_sim::~CVT_sim()
{
	if( gl != NULL ){
		delete [] gl;
	}
	if( vt != NULL ){
		delete [] vt;
	}
	if( tr != NULL ){
		delete [] tr;
	}
	if( nt != NULL ){
		delete [] nt;
	}
	if( pfl != NULL ){
		delete [] pfl;
	}
	if( pfr != NULL ){
		delete [] pfr;
	}

	
	delete [] A;
	delete [] X;
	delete [] b;
}

void CVT_sim::matrixTran( double *a, double *a_tran )
{
	// transform sparse matrix into intel MKL format 将稀疏矩阵转换为国际MKL格式
}

void CVT_sim::Init( double *a, double *x, int nTube )
{
	// calculate the coefficient for calculating the time domain acoustic elements of short tubes 计算 计算短管时域声学元素用到的系数
	Krs = 8 * pi * mu_air;// adopted from "Acoustic Phonetics" P27 阻力 粘滞系数

   // for circulkar cross-sectional tube 横断面管
	Kls = ro_air;//kls=空气密度
	Kcp = 1.0 / (ro_air * c_air * c_air);

	// for glottis 声门部分 Ｒ　Ｄ　Ｌ
	Krg = 12 * mu_air * dg * pow(lg,2);
	Kdg = 0.5 * kg * ro_air;
	Klg = ro_air * dg;

	// for radiation at lips/nostrils end  嘴唇／鼻孔处的辐射 P205式（8）
	Kgr = 9 * pow(pi,2) / ( 128 * ro_air * c_air ); //Gr 电导
	Ksr = 3 * pi * sqrt( pi ) / ( 8 * ro_air ); //Sr 电纳

	// for wall vibration 振幅
	Kb = kg * ro_air;
	Kr = ro_air / (2.0*pi); 

	//// set the configuration of the system
	//turnOntract( false, TR );
	//turnOntract( false, PFL );
	//turnOntract( false, PFR );

	// load the configuration of tract system
	// Initialize short tubes of vocal tract
	nvt = nTube;
	vt  = new Tube[nvt];
	for( int i = 0; i < nvt; i++ ){
		vt[i].a[1] = vt[i].a[0] = nonZero( a[i] );
		vt[i].x[1] = vt[i].x[0] = nonZero( x[i] );
		vt[i].s = sqrt( vt[i].a[0]/pi ) * 2 * pi * 2.0;

		vt[i].bRigid = false;
	}
	clearUP( vt, nvt );
	clearResidue( vt, nvt );
	clearNoiseSource( vt, nvt );
	Qls_g = 0.0;
	Vrl = 0.0;
	Ul[0] = Ul[1] = 0.0;

	// Initialize short tubes of nasal tract
	if( nt_on == true ){
		loadconfig( "../config/Nasal_Maeda", NT );
		clearUP( nt, nnt );
		clearResidue( nt, nnt );
		clearNoiseSource( nt, nnt );
		sinusImpedance( vol_sin, len_sin, apt_sin, R_sin, L_sin, C_sin, Y_sin );
		Vrn = 0.0;
		Un[0] = Un[1] = 0.0;
	}
	else{
		nnt = 0;
	}

	// Initialize the sum of volume velocity introduced by wall vibration
	Uw[0] = Uw[1] = 0.0;

	// determin constriction within the vocal tract
	findConstriction( );

	if( tr_on == true ){
		loadconfig( "..\\config\\tr.txt", TR );
		acoust_elm_t( tr, ntr, m_tr, b_tr, k_tr );
		clearUP( tr, ntr );
		clearResidue( tr, ntr );
		clearNoiseSource( tr, ntr );
	}
	if( pfl_on == true ){
		loadconfig( ".\\config\\pfl.txt", PFL );
		acoust_elm_t( pfl, npfl, m_vt, b_vt, k_vt );
		clearUP( pfl, npfl );
		clearResidue( pfl, npfl );
		clearNoiseSource( pfl, npfr );
	}
	if( pfr_on == true ){
		loadconfig( ".\\config\\pfr.txt", PFR );
		acoust_elm_t( pfr, npfr, m_vt, b_vt, k_vt );
		clearUP( pfr, npfr );
		clearResidue( pfr, npfr );
		clearNoiseSource( pfr, npfr );
	} 

	// initialize the random number generator, and set the rand
	srand( (unsigned)time(NULL) );

	// prepare for linear equations
	// Ps  = 8.8e-2 * H2O;
	nEq = nvt+1;                                // nvt+1 equations for vocal tract
	if( nt_on == true )  nEq = nEq + (nnt + 1); // nnt+1 equations for nasal tract
	if( tr_on == true )  nEq = nEq + ntr;
	if( pfl_on == true ) nEq = nEq + npfl;
	if( pfr_on == true ) nEq = nEq + npfr;

	X = new double[nEq];
	b = new double[nEq];
	A = new double[nEq*nEq];

	memset( X, 0, nEq*sizeof(double) );
	memset( b, 0, nEq*sizeof(double) );
	memset( A, 0, nEq*nEq*sizeof(double) );

	// determin the location of the first equation for each component of the vocal system
	// start position of trachea in linear equations
	nS_tr = 0;

	// start position of vocal tract in linear equations
	if( tr_on == true ){
		nS_vt = nS_tr + ntr;
	}
	else{
		nS_vt = 0;
	}

	// start position of the nasal branch in linear equations
	if( nt_on == true ){
		nS_nt = nS_vt + nvt + 1;
	}

	// start position of the left piriform fossa in linear equations
	if( pfl_on == true && nt_on == true ){
		nS_pfl = nS_nt + nnt + 1;
	}
	else if( pfl_on == true && nt_on == false ){
		nS_pfl = nS_vt + nvt + 1;
	}

	// start position of the right piriform fossa in linear equations
	if( nt_on == true && pfl_on == false && pfr_on == true){
		nS_pfr = nS_nt + nnt + 1;
	}
	else if( nt_on == false && pfl_on == false && pfr_on == true ){
		nS_pfr = nS_vt + nvt + 1;
	}
	else if( pfl_on == true && pfr_on == true ){
		nS_pfr = nS_pfl + npfl;
	}

	return;
}


void CVT_sim::InitResidue( double *ag, int nAg, double Ps )
//function: This function is to set the intitial pressure in the tract and corresponding residues
//          accroding the to configuration of the vocal tract
{
	int i;
	double ag_min = 100;
	double a_th = 0.2e-4;  // threshold for the cross-sectional area at constriction
	
	// find the constriction within the glottis
	for( i = 0; i < nAg; i++ ){
		if( ag_min > ag[i] ){
			ag_min = nonZero( ag[i] );
		}
	}

	// find the constriction within the vocal tract
	findConstriction( );

	// Initialize the residue in trachea
	if( ag_min < a_th ){
		if( tr_on == true ){
			for( i = 0; i < ntr; i++ ){
				tr[i].P   = Ps;
				tr[i].Qcp = 2 * fs * tr[i].Cp * tr[i].P;
			}
		}
	}

	// initialzie the residue in vocal tract and piriform fossa
	if( Ac.a < ag_min ){
		for( i = 0; i <= Ac.nPos; i++ ){
			vt[i].P = Ps; 
			vt[i].Qcp = 2 * fs * vt[i].Cp * vt[i].P;
		}

		if( pfl_on == true && Ac.nPos > nPos_pf ){
			for ( i = 0; i < npfl; i++ ){
				pfl[i].P = Ps;
				pfl[i].Qcp = 2 * fs * pfl[i].Cp * pfl[i].P;
			}
		}

		if( pfr_on == true && Ac.nPos > nPos_pf ){
			for( i = 0; i < npfr; i++ ){
				pfr[i].P = Ps;
				pfr[i].Qcp = 2 * fs * pfr[i].Cp * pfr[i].P;
			}
		}
	}

	return;
}

double CVT_sim::td_sim( double *ag, double *a, double *x, double *a_nv, double Ps )
// time domain simulation of the vocal system
{
	int i = 0;
	double ag_min = 100.0;
	double Pr = 0;

	//************************* update the cross-sectional area and length of glottis, vocal tract, and coupling between **********************//
	// vocal tract and nasal tract
	refresh_ax(a, x, a_nv);
	findConstriction( );

	// calculate the DC component of air flow
	for( i = 0; i < ngl; i++ ){
		if( ag_min > ag[i] )
			ag_min = ag[i];
	}

	//********************************************** calculate the DC component of the air flow *************************************************//
	evalUdc( ag_min, Ps );

	//****************************** computing the acoustic elements of the components of the vocal system **************************************//
	// acoustic element for trachea
	if( tr_on == true ){
		acoust_elm_t( tr, ntr, m_tr, b_tr, k_tr );
	}

	// acoustic element for glottis
	glottalImpedance( ag, ngl );

	// acoustic element for vocal tract
	acoust_elm_t( vt, nvt, m_vt, b_vt, k_vt );
	
	// Add Bernoulli resistance to acoustic elements of vocal tract    
	addBernoulli( vt, nvt );
//	vt[0].Rs_l = vt[0].Rs_l + abs(vt[0].Ui) * ro_air / (2*pow(vt[0].a[0],2));

	// acoustic element for nasal tract 
	if( nt_on == true ){
		acoust_elm_t( nt, nnt, m_vt, b_vt, k_vt );
		modifyNT_sin( nPos_sin, R_sin, L_sin, C_sin, Y_sin );
	}

	// acoustic elements for piriform fossa
	if( pfl_on == true ){
		acoust_elm_t( pfl, npfl, m_vt, b_vt, k_vt ); 
	}
	if( pfr_on == true ){
		acoust_elm_t( pfr, npfr, m_vt,b_vt, k_vt );
	}

	// radiation impedance at the lips end
	radiationImpedance( vt[nvt-1].a[0], Grl, Srl );

	// radiation impedance at the nostrils end
	if( nt_on == true ){
		radiationImpedance( nt[nnt-1].a[0], Grn, Srn );
	}

	//******************************************* calculate the monopole and dipole noise sources *********************************************//
	//evalNoiseSource( vt, nvt );

	//********************************************** construct linear equations for vocal system *********************************************//
	makeMatrix( Ps );

	// solve linear equations;
	matrix.ForwardElimination( A, b, nEq );

#ifdef CHECK_F
	FILE *fp_f;
	fp_f = fopen( ".\\forward_c.txt", "w" );
	for( int j = 0; j < nEq; j++ ){
		fprintf( fp_f, "%d:\t", j );
		for( int k = 0; k < nEq; k++ ){
			fprintf( fp_f, "%e\t", A[j*nEq+k] );
		}
		fprintf( fp_f, "%e\n", b[j] );
	}
	fclose( fp_f );
#endif

	matrix.BackwardSubstitution( A, b, X, nEq );
 //   int info;
	//info = LAPACKE_dgetrf( LAPACK_ROW_MAJOR, nEq, nEq, A, nEq, ipiv ); 
	//info = LAPACKE_dgetrs( LAPACK_ROW_MAJOR, 'N', 2, 1, A, nEq, ipiv, b, nEq );

#ifdef CHECK_R
	FILE *fp_r;
	fp_r = fopen("root_c.txt", "w");
	for( int j = 0; j < nEq; j ++ ){
		//fprintf( fp_r, "%e\n", x[j] );
		fprintf( fp_r, "%e\n", X[j] );
	}
	fclose( fp_r );
#endif

	// calculate the air stream flow into and out of each short tube
	// evalFlow( x );  
	evalFlow( X );

	//***************************************************** refresh the vocal system ********************************************************//
	// refresh the residue of glottis, vocal tract, nasal tract, and the noise source in the vocal tract
	refresh( );

	//*************************************************** calculate radiated air pressure ***************************************************//
	// radiation from lips
	Ul[1] = Ul[0]; Ul[0] = vt[nvt-1].Uo;
	Pl = ( Ul[0] - Ul[1] ) * fs;

	// radiation from nostrils
	if( nt_on == true ){
		Un[1] = Un[0]; Un[0] = nt[nnt-1].Uo;
		Pn = ( Un[0] - Un[1] ) * fs;
	}

	//Pl = (vt[nvt-1].Uo + constPl) / (Grl + Srl/(2*fs));
	//Pn = (nt[nnt-1].Uo + constPn) / (Grn + Srn/(2*fs));

	// radiation from wall vibration
	Uw[1] = Uw[0]; 
	Uw[0] = 0;
	for( i = 0; i < nvt; i++ ){
		if( vt[i].bRigid == false ){
			Uw[0] = Uw[0] + vt[i].u3;
		}
	}
	Pw = ( Uw[0] - Uw[1] ) * fs;


	if( nt_on == true && wall_rad == true ){
		Pr = ( Pw + Pn + Pl ) * Kr;
	}
	else if( nt_on == false && wall_rad == true ){
		Pr = ( Pw + Pl ) * Kr;
	}
	else if( nt_on == false && wall_rad == false ){
		Pr = Pl * Kr;
	}

	return Pr;
}

void CVT_sim::evalFlow( double *x )
{
	int i;

	// air flow in trachea
	if( tr_on == true ){
		for( i = 0; i <  ntr; i++ ){
			tr[i].Ui = x[nS_tr+i];
			tr[i].Uo = x[nS_tr+i+1];
		}
	}

	// air flow in left branch of piriform fossa
	if( pfl_on == true ){
		for( i = npfl-1; i >= 0; i-- ){
			pfl[i].Ui = x[nS_pfl+i];
			if( i < npfl-1 ){
				pfl[i].Uo = pfl[i+1].Ui;
			}
			else{
				pfl[i].Uo = 0;
			}
		}
	}

	// air flow in right piriform fossa
	if( pfr_on == true ){
		for( i = npfr-1; i >= 0; i-- ){
			pfr[i].Ui = x[nS_pfr+i];
			if( i < npfr-1 ){
				pfr[i].Uo = pfr[i+1].Ui;
			}
			else{
				pfr[i].Uo = 0;
			}
		}
	}

	// air flow in nasal tract
	if( nt_on == true ){
		for( i = 0; i< nnt; i++ ){
			nt[i].Ui = x[nS_nt+i];
			nt[i].Uo = x[nS_nt+i+1];
		}
	}

	// air flow in vocal tract
	vt[nvt-1].Uo = x[nS_vt+nvt];
	for( i = nvt-1; i >= 0; i-- ){ 
		vt[i].Ui = x[nS_vt+i];

		if( i+1 <= nvt-1 ){
			vt[i].Uo = vt[i+1].Ui;
		}

		if( nt_on == true && i == nPos_nt-1 ){
			// air flow at the nasal-pharyngeal port
			vt[i].Uo = vt[i+1].Uo + nt[0].Ui;
		}
		else if( i == nPos_pf-1 ){
			// air flow at the piriform fossa port
			if( pfl_on == true ){
				vt[i].Uo = vt[i].Uo + pfl[0].Ui;
			}
			if( pfr_on == true ){
				vt[i].Uo = vt[i].Uo + pfr[0].Ui;
			}
		}
	}

	return;
}


void CVT_sim::makeMatrix( double Ps )
{
	int i, j;
	int nIndx;

	double h;

	i = 0;

	//***************************************** linear equations for trachea ****************************************************//
	if( tr_on == true ){
		tr[0].Rs_l = tr[0].Rs_l + Rlung;
		for( j = 0; j < ntr; j++ ){
			if( j == 0 ){
				h = 2 * fs * tr[j].Ls + tr[j].Rs_l; 

				A[i*nEq + (i+1)] = 1/ ( 2 * fs * tr[j].Cp + tr[j].Yw);
				A[i*nEq + i ]    = -1.0 * ( h + A[i*nEq + (i+1)] );

				tr[j].constP = tr[j].Qcp - tr[j].Ud - tr[j].Yw *( tr[j].Qlw - tr[j].Vcw );
				tr[j].constU = tr[j].Qls_l;

				b[i] = A[i*nEq + (i+1)] * tr[j].constP - Ps - tr[j].constU;
			}
			else{
				h = 2 * fs * (tr[j-1].Ls + tr[j].Ls) + (tr[j-1].Rs_r + tr[j].Rs_l); 

				A[i*nEq + (i-1)] = 1/ ( 2 * fs * tr[j-1].Cp + tr[j-1].Yw );
				A[i*nEq + (i+1)] = 1/ ( 2 * fs * tr[j].Cp   + tr[j].Yw );
				A[i*nEq + i]     = -1.0 * ( h + A[i*nEq + (i-1)] + A[i*nEq + (i+1)] );

				tr[j].constP = tr[j].Qcp - tr[j].Ud - tr[j].Yw *( tr[j].Qlw - tr[j].Vcw );
				tr[j].constU = tr[j-1].Qls_r + tr[j].Qls_l;

				b[i] =  A[i*nEq + (i+1)] * tr[j].constP -  A[i*nEq + (i-1)] * tr[j-1].constP - tr[j].constU;
			}

			i = i+1;
		}
	}

	//********************************************* linear equations for vocal tract *****************************************************//
#ifdef CHECK
	FILE *fp_log;
	fp_log = fopen( ".\\matrix_log_c.txt", "w" );
#endif

	for( j = 0; j <= nvt; j++ ){ 
		if( j == 0 ){
			if( tr_on == true ){
				h = 2 * fs * (tr[ntr-1].Ls + + Lg + vt[j].Ls) + (tr[ntr-1].Rs_r + Rg + vt[j].Rs_l);

				A[i*nEq+(i-1)] = 1 / ( 2 * fs * tr[ntr-1].Cp + tr[ntr-1].Yw ); 
				A[i*nEq+(i+1)] = 1 / ( 2 * fs * vt[j].Cp + vt[j].Yw );
				A[i*nEq+i]     = -1.0 * ( h + A[i*nEq+(i-1)] + A[i*nEq+(i+1)] );

				vt[j].constP = vt[j].Qcp - vt[j].Ud - vt[j].Yw * ( vt[j].Qlw - vt[j].Vcw );
				vt[j].constU = tr[ntr-1].Qls_r + Qls_g + vt[j].Qls_l;

				b[i] = A[i*nEq+(i+1)] * vt[j].constP - A[i*nEq+(i-1)] * tr[ntr-1].constP - vt[j].constU;
			}
			else{
				Rg = Rg + Rlung;
				h = 2 * fs * ( Lg + vt[j].Ls ) + ( Rg + vt[j].Rs_l ); 

				A[i*nEq + (i+1)] = 1/ ( 2 * fs * vt[j].Cp + vt[j].Yw);
				A[i*nEq + i ]    = -1.0 * ( h + A[i*nEq + (i+1)] );
				
				vt[j].constP = vt[j].Qcp - vt[j].Ud - vt[j].Yw *( vt[j].Qlw - vt[j].Vcw );
				vt[j].constU = vt[j].Qls_l + Qls_g;

				b[i] = A[i*nEq + (i+1)] * vt[j].constP - Ps - vt[j].constU;
			}
		}
		else if( j == nvt ){
			h = 2 * fs * vt[j-1].Ls + vt[j-1].Rs_r;

			A[i*nEq+(i-1)] = 1.0 / ( 2 * fs * vt[j-1].Cp + vt[j-1].Yw);
			A[i*nEq+i]     = -1.0 * ( h + A[i*nEq+(i-1)] + 1.0/(Grl + Srl/(2*fs)) );

			constPl = -Vrl;
			constUl = vt[j-1].Qls_r;

			b[i] = constPl / (Grl + Srl/(2*fs)) - A[i*nEq+(i-1)] * vt[j-1].constP - constUl;
		}
		else{
			h = 2 * fs * (vt[j-1].Ls + vt[j].Ls) + (vt[j-1].Rs_r + vt[j].Rs_l); 

			A[i*nEq + (i-1)] = 1/ ( 2 * fs * vt[j-1].Cp + vt[j-1].Yw );
			A[i*nEq + (i+1)] = 1/ ( 2 * fs * vt[j].Cp   + vt[j].Yw );
			A[i*nEq + i]     = -1.0 * ( h + A[i*nEq + (i-1)] + A[i*nEq + (i+1)] );

			vt[j].constP = vt[j].Qcp - vt[j].Ud - vt[j].Yw *( vt[j].Qlw - vt[j].Vcw );
			vt[j].constU = vt[j-1].Qls_r + vt[j].Qls_l;

			b[i] =  A[i*nEq + (i+1)] * vt[j].constP -  A[i*nEq + (i-1)] * vt[j-1].constP - vt[j].constU;

			// add coefficients introduced by side branches
			// side branch of nasal tract
			if( nt_on == true ){
				if( j == nPos_nt-1 ){
					A[i*nEq+nS_nt] = A[i*nEq+(i+1)];
				}
				else if( j == nPos_nt ){
					A[i*nEq+nS_nt] = -1.0 * ( 2 * fs * vt[j-1].Ls + vt[j-1].Rs_r + A[i*nEq+(i-1)] );
				}
			}

			// side branch of left piriform fossa
			if( pfl_on == true ){
				if( j == nPos_pf-1 ){
					A[i*nEq+nS_pfl] = A[i*nEq+(i+1)];
				}
				else if( j == nPos_pf ){
					A[i*nEq+nS_pfl] = -1.0 * ( 2 * fs * vt[j-1].Ls + vt[j-1].Rs_r + A[i*nEq+(i-1)] ); 
				}
			}

			// side branch of right piriform fossa
			if( pfr_on == true ){
				if( j == nPos_pf-1 ){
					A[i*nEq+nS_pfr] = A[i*nEq+(i+1)];
				}
				else if( j == nPos_pf ){
					A[i*nEq+nS_pfr] = -1.0 * ( 2*fs*vt[j-1].Ls + vt[j-1].Rs_r + A[i*nEq+(i-1)] ); 
				}
			}
		}

#ifdef CHECK
		//fprintf( fp_log, "%d\t", j );
		//for( int k = 0; k < nEq; k++ ){
		//	fprintf( fp_log, "%e\t", A[i*nEq+k] );
		//}
		//fprintf( fp_log, "%e\n", b[i] );

		if( j == 0 ){
			fprintf( fp_log, "%d\t%f\t%f\t%f\n", j,  A[i*nEq + i], A[i*nEq + i+1], b[i]);
			//printf("%d\t%f\t%f\t%f\n", j,  A[i*nEq + i], A[i*nEq + i+1], b[i]);
		}
		else if( j > 0 && j < nvt ){
			fprintf( fp_log, "%d\t%f\t%f\t%f\t%f\n", j, A[i*nEq + i-1], A[i*nEq + i], A[i*nEq + i+1], b[i]);
			//printf( "%d\t%f\t%f\t%f\t%f\n", j, A[i*nEq + i-1], A[i*nEq + i], A[i*nEq + i+1], b[i]);
		}
		else if( j == nvt ){
			fprintf( fp_log, "%d\t%f\t%f\t%f\n", j, A[i*nEq + i-1], A[i*nEq + i], b[i]);
			//printf( "%d\t%f\t%f\t%f\n", j, A[i*nEq + i-1], A[i*nEq + i], b[i]);
		}
		//fprintf( fp_log, "\n" );
#endif
		i = i+1;
	}

#ifdef CHECK
	fclose( fp_log );
#endif

	//************************************************ linear equations for nasal tract *************************************************//
	// the nasal tract is between the section No. nPos_nt and nPos_nt+1  of vocal tract
	if( nt_on == true ){
		for( j = 0; j <= nnt; j++ ){
			if( j == 1 ){
				nIndx = nS_vt + nPos_nt - 1; 

				A[i*nEq+nIndx]   = 1.0/( 2 * fs * vt[nPos_nt-1].Cp + vt[nPos_nt-1].Yw);
				A[i*nEq+nIndx+1] = -1.0 * ( 2 * fs * vt[nPos_nt-1].Ls + vt[nPos_nt-1].Rs_r + A[i*nEq+nIndx] );

				h = 2 * fs * ( vt[nPos_nt-1].Ls + nt[j].Ls ) + vt[nPos_nt-1].Rs_r + nt[j].Rs_l;

				A[i*nEq+(i+1)] = 1.0/( 2 * fs * nt[j].Cp + nt[j].Yw); 
				A[i*nEq+i]     = -1.0 * (h + A[i*nEq+nIndx] +  A[i*nEq+(i+1)] );

				nt[j].constP = nt[j].Qcp - nt[j].Ud - nt[j].Yw * ( nt[j].Qlw - nt[j].Vcw );
				nt[j].constU = vt[nPos_nt-1].Qls_r + nt[j].Qls_l;

				b[i] = A[i*nEq+(i+1)] * nt[j].constP - A[i*nEq+nIndx] * vt[nPos_nt-1].constP - nt[j].constU;
			}
			else if( j == nnt ){
				h = 2 * fs * nt[j-1].Ls + nt[j-1].Rs_r;

				A[i*nEq+(i-1)] = 1.0 / ( 2 * fs * nt[j-1].Cp + nt[j-1].Yw);
				A[i*nEq+i]     = -1.0 * ( h + A[i*nEq+(i-1)] + 1.0/(Grn + Srn/(2*fs)) );

				constPn = -Vrn;
				constUn = nt[j-1].Qls_r;

				b[i] = 1.0/(Grn + Srn/(2*fs)) * constPn  - A[i*nEq+(i-1)] * nt[j-1].constP - constUn;
			}
			else{
				h = 2 * fs * (nt[j-1].Ls + nt[j].Ls) + (nt[j-1].Rs_r + nt[j].Rs_l); 

				A[i*nEq + (i-1)] = 1/ ( 2 * fs * nt[j-1].Cp + nt[j-1].Yw);
				A[i*nEq + (i+1)] = 1/ ( 2 * fs * nt[j].Cp   + nt[j].Yw);
				A[i*nEq + i]     = -1.0 * ( h + A[i*nEq + (i-1)] + A[i*nEq + (i+1)] );

				nt[j].constP = nt[j].Qcp - nt[j].Ud - nt[j].Yw *( nt[j].Qlw - nt[j].Vcw );
				nt[j].constU = nt[j-1].Qls_r + nt[j].Qls_l;

				b[i] =  A[i*nEq + (i+1)] * nt[j].constP -  A[i*nEq + (i-1)] * nt[j-1].constP - nt[j].constU;
			}

			i = i+1;
		}
	}

	//************************************ linear equations for left branch of piriform fossa ***************************************//
	double tmp;
	if( pfl_on == true ){
		for( j = 0; j < npfl; j++ ){
			if( j == 1 ){
				nIndx = nS_vt + nPos_pf - 1;

				A[i*nEq+nIndx]   = 1.0/( 2 * fs * vt[nPos_pf-1].Cp + vt[nPos_pf-1].Yw);
				A[i*nEq+nIndx+1] = -1.0 * ( 2 * fs * vt[nPos_pf-1].Ls + vt[nPos_pf-1].Rs_r + A[i*nEq+nIndx] );

				h = 2 * fs * ( vt[nPos_pf-1].Ls + pfl[j].Ls ) + vt[nPos_pf-1].Rs_r + pfl[j].Rs_l;

				A[i*nEq+(i+1)] = 1.0/( 2 * fs * pfl[j].Cp + pfl[j].Yw); 
				A[i*nEq+i]     = -1.0 * ( h + A[i*nEq+nIndx] +  A[i*nEq+(i+1)] );

				pfl[j].constP = pfl[j].Qcp - pfl[j].Ud - pfl[j].Yw * ( pfl[j].Qlw - pfl[j].Vcw );
				pfl[j].constU = vt[nPos_pf-1].Qls_r + pfl[j].Qls_l;

				b[i] = A[i*nEq+(i+1)] * pfl[j].constP - A[i*nEq+nIndx] * vt[nPos_pf-1].constP - pfl[j].constU;
			}
			else if( j == npfl-1 ){
				h = 2 * fs * (pfl[j-1].Ls + pfl[j].Ls) + (pfl[j-1].Rs_r + pfl[j].Rs_l); 

				A[i*nEq + (i-1)] = 1/ ( 2 * fs * pfl[j-1].Cp + pfl[j-1].Yw );
				tmp              = 1/ ( 2 * fs * pfl[j].Cp   + pfl[j].Yw );
				A[i*nEq + i]     = -1.0 * ( h + A[i*nEq + (i-1)] + tmp );

				pfl[j].constP = pfl[j].Qcp - pfl[j].Ud - pfl[j].Yw *( pfl[j].Qlw - pfl[j].Vcw );
				pfl[j].constU = pfl[j-1].Qls_r + pfl[j].Qls_l;

				b[i] = tmp * pfl[j].constP -  A[i*nEq + (i-1)] * pfl[j-1].constP - pfl[j].constU;
			}
			else{
				h = 2 * fs * (pfl[j-1].Ls + pfl[j].Ls) + (pfl[j-1].Rs_r + pfl[j].Rs_l); 

				A[i*nEq + (i-1)] = 1/ ( 2 * fs * pfl[j-1].Cp + pfl[j-1].Yw );
				A[i*nEq + (i+1)] = 1/ ( 2 * fs * pfl[j].Cp   + pfl[j].Yw );
				A[i*nEq + i]     = -1.0 * ( h + A[i*nEq + (i-1)] + A[i*nEq + (i+1)] );

				pfl[j].constP = pfl[j].Qcp - pfl[j].Ud - pfl[j].Yw *( pfl[j].Qlw - pfl[j].Vcw );
				pfl[j].constU = pfl[j-1].Qls_r + pfl[j].Qls_l;
				
				b[i] =  A[i*nEq + (i+1)] * pfl[j].constP -  A[i*nEq + (i-1)] * pfl[j-1].constP - pfl[j].constU;
			}
			i = i+1;
		}
	}

	//*********************************** linear equations for right branch of piriform fossa ***************************************//
	if( pfr_on == true ){
		for( j = 0; j < npfr; j++ ){
			if( j == 1 ){
				nIndx = nS_vt + nPos_pf - 1; 

				A[i*nEq+nIndx]   = 1.0/( 2 * fs * vt[nPos_pf-1].Cp + vt[nPos_pf-1].Yw);
				A[i*nEq+nIndx+1] = -1.0 * ( 2 * fs * vt[nPos_pf-1].Ls + vt[nPos_pf-1].Rs_r + A[i*nEq+nIndx] );

				h = 2 * fs * ( vt[nPos_pf-1].Ls + pfr[j].Ls ) + vt[nPos_pf-1].Rs_r + pfr[j].Rs_l;
				
				A[i*nEq+(i+1)] = 1.0/( 2 * fs * pfr[j].Cp + pfr[j].Yw);
				A[i*nEq+i]     = -1.0 * (h + A[i*nEq+nIndx] +  A[i*nEq+(i+1)] );

				pfr[j].constP = pfr[j].Qcp - pfr[j].Ud - pfr[j].Yw * ( pfr[j].Qlw - pfr[j].Vcw );
				pfr[j].constU = vt[nPos_pf-1].Qls_r + pfr[j].Qls_l;
				
				b[i] = A[i*nEq+(i+1)] * pfr[j].constP - A[i*nEq+nIndx] * vt[nPos_pf-1].constP - pfr[j].constU;
			}else if( j == npfr-1 ){
				h = 2 * fs * (pfr[j-1].Ls + pfr[j].Ls) + (pfr[j-1].Rs_r + pfr[j].Rs_l); 
				
				A[i*nEq + (i-1)] = 1/ ( 2 * fs * pfr[j-1].Cp + pfr[j-1].Yw);
				tmp              = 1/ ( 2 * fs * pfr[j].Cp   + pfr[j].Yw);
				A[i*nEq + i]     = -1.0 * ( h + A[i*nEq + (i-1)] + tmp );

				pfr[j].constP = pfr[j].Qcp - pfr[j].Ud - pfr[j].Yw *( pfr[j].Qlw - pfr[j].Vcw );
				pfr[j].constU = pfr[j-1].Qls_r + pfr[j].Qls_l  + pfr[j].Ns[0];

				b[i] =  tmp * pfr[j].constP -  A[i*nEq + (i-1)] * pfr[j-1].constP - pfr[j].constU;
			}else{
				h = 2 * fs * (pfr[j-1].Ls + pfr[j].Ls) + (pfr[j-1].Rs_r + pfr[j].Rs_l); 
				
				A[i*nEq + (i-1)] = 1/ (2*fs*pfr[j-1].Cp + pfr[j-1].Yw);
				A[i*nEq + (i+1)] = 1/ (2*fs*pfr[j].Cp   + pfr[j].Yw);
				A[i*nEq + i]     = -1.0 * ( h + A[i*nEq + (i-1)] + A[i*nEq + (i+1)] );

				pfr[j].constP = pfr[j].Qcp - pfr[j].Ud - pfr[j].Yw *( pfr[j].Qlw - pfr[j].Vcw );
				pfr[j].constU = pfr[j-1].Qls_r + pfr[j].Qls_l  + pfr[j].Ns[0];

				b[i] =  A[i*nEq + (i+1)] * pfr[j].constP -  A[i*nEq + (i-1)] * pfr[j-1].constP - pfr[j].constU;
			}
			i = i+1;
		}
	}
	return;
}

void CVT_sim::clearResidue( Tube *pt, const int nTube )
{
	int i;

	for( i = 0; i < nTube; i++ ){
		pt[i].Qls_l = 0.0;
		pt[i].Qls_r = 0.0;
		pt[i].Qcp   = 0.0;

		if( pt[i].bRigid == false ){
			pt[i].Qlw   = 0.0;
			pt[i].Vcw   = 0.0;
		}
	}

	return;
}

void CVT_sim::clearUP( Tube *pt, const int nTube )
{
	int i;
	for( i = 0; i < nTube; i++ ){
		pt[i].Ui = 0;
		pt[i].Uo = 0;
		pt[i].P  = 0;
		pt[i].Ud = 0;

		if( pt[i].bRigid == false ){
			pt[i].u3 = 0;
		}	
	}
	return;
}

void CVT_sim::clearNoiseSource( Tube *pt, const int nTube )
{
	int i, j;
	for( i = 0; i < nTube; i++ ){
		for( j = 0; j < 3; j++ ){
			pt[i].Ns[j] = 0.0;
			pt[i].Pn[j] = 0.0;
			pt[i].Un[j] = 0.0;
			pt[i].ru[j] = 0.0;
			pt[j].rp[j] = 0.0;
			pt[i].Amp_u[j] = 0.0;
			pt[i].Amp_p[j] = 0.0;
		}
	}
	return;
}

void CVT_sim::refresh_ax( double *a, double *x, double *a_nv )
// refresh the cross-sectional area and length of the short tubes comprise the vocal tract
// a     cross-section area function of vt
// x     length for short tubes of vt
// a_nv  cross-sectinal area of the coupling part between vocal tract and nasal tract
{
	int i;
	double r;

	//// refresh the cross-sectional area and length of glottis
	//r = 2.0;
	//for( i = 0; i < ngl; i++ ){
	//gl[i].a[1] = gl[i].a[0];
	//gl[i].a[0] = nonZero( ag[i] );

	//gl[i].s = sqrt( gl[i].a[0]/pi ) * 2 * pi * r;
	//}

	// refresh the cross-sectional area and length of vocal tract
	r = 2.0;
	for( i = 0; i < nvt; i++ ){
		vt[i].a[1] = vt[i].a[0];
		vt[i].a[0] = nonZero( a[i] );
		vt[i].s = sqrt( vt[i].a[0]/pi ) * 2 * pi * r;
		vt[i].x[1] = vt[i].x[0];vt[i].x[0] = x[i];
	}

	// refresh the coupling between nasal tract and vocal tract
	if( nt_on == true ){
		r = 5.0;
		for( i = 0; i < nvc; i++ ){
			nt[i].a[1] = nt[i].a[0];
			nt[i].a[0] = nonZero( a_nv[i] );
			nt[i].s = sqrt( nt[i].a[0]/pi ) * 2 * pi * r;
		}
	}

	return;
}

void CVT_sim::refresh( )
//Function: refresh the redisure of the vocal tract, nasal tract, glottis, and noise source etc.
{
	// refresh the residues within the vocal system
	Qls_g = 4 * fs * Lg * vt[0].Ui - Qls_g;
	refresh_residue( vt, nvt );
	if( nt_on == true )
		refresh_residue( nt, nnt );
	if( tr_on == true )
		refresh_residue( tr, ntr );
	if( pfl_on == true )
		refresh_residue( pfl, npfl );
	if( pfr_on == true )
		refresh_residue( pfr, npfr );

	// refresh the residues of radiation
	Pl = ( vt[nvt-1].Uo + constPl ) / ( Grl + Srl/(2*fs) );
	Vrl = Srl * Pl/fs + Vrl;

	if( nt_on == true ){
		Pn = ( nt[nnt-1].Uo + constPn ) / ( Grn + Srn/(2*fs) );
		Vrn = Srn * Pn/fs + Vrn; 
	}

	// refresh the noise source;
	// refresh_noise( vt, nvt );

	return;
}


void CVT_sim::refresh_noise( Tube *pt, const int nTube )
//Function: refresh the noise source along the vocal tract
{ 
	int i;
	for( i = 0; i < nTube; i++ ){
		pt[i].rp[2] = pt[i].rp[1]; 
		pt[i].rp[1] = pt[i].rp[0]; 
		pt[i].ru[2] = pt[i].ru[1]; 
		pt[i].ru[1] = pt[i].ru[0];

		pt[i].Amp_p[2] = pt[i].Amp_p[1]; 
		pt[i].Amp_p[1] = pt[i].Amp_p[0];
		pt[i].Amp_u[2] = pt[i].Amp_u[1];
		pt[i].Amp_u[1] = pt[i].Amp_u[0];

		pt[i].Un[2] = pt[i].Un[1]; 
		pt[i].Un[1] = pt[i].Un[0];
		pt[i].Pn[2] = pt[i].Pn[1];
		pt[i].Pn[1] = pt[i].Pn[0];
	}
	return;
}

void CVT_sim::refresh_residue( Tube *pt, const int nTube )
// function: refresh P, u3 and corresponding residue within the short tubes
{
	int i;

#ifdef CHECK_RE
	FILE *fp_re;
	fp_re = fopen( "./residue_c.txt", "w" );
#endif


	for( i = 0; i< nTube; i++ ){
		if( pt[i].bRigid == false ){
			pt[i].P  = ( pt[i].Ui - pt[i].Uo + pt[i].constP ) / ( 2 * fs * pt[i].Cp + pt[i].Yw );
		}else{
			pt[i].P  = ( pt[i].Ui - pt[i].Uo + pt[i].constP ) / ( 2 * fs * pt[i].Cp );
		}
		pt[i].Qls_l = 4 * fs * pt[i].Ls * pt[i].Ui - pt[i].Qls_l;
		pt[i].Qls_r = 4 * fs * pt[i].Ls * pt[i].Uo - pt[i].Qls_r;
		pt[i].Qcp   = 4 * fs * pt[i].Cp * pt[i].P  - pt[i].Qcp;



		if( pt[i].bRigid == false ){
			pt[i].u3  = pt[i].Yw * ( pt[i].P - pt[i].Vcw + pt[i].Qlw );
			pt[i].Qlw = 4 * fs * pt[i].Lw * pt[i].u3 - pt[i].Qlw;
			pt[i].Vcw = pt[i].u3 / (fs * pt[i].Cw)   + pt[i].Vcw;
		}

#ifdef CHECK_RE
		fprintf( fp_re, "%d\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n", i, pt[i].P, pt[i].Qls_l, pt[i].Qls_r, pt[i].Qcp, pt[i].u3, pt[i].Qlw, pt[i].Vcw );
#endif
	}

#ifdef CHECK_RE
	fclose( fp_re );
#endif


	return;
}

void CVT_sim::evalUdc( const double ag, const double Ps )
// function: evaluate the dc flow in the vocal tract
{
	double Rg_v;
	double Rc_v;
	double rv;
	double rk;

	// calculate the total viscous resistane of the two constrictions
	Rg_v = Krg / pow(ag,3);
	Rc_v = Krs * Ac.x / pow(Ac.a,2);

	// calculate the coefficients for the first and  second order term of equation
	// rk*Udc^2 + rv*Udc - Ps = 0; 
	// for detail please referenc "Acoustic Phonetics", Page 29
	rv = Rg_v + Rc_v;
	rk = Kb * ( 1.0/pow(ag,2) + 1.0/pow(Ac.a,2) ) / 2;

	// calculate the dc flow in the vocal tract
	Udc = ( -rv + sqrt( pow(rv,2) + 4*rk* abs(Ps) ) ) / (2*rk);

	Udc = nonZero( Udc );

	return;
}

void CVT_sim::findConstriction( )
// Function: determine the location of maximal constriction along the vocal tract,
// and the location of monopole and dipole sources
{
	int i;
	int nPos;
	double l = 0.0;
	double a_min = 100.0;

	for( i = 0; i < nvt; i++ ){
		if( vt[i].a[0] < a_min ){
			a_min = vt[i].a[0];
			nPos  = i;
		}
	}
	Ac.nPos = nPos;
	Ac.a    = vt[nPos].a[0];
	Ac.x    = vt[nPos].x[0];

	// determine the position of the monople flow source
	for( i = nPos;  i < nvt; i++ ){
		if( (vt[i].a[0] - vt[nPos].a[0]) > 0.1e-4 ){
			Ac.nPos_u = i+1;
			break;
		}
	}

	// determine the distance from outlet of constriction to the obstacle
	if( Ac.nPos_u >= nPos_th ){
		Ac.dist  = 0.5e-2;
		Ac.eta_n = 0.5;
		if( Ac.nPos_u == nvt ){
			Ac.nPos_p = Ac.nPos_u;
		}else{
			l = 0;
			for( i = Ac.nPos_u-1; i < nvt; i++ ){
				l = l + vt[i].x[0];
				if( l >= 0.5e-2 ){
					Ac.nPos_p = i+1;
					break;
				}
			}if( l < 0.5e-2 ){
				Ac.nPos_p = nvt;
			}
		}
	}else{
		l = 0;
		for( i = Ac.nPos_u-1; i < nPos_th; i++ ){
			l = l + vt[i].x[0];
			if( l >= 0.5e-2 ){
				// obstacle is the vocal tract wall
				Ac.dist   = 0.5e-2;
				Ac.eta_n  = 0.5;
				Ac.nPos_p = i+1;
				break;
			}
		}
		if( l < 0.5e-2 ){
			// obstacle is the teeth
			Ac.dist  = l;
			Ac.eta_n = 1.0;
			Ac.nPos_p = nPos_th;
		}
	}

	return;
}

void CVT_sim::evalNoiseSource( Tube *pt, int nTube )
{
	int i = 0, j = 0;
	double d;            // characteristic dimension of the tube 
	double u;            // partical vecolcity of air stream in the tube
	double fc_u = 1100;  // cutoff frequency for flow source
	double fc_p;         // cutoff frequency for pressure source
	
	// calculate the monople and dipole source for noise
	for( i = 0; i < nTube; i++ ){
		pt[i].ru[0] = pt[i].rp[0] = 0;
		pt[i].Amp_u[0] = 0;
		pt[i].Amp_p[0] = 0;
	}
	u = Udc / Ac.a;
	d = sqrt( 4 * pt[Ac.nPos].a[0] / pi );

	pt[Ac.nPos].Re = u * d / mu_air;
	fc_p = 0.4 * u/d;
	
	if( pt[Ac.nPos].Re > 1800 ){
		pt[Ac.nPos_u].ru[0] = pt[Ac.nPos_p].rp[0] = randn();
		pt[Ac.nPos_u].Amp_u[0] = alpha_n * ( pow(pt[i].Re,2) - pow(1.8e3, 2) ); 
		pt[Ac.nPos_p].Amp_p[0] = alpha_n * ( pow(pt[i].Re,2) - pow(1.8e3, 2) ) * Ac.eta_n * exp( -1.0 * Ac.dist/tau_n );
	}

	for( i = 0; i < nTube; i++ ){
		// low pass filtering of the amplitude of random noise with single pole low pass filter
		// where a[0] = 0.15; b = {1.0, -0.85}
		pt[i].Amp_u[0] = 0.15 * pt[i].Amp_u[0] + 0.85 * pt[i].Amp_u[1];
		pt[i].Amp_u[1] = pt[i].Amp_u[0];

		pt[i].Amp_p[0] = 0.15 * pt[i].Amp_p[0] + 0.85 * pt[i].Amp_p[1];
		pt[i].Amp_p[1] = pt[i].Amp_p[0];

		// low pass filter of the random noise;
		dsp.butterworth( 2, &fc_u, 1, fs, 1, Au, Bu );
		for( j = 0; j < 3; j++ ){
			Bu[j] = pt[i].Amp_u[0] * Bu[i];
		}
		dsp.filter( Au, 2, Bu, 2, pt[i].ru, pt[i].Un );

		//pt[i].Un[0] = pt[i].Amp_u[0] * pt[i].Un[0];
		dsp.butterworth( 2, &fc_p, 1, fs, 1, Ap, Bp );
		for( j = 0; j < 3; j++ ){
			Bp[j] = pt[i].Amp_p[0] * Bp[i];
		}
		dsp.filter( Ap, 2, Bp, 2, pt[i].rp, pt[i].Pn );
		//pt[i].Pn[0] = pt[i].Amp_p[0] * pt[i].Pn[0];
	}
	
	return;
}

void CVT_sim::turnOntract( bool bOn, int nType )
{
	switch( nType ){
		case NT:
			nt_on = bOn;
			break;
		case TR:
			tr_on = bOn;
			break;
		case PFL: 
			pfl_on = bOn;
			break;
		case PFR:
			pfr_on = bOn; 
			break;
	}
	return;
}

void CVT_sim::loadconfig( char *fn, int nType )
{
	int nTube;int i;
	double a, x;  // cross-sectional area and length
	double r;     // the ration between area and perimeter
	bool bRigid;
	FILE *fp;
	
	fp = fopen(fn, "r" );
	if( fp == NULL ){
		cout<<"Cann't open file: "<<fn<<endl;exit( 1 );
	}
	fscanf( fp, "%d", &nTube );

	if( nType == NT ){
		// load the # of short tubes in the nasal tract which change their cross-sectional areas 
		// when coupling with the vocal tract
		fscanf( fp, "%d", &nPos_nt);
		fscanf( fp, "%d", &nvc );
		fscanf( fp, "%d", &nPos_sin );  // location of sinus
	}
	else{
		nvc = 0;
	}

	if( nType == VT ){
		// load the location where nasal tract and piriform fossa coupling with vocal tract
		fscanf( fp, "%d", &nPos_pf );
		fscanf( fp, "%d", &nPos_nt );
	}

	Tube *pt = new Tube[nTube];

	switch( nType ){
		case TR:
			tr = pt;
			ntr = nTube; 
			r = 1.0; 
			bRigid = true;  
			break;
		case VT: 
			vt = pt; 
			nvt = nTube; 
			r = 2.0; 
			bRigid = false;
			break;
		case NT: nt = pt;
			nnt = nTube; 
			r = 5.0; 
			bRigid = false; 
			break;
		case PFL: 
			pfl = pt; 
			npfl = nTube; 
			r = 1.2; 
			bRigid = false;
			break;
		case PFR: 
			pfr = pt; 
			npfr = nTube; 
			r = 1.2; 
			bRigid = false; 
			break; 
	}

	for( i = 0; i < nTube; i++ ){
		fscanf( fp, "%lf %lf", &a, &x );
		pt[i].a[0] = pt[i].a[1] = nonZero(a);
		pt[i].x[0] = pt[i].x[1] = x;
		pt[i].s = sqrt( pt[i].a[0]/pi ) * 2 * pi * r;
		pt[i].bRigid = bRigid;
	}
	fclose( fp );
	return;
} 

double CVT_sim::nonZero( double x )
{ 
	double r;
	if( x < AREA_MIN )
		r = AREA_MIN;
	else
		r = x;

	return r;
}


void CVT_sim::acoust_elm_t( Tube *pt, int nTube, double m, double b, double k )
//Function: calcualte the acoustic elemenst Rs, Ls, Cp, Lw, Cw, Rw, etc.
{ 
	int i;
	double x_mp_s;

	Klw = m;
	Krw = b;
	Kcw = 1.0 / k;

	for( i = 0; i < nTube; i++ ){
		pt[i].Ls = 0.5 * Kls *  pt[i].x[0] / pt[i].a[0] ;
		pt[i].Rs = 0.5 * Krs * pt[i].x[0] / pow( pt[i].a[0], 2 );

		pt[i].Rs_l = pt[i].Rs;
        pt[i].Rs_r = pt[i].Rs;

        pt[i].Cp = Kcp * pt[i].x[0] * pt[i].a[0];
		pt[i].Ud = fs * (pt[i].x[0]*pt[i].a[0] - pt[i].x[1]*pt[i].a[1]);

		// calculate the acoustic introduced by wall vibration
		if( pt[i].bRigid == false ){
			x_mp_s =  pt[i].x[0] * pow(pt[i].s,2);
			pt[i].Lw = Klw / x_mp_s;
			pt[i].Rw = Krw / x_mp_s;
			pt[i].Cw = Kcw * x_mp_s;

			pt[i].Yw = 1 / ( 2 * fs * pt[i].Lw + pt[i].Rw + 1.0/(2*fs*pt[i].Cw) );
		}
        else{
            pt[i].Lw = 0.0;
            pt[i].Rw = 0.0;
            pt[i].Cw = 0.0;
            
            pt[i].Yw = 0.0;
        }
	}

	return;
}

void CVT_sim::addBernoulli( Tube *pt, int nTube )
{
	int i;

	for( i = 1; i < nTube-1; i++ ){

		if( pt[i].a[0] < pt[i-1].a[0] ){
			pt[i].Rs_l = pt[i].Rs_l + abs(pt[i].Ui) * ro_air * 1.42 / (2*pow(pt[i].a[0],2));
		}

		if( pt[i].a[0] > pt[i+1].a[0] ){
			pt[i].Rs_r = pt[i].Rs_r - abs(pt[i].Uo) * ro_air * 1.42 / (2*pow(pt[i].a[0],2));
		}
	}

	return;
}

void CVT_sim::modifyNT_sin( int nPos_sin, double r_sin, double l_sin, double c_sin, double y_sin )
// modify the acoustic elements of nasal tract by considering nasal sinus
{
	if( bSubstitute == true ){
		nt[nPos_sin].Rw = r_sin;
		nt[nPos_sin].Lw = l_sin;
		nt[nPos_sin].Cw = c_sin;
		nt[nPos_sin].Yw = y_sin;
		nt[nPos_sin].Cp = 0;
	}
	else{
		nt[nPos_sin].Yw = nt[nPos_sin].Yw + y_sin;
	}
	return;
}

void CVT_sim::glottalImpedance(double *ag, int nTube)
{
	int i;

	Rg = 0;
	Lg = 0;

	for( i = 0; i < nTube; i++ ){
		Rg = Krg / pow(ag[i],3) + Kdg * Udc / pow(ag[i],2);
		Lg = Klg / ag[i];
	}

	return;
}

void CVT_sim::radiationImpedance( double a, double &gr, double &sr )
{ 
	gr = Kgr * a;
	sr = Ksr * sqrt(a);
	
	return;
}

void CVT_sim::sinusImpedance(const double vol, const double len, const double apt, double &R_sin, double &L_sin, double &C_sin, double &Y_sin)
// calculate the Impedance of nasal sinus by using Helmholtz resonator model
{
	R_sin = Krs * len / pow(apt,1.5);
	L_sin = Kls * len / apt;
	C_sin = Kcp * vol;
	Y_sin = 1 / ( 2*L_sin + R_sin + 1.0/(2.0*C_sin) );

	return;
}

double CVT_sim::glottalArea_dy(double T0, double t, double OQ, double SQ, double Ap )
//function: caclualte dynamic component of the glottal area according
//          to the input parameters
// T0   peroid of vibration
// t    time instance
// OQ   open quotient
// SQ   speed quotient
// Ap   miximum openning of this period
{
	double t1;   // time instance where glottal area achieve the maximal openning
	double t2;   // time instance where glottal area become zero from nonzero
	double K;
	double Ag_dy;

	t2 = OQ * T0;
	t1 = SQ/(1.0+SQ) * t2;
	K = 1.0 / (1.0 - cos(pi*(t2-t1)/t1));
	if( t <= t1 ){
		Ag_dy = 0.5 * Ap * ( 1 - cos(pi*t/t1) );
	}
	else if( t > t1 && t <= t2 ){
		Ag_dy = Ap * ( 1 - K + K * cos(pi*(t-t1)/t1) );
	}
	else{
		Ag_dy = 0.0;
	}
	Ag_dy = nonZero( Ag_dy );
	
	return Ag_dy;
}

double CVT_sim::glottalArea_dc( )
//function: calculate the dc component of the glottal area
//need to be implement
{
	double Ag_dc;

	Ag_dc = 0.0;

	return Ag_dc;
}

double CVT_sim::randn( )
// Function: generate Gaussian distributed random number
{
	int i;
	double x;
	double t[100];

	// mean and variance of equal-distribution span in [a b], a = 0.0, b = 1.0;
	// u = (a+b)/2; var = (b-a)^2/12;
	double u   = 0.5;
	double var = 1.0/12;

	// generate equal distributed rand number between 0 and 1.0
	x = 0.0;
	for( i = 0; i < 100; i++ ){
		t[i] = (double)rand() / RAND_MAX;
		x = x + t[i];
	}

	// generate standard normal distributed rand number using central-limit theory    
	x = ( x - 100 * u ) / sqrt( 100 * var );
	return x;
}

//void CVT_sim::filter( double *x, double *y, double *a, double *b )
//{
//y[2] = y[1];
//y[1] = y[0];
//
//y[0] = b[0] * x[0] + b[1] * x[1] + b[2] * x[2] - a[1] * y[1] - a[2] * y[2];
//
//x[2] = x[1];
//x[1] = x[0];
//
//return;
//}

//double CVT_sim::frictionNoise( int Entry, double inMem[4][3], double outMem[4][3], double *mem )
//// function: calculate friction noise source within the vocal tract 
//// generat filtered noise to mimic turbulent noisee in the vocaltract.
//// Computer generated random numbers undergo two filters: 
//// a high-pass to cutoff low frequency components below 200Hz,
//// and a low-pass (-6dB/oct) to shape noise at high frequencies//// The low-pass filter is a simple first order integrater with alpha = 0.8.
//// The 3rd order high-pass filter is used for the low cut.
//// The cutoff frequency is 200Hz when the sampling frequency is 10kHz//// Arguments:
//// Entry           if == 0, initialization//// inMem        fixed working array ( input memories )
//// outMem        fixed working array ( output memories )
//// *mem        fixed memory fr low pss filter//{
//int i, j;
//double s;
//
//static double a[2][3] = { {1.0, -1.8977, 0.9069}, {1.0, -0.9067, 0.0} };
//static double b[2][3] = { {0.9546, -1.9093, 0.9546}, {0.9974, -0.9974, 0.0} };
//
//static double max, alpha = 0.8;
//
//if( Entry == 0 ){
//max = 2.0 / RAND_MAX; // RAND_MAX 


//for( j = 0; j < 3; j++ ){
//for( i = 0; i < 2; i++ ){
//if( i == 0 ){
//inMem[i][0] = max * rand() - 1;
//}
//else{
//inMem[i][0] = outMem[i-1][0];
//}
//filter( inMem[i], outMem[i], a[i], b[i] );
//}
//s = outMem[1][0];
//s = s + alpha * (*mem);
//*mem = s;
//}
//}
//
//for( i = 0; i < 2; i++ ){
//if( i == 0 ){
//inMem[i][0] = max * rand() - 1;
//}
//else{
//inMem[i][0] = outMem[i-1][0];
//}
//filter( inMem[i], outMem[i], a[i], b[i] );
//}
//
//s = outMem[1][0];
//s = s + alpha * (*mem);
//*mem = s;
//
//return s;
//}


// There are three kind of turbulent noise sources:
// monopole source:  random velocity fluctuation in constriction ( volum source )
//                   velocity fluctionation in the jet leaves the constriction//                   the cutoff frequency of the monopole source is 1100Hz// dipole source:    air flow impinge on the obstacle  ( pressure source )
//                   the cut off frequency for dipole source is between 2500Hz and 6000 Hz