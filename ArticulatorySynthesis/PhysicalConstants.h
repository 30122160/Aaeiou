#ifndef _PHYSICAL_CONSTANTS_H_
#define _PHYSICAL_CONSTANTS_H_

#include <cmath>

const double pi = 4 * atan( 1.0 );     // 

//******************************* Physical constancts for transmition line vocal tract model ***********************************/
const double ro_air = 1.2929;  // air density (Unit: kg/m^3)
const double c_air  = 331.45;  // sound velocity in air (Unit: m/s)
const double mu_air = 1.86e-5; // viscous coefficient (Unit:N*s/m)

const double eta    = 1.4;     // adiabatic constant ??
const double cp     = 0.24e6;  // specific heat (Unit: cal/(kg*degree))
const double lamda  = 5.5e-3;  // head conduction (Unit: cal/(m*s*degree)) 

const double m_vt = 15;        // mass of vocal tract wall per square meter (Unit: kg/m^2)
const double b_vt = 1.6e4;     // machanial resistance of vocal tract wall (Unit: kg/(m^2*s))
const double k_vt = 3.0e2;    // stiffness of vocal tract wall (Unit: kg/(s^2))

const double m_tr = 5;         //
const double b_tr = 0.3e4;     // 
const double k_tr = 1.0e7;     // ??? the value need to be checked

const double H2O   = 9.8e3;    // pressure produced by water column with the height of 1 meter
const double Re_c  = 1800.0;   // critical Reynold number

const double extr_loss = 50;   // extra heat loss for the nasal tract  
// configuration for glottis
const double lg = 1.2e-2;     // lenght of vocal folds
const double dg = 0.3e-2;     // thickness of glottis
const double kg = 1.40;       // coefficient related to the kinetic component of glottal resistance

const double vol_lung = 4.0;    // the volume of lung (Unit: m^3)
const double vol_sin  = 2.5e-2; // the volume of nasal side cavity (Unit: m^3)
const double apt_sin  = 1.0e-5; // aperture of sinus canal (Unit:m^2)
const double len_sin  = 0.5e-2; // length of sinus canal

const double Rlung = 15/sqrt(600.0); 

const double noiseAmp = 100.0;  // amplitude bound for noise source
// constancts for generating turbulent sources
// see "Simulation of losses due to turbulence in the time varying vocal system", 
// IEEE Trans. on Audio, speech, and language processing Vol15, No4
const double alpha_n = 4.0e-6;                 //  (Unit: N/m^2 )
const double beta_n  = 5.0e-8;                 //  (Unit: m^5 / (N*s) ) 
const double tau_n   = 1.23;                   // reference for generating dipole source ( voltage source )

#endif