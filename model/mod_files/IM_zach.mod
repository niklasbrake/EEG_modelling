TITLE Cortical Muscarinic current

COMMENT
-----------------------------------------------------------------------------
	Muscarinic potassium current
	----------------------------

  - noninactivating K+ current responsible for the adaptation of firing rate 
    and the afterhyperpolarization (AHP) of cortical pyramidal cells

  - differential equations (HH-like)

  Model from Mainen, Z.M. and Sejnowski, T.J., 1996

  Written by Alain Destexhe, Laval University, 1996
-----------------------------------------------------------------------------
ENDCOMMENT


INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX imZ
	USEION k READ ek WRITE ik
        RANGE gkbar, m_inf, tau_m
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}


PARAMETER {
	v		(mV)
	celsius = 36    (degC)
	ek	= -90	(mV)
	gkbar	= 1e-6	(mho/cm2)
}



STATE {
	m
}

ASSIGNED {
	ik	(mA/cm2)
	m_inf
	tau_m	(ms)
	tadj
}

BREAKPOINT {
	SOLVE states METHOD euler
	ik = gkbar * m * (v - ek)
}

DERIVATIVE states { 
	evaluate_fct(v)

	m' = (m_inf - m) / tau_m
}

UNITSOFF
INITIAL {
:
:  The Q10 value is assumed to be 2.3
:
        tadj = 2.3 ^ ((celsius-23)/10)
	evaluate_fct(v)
	m = m_inf
}

PROCEDURE evaluate_fct(v(mV)) {  LOCAL a,b

	a =  1e-4 * (v+30) / ( 1 - exp(-(v+30)/9) )
	b = -1e-4 * (v+30) / ( 1 - exp( (v+30)/9) )

	tau_m = 1 / (a + b) / tadj
	m_inf = a / (a + b)
}
UNITSON


