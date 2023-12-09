TITLE Delayer rectifier

COMMENT
-----------------------------------------------------------------------------
	"delayer-rectifier" K current for action potentials
	---------------------------------------------------

  - potassium current, voltage-dependent
  - iterative equations

  Model of IKd for hippocampal pyramidal cells, from
  Traub & Miles, Neuronal Networks of the Hippocampus, Cambridge, 1991

  Written by Alain Destexhe, Laval University, 1996
-----------------------------------------------------------------------------
ENDCOMMENT


INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX ikdT
	USEION k READ ek WRITE ik
	RANGE gkbar, vtraub
	RANGE n_inf
	RANGE tau_n
	RANGE n_exp
}


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

PARAMETER {
	gkbar	= .005 	(mho/cm2)
	vtraub	= -55	(mV)		: adjusts threshold
	ek	= -90	(mV)
	celsius = 36    (degC)
	dt              (ms)
	v               (mV)
}

STATE {
	n
}

ASSIGNED {
	ik	(mA/cm2)
	n_inf
	tau_n
	n_exp
	tadj
}


BREAKPOINT {
	SOLVE states
	ik  = gkbar * n*n*n*n * (v - ek)
}


:DERIVATIVE states {
:	evaluate_fct(v)
:	n' = (n_inf - n) / tau_n
:}

PROCEDURE states() {	: exact when v held constant
	evaluate_fct(v)
	n = n + n_exp * (n_inf - n)
}

UNITSOFF
INITIAL {
:
:  Q10 was assumed to be 2.3 for both currents
:
:  original measurements at room temperature

	tadj = 3.0 ^ ((celsius-36)/ 10 )
	evaluate_fct(v)
	n = n_inf
}

PROCEDURE evaluate_fct(v(mV)) { LOCAL a,b,v2

	v2 = v - vtraub : convert to traub convention

	a = 0.032 * (15-v2) / ( exp((15-v2)/5) - 1)
	b = 0.5 * exp((10-v2)/40)

	tau_n = 1 / (a + b) / tadj
	n_inf = a / (a + b)

	n_exp = 1 - exp(-dt/tau_n)
}

UNITSON
