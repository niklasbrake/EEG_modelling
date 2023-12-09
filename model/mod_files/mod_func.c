#include <stdio.h>
#include "hocdec.h"
#define IMPORT extern __declspec(dllimport)
IMPORT int nrnmpi_myid, nrn_nobanner_;

extern void _IKd_traub_reg();
extern void _IM_zach_reg();
extern void _INa_traub_shifted_reg();
extern void _Leak_reg();

void modl_reg(){
	//nrn_mswindll_stdio(stdin, stdout, stderr);
    if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
	fprintf(stderr, "Additional mechanisms from files\n");

fprintf(stderr," IKd_traub.mod");
fprintf(stderr," IM_zach.mod");
fprintf(stderr," INa_traub_shifted.mod");
fprintf(stderr," Leak.mod");
fprintf(stderr, "\n");
    }
_IKd_traub_reg();
_IM_zach_reg();
_INa_traub_shifted_reg();
_Leak_reg();
}
