#include <stdio.h>
#include "hocdec.h"
extern int nrnmpi_myid;
extern int nrn_nobanner_;

extern void _IKd_traub_reg(void);
extern void _IM_zach_reg(void);
extern void _INa_traub_shifted_reg(void);
extern void _Leak_reg(void);

void modl_reg(){
  if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
    fprintf(stderr, "Additional mechanisms from files\n");

    fprintf(stderr," \"IKd_traub.mod\"");
    fprintf(stderr," \"IM_zach.mod\"");
    fprintf(stderr," \"INa_traub_shifted.mod\"");
    fprintf(stderr," \"Leak.mod\"");
    fprintf(stderr, "\n");
  }
  _IKd_traub_reg();
  _IM_zach_reg();
  _INa_traub_shifted_reg();
  _Leak_reg();
}
