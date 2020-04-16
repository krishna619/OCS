#include <octave/oct.h>

#ifdef __cplusplus
extern "C"
{
#endif
#include "ngspice/sharedspice.h"
#ifdef __cplusplus
}  /* end extern "C" */
#endif

DEFUN_DLD(ngspice_init, args,,
          "Start NgSpice Simulator")
{
    ngSpice_Init((SendChar*)NULL, (SendStat*)NULL, (ControlledExit*)NULL,
                 (SendData*)NULL, (SendInitData*)NULL, (BGThreadRunning*)NULL, NULL);
    return octave_value_list();
}
