#include <octave/oct.h>
#include <string.h>

#ifdef __cplusplus
extern "C"
{
#endif
#include "ngspice/sharedspice.h"
#ifdef __cplusplus
}  /* end extern "C" */
#endif

DEFUN_DLD(ngspice_command, args,,
          "Issue ngspice command\nUsage: ngspice_command(string array)\n")
{
    int nargin = args.length();

    if(nargin != 1) {
        print_usage();
    } else {
        charMatrix str_vect = args(0).char_matrix_value();
        Matrix m(str_vect.rows(), 1);
        
        for(octave_idx_type i = 0; i < str_vect.rows(); i++) {
            char *cmd_str = strdup(str_vect.row_as_string(i).c_str());
            
            if(cmd_str != (char*)NULL) {
                // octave_stdout << cmd_str << "\n";
                m(i, 0) = ngSpice_Command(cmd_str);
                free(cmd_str);
            } else {
                m(i, 0) = 0;
            }
        }
        return octave_value(m);
    }
    return octave_value_list();
}
