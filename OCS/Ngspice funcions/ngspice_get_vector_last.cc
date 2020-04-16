#include <octave/oct.h>
#include <string.h>

#ifdef __cplusplus
extern "C"
{
#endif
#include "ngspice/sharedspice.h"
   /* Dvec flags. */
   enum dvec_flags {
      VF_REAL = (1 << 0),           /* The data is real. */
      VF_COMPLEX = (1 << 1),        /* The data is complex. */
      VF_ACCUM = (1 << 2),          /* writedata should save this vector. */
      VF_PLOT = (1 << 3),           /* writedata should incrementally plot it. */
      VF_PRINT = (1 << 4),          /* writedata should print this vector. */
      VF_MINGIVEN = (1 << 5),       /* The v_minsignal value is valid. */
      VF_MAXGIVEN = (1 << 6),       /* The v_maxsignal value is valid. */
      VF_PERMANENT = (1 << 7)       /* Don't garbage collect this vector. */
   };

#define isreal(v)   ((v)->v_flags & VF_REAL)
#define iscomplex(v)    ((v)->v_flags & VF_COMPLEX)

#ifdef __cplusplus
}  /* end extern "C" */
#endif

DEFUN_DLD(ngspice_get_vector_last, args, nargout,
          "Get last data of the named vector\n"
          "Usage: ngspice_get_vector(string array of vector names)\n")
{
   int nargin = args.length();
   octave_value_list retval;

   if(nargin != 1) {
      print_usage();
   } else {
      const charMatrix str_vect = args(0).char_matrix_value();

      for(octave_idx_type i = 0; (i < str_vect.rows()) && (i < nargout); i++) {
         char* const vect_name = strdup(str_vect.row_as_string(i).c_str());

         if(vect_name != (char*)NULL) {
            const pvector_info info = ngGet_Vec_Info(vect_name);

            // octave_stdout << vect_name << "\n";
            // octave_stdout << info << "\n";

            /* check if vector info could be retrieved */
            if(info != (pvector_info)NULL) {
               /* check if real or complex */
               if(isreal(info)) {
                  retval(i) = info->v_realdata[info->v_length - 1];
               } else if(iscomplex(info)) {
                  Complex c(info->v_compdata[info->v_length - 1].cx_real,
                            info->v_compdata[info->v_length - 1].cx_imag);
                  retval(i) = c;
               } else {
                  retval(i) = octave_value(0);
               }
            } else {
               warning("No vector with name %s", vect_name);
               retval(i) = octave_value(0);
            }
            /* free vect_name string */
            free(vect_name);
         } else {
            warning("Out of memory!");
            retval(i) = octave_value(0);
         }
      }
      return retval;
   }
   return octave_value_list();
}
