/*
  
  Copyright (C) 2009 Massimiliano Culpo
  
  This file is part of:
  OCS - A Circuit Simulator for Octave

  OCS is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program (see the file LICENSE); if not,
  see <http://www.gnu.org/licenses/>.
  
  author: Massimiliano Culpo <culpo@math.uni-wuppertal.de>
  
*/

#include <octave/oct.h>
#include <octave/parse.h>
#include <octave/str-vec.h>
#include <string>
#include <octave/pager.h>

using namespace std;

/* Print parameters */
void print_parameters(double rd,  double W,   double L,   double mu0, double Vth,
		      double Cox, double Cgs, double Cgd, double Cgb, double Csb,
		      double Cdb, double Tshift)
{
  octave_stdout << "PARAMETER TABLE: Simplified Shichman-Hodges MOS-FET\n\n";
  octave_stdout << "Name\t\tValue\n";
  octave_stdout << "rd\t\t"      << rd     << "\n" ;
  octave_stdout << "W\t\t"       << W      << "\n" ;
  octave_stdout << "L\t\t"       << L      << "\n" ;
  octave_stdout << "mu0\t\t"     << mu0    << "\n" ;
  octave_stdout << "Vth\t\t"     << Vth    << "\n" ;
  octave_stdout << "Cox\t\t"     << Cox    << "\n" ;
  octave_stdout << "Cgs\t\t"     << Cgs    << "\n" ;
  octave_stdout << "Cgd\t\t"     << Cgd    << "\n" ;
  octave_stdout << "Cgb\t\t"     << Cgb    << "\n" ;
  octave_stdout << "Csb\t\t"     << Csb    << "\n" ;
  octave_stdout << "Cdb\t\t"     << Cdb    << "\n" ;
  octave_stdout << "Tshift\t\t"  << Tshift << "\n" ;
  octave_stdout << "\n\n";
}

/* Print computed values */
void print_values(double gm, double gd,   double ids,    double didT,
		  double P,  double dPdT, double dPdvgs, double dPdvds)
{
  octave_stdout << "COMPUTED VALUES TABLE: Simplified Shichman-Hodges MOS-FET\n\n";
  octave_stdout << "Name\t\tValue\n";
  octave_stdout << "gm\t\t"      << gm     << "\n" ;
  octave_stdout << "gd\t\t"      << gd     << "\n" ;
  octave_stdout << "ids\t\t"     << ids    << "\n" ;
  octave_stdout << "didT\t\t"    << didT   << "\n" ;
  octave_stdout << "P\t\t"       << P      << "\n" ;
  octave_stdout << "dPdT\t\t"    << dPdT   << "\n" ;
  octave_stdout << "dPdvgs\t\t"  << dPdvgs << "\n" ;
  octave_stdout << "dPdvds\t\t"  << dPdvds << "\n" ;
  octave_stdout << "\n\n";
}

/* Set the parameters of the simplified Shichman-Hodges MOS-FET */
void set_parameters(ColumnVector parameters,  string_vector parameternames,
		    double *rd,  double* W,   double* L,   double* mu0, double* Vth,
		    double* Cox, double* Cgs, double* Cgd, double* Cgb, double* Csb,
		    double* Cdb, double* Tshift)
{
  octave_idx_type nnames = parameternames.numel ();
  octave_idx_type niter  = 0;
  //FIXME: it should be better to use Octave_map<string,value> for parameters
  while (niter < nnames)
    {
      if (parameternames[niter] == "rd")
	*rd = parameters(niter);
      else if (parameternames[niter] == "W")
	*W = parameters(niter);
      else if (parameternames[niter] == "L")
	*L = parameters(niter);
      else if (parameternames[niter] == "mu0")
	*mu0 = parameters(niter);
      else if (parameternames[niter] == "Vth")
	*Vth = parameters(niter);
      else if (parameternames[niter] == "Cox")
	*Cox = parameters(niter);
      else if (parameternames[niter] == "Cgs")
	*Cgs = parameters(niter);
      else if (parameternames[niter] == "Cgd")
	*Cgd = parameters(niter);
      else if (parameternames[niter] == "Cgb")
	*Cgb = parameters(niter);
      else if (parameternames[niter] == "Csb")
	*Csb = parameters(niter);
      else if (parameternames[niter] == "Cdb")
	*Cdb = parameters(niter);
      else if (parameternames[niter] == "Tshift")
	*Tshift = parameters(niter);
      else
	warning ("%s",(string("Mshichmanhodgesmosfet: unknown parameter").append (parameternames[niter])).c_str ());
    
      niter++;
    }
}

/* Compute values for n-mos model*/
void nmos(ColumnVector extvar, double mu0,   double Cox,     double W,
	  double L,            double Vth,   double rd,      double Tshift,
	  double *gm,          double *gd,   double *ids,    double *didT,
	  double *P,           double *dPdT, double *dPdvgs, double *dPdvds)
{
  double vg   = extvar(0); // V-gate
  double vs   = extvar(1); // V-source
  double vd   = extvar(2); // V-drain
  double vb   = extvar(3); // V-bulk
  double T    = extvar(4); // Temperature

  double k    = mu0*Cox*pow((T + Tshift)/300.0,-3.0/2.0)*W/L;
  double dkdT = mu0*Cox*W*(-3.0/2)*pow((T + Tshift)/300.0,-5.0/2.0 )*(1.0/300.0)/L;
  
  double vgs  = vg - vs;
  double vds  = vd - vs;

  if (vgs < Vth)
    {
      *gm   = 0;
      *gd   = 1/rd;
      *ids  = vds*(*gd);
      *didT = 0;
    }  
  else if ( ( (vgs-Vth)>= vds ) && (vds>=0))
    {
      *ids  = k*((vgs-Vth)*vds - pow(vds,2)/2 ) + vds/rd;
      *gm   = k*vds;
      *gd   = k*(vgs-Vth-vds) + (1/rd);
      *didT = dkdT*((vgs-Vth)*vds-(pow(vds,2))/2);
    }  
  else if (((vgs-Vth)>=(vds))&&(vds<0))
    {  
      *gm  = 0;
      *gd  = 1/rd;
      *ids = vds*(*gd);
      *didT= 0;
    }  
  else // (i.e. if 0 <= vgs-vth <= vds)
    {  
      *ids = (k/2)*pow((vgs-Vth),2) + vds/rd;
      *gm  = k*(vgs-Vth);
      *gd  = 1/rd;
      *didT= (dkdT/(2))*pow((vgs-Vth),2);
    }  

  *P       = -(*ids)*vds;
  *dPdT    = -(*didT)*vds;
  *dPdvgs  = -(*gm)*vds;
  *dPdvds  = -((*gd)*vds + (*ids));

}

/* Compute values for p-mos model*/
void pmos(ColumnVector extvar, double mu0,   double Cox,     double W,
	  double L,            double Vth,   double rd,      double Tshift,
	  double *gm,          double *gd,   double *ids,    double *didT,
	  double *P,           double *dPdT, double *dPdvgs, double *dPdvds)
{
  double vg   = extvar(0); // V-gate
  double vs   = extvar(1); // V-source
  double vd   = extvar(2); // V-drain
  double vb   = extvar(3); // V-bulk
  double T    = extvar(4); // Temperature

  double k    = - mu0*Cox*pow((T + Tshift)/300.0,-3.0/2.0)*W/L;
  double dkdT = - mu0*Cox*W*(-3.0/2.0)*pow((T + Tshift)/300.0,-5.0/2.0 )*(1.0/300.0)/L;

  double vgs  = vg - vs;
  double vds  = vd - vs;

  if (vgs > Vth)
    {
      *gm   = 0;
      *gd   = 1/rd;
      *ids  = vds*(*gd);
      *didT = 0;
    }  
  else if ( ( (vgs-Vth)<= vds ) && (vds<=0))
    {
      *ids  = k*((vgs-Vth)*vds - pow(vds,2)/2 ) + vds/rd;
      *gm   = k*vds;
      *gd   = k*(vgs-Vth-vds) + (1/rd);
      *didT = dkdT*((vgs-Vth)*vds-(pow(vds,2))/2);
    }  
  else if (((vgs-Vth)<=(vds))&&(vds>0))
    {  
      *gm  = 0;
      *gd  = 1/rd;
      *ids = vds*(*gd);
      *didT= 0;
    }  
  else // (i.e. if 0 <= vgs-vth <= vds)
    {  
      *ids = (k/2)*pow((vgs-Vth),2) + vds/rd;
      *gm  = k*(vgs-Vth);
      *gd  = 1/rd;
      *didT= (dkdT/(2))*pow((vgs-Vth),2);
    }  

  *P       = -(*ids)*vds;
  *dPdT    = -(*didT)*vds;
  *dPdvgs  = -(*gm)*vds;
  *dPdvds  = -((*gd)*vds + (*ids));

}

DEFUN_DLD(Mshichmanhodgesmosfet,args,nargout,
"-*- texinfo -*-\n\
\n\
@deftypefn{Loadable Function} @\n\
{[@var{a},@var{b},@var{c}]=} Mshichmanhodgesmosfet@\n\
(@var{string}, @var{parameters}, @var{parameternames}, @\n\
@var{extvar},@var{intvar},@var{t})\n\
\n\
SBN file implementing Schichman-Hodges MOSFETs model.\n\
\n\
@var{string} is used to select among models. Possible models are:\n\
\n\
@enumerate\n\
@item @var{string} = NMOS (Simplified Shichman-Hodges n-MOSFET)\n\
@item @var{string} = PMOS (Simplified Shichman-Hodges p-MOSFET)\n\
@end enumerate\n\
\n\
Parameters for all the above models are:\n\
@itemize\n\
@item rd     -> parasitic resistance between drain and source\n\
@item W      -> MOSFET width\n\
@item L      -> channel length\n\
@item mu0    -> reference value for mobility\n\
@item Vth    -> threshold voltage\n\
@item Cox    -> oxide capacitance\n\
@item Cgs    -> gate-source capacitance\n\
@item Cgd    -> gate-drain capacitance\n\
@item Cgb    -> gate-bulk capacitance\n\
@item Csb    -> source-bulk capacitance\n\
@item Cdb    -> drain-bulk capacitance\n\
@item Tshift -> shift for reference temperature on MOSFETs\n\
@end itemize\n\
See the @cite{IFF file format specifications} for details about\n\
the output structures.\n\
\n\
@seealso{prs_iff,asm_initialize_system,asm_build_system}\n\
\n\
@end deftypefn\n	  \
")
{
  
  octave_value_list retval; // Contain returned values
  octave_idx_type   nargin = args.length();

  /* Input parameters */
  string        eltype;
  ColumnVector  parameters;
  string_vector parameternames;
  ColumnVector  extvar;
  ColumnVector  intvar(5,0.0);
  double        t;

  /* Model parameters */
  double rd;
  double W,L;
  double mu0,Vth,Cox;
  double Cgs,Cgd,Cgb,Csb,Cdb;
  double Tshift;
  
  /* Model variables */
  double vg; // V-gate
  double vs; // V-source
  double vd; // V-drain
  double vb; // V-bulk
  double T ; // Temperature

  double Qgb;// Charges
  double Qgs;
  double Qgd;
  double Qsb;
  double Qdb;

  double gm; // Conductances et similia
  double gd;
  double ids;
  double didT;
  double P;
  double dPdT;
  double dPdvgs;
  double dPdvds;

  /* Output values */
  Matrix       a(10,10,0.0), b(10,10,0.0);
  ColumnVector c(10,0.0);

  /* Check and retrieve input */
  if (nargin != 6)
    error("Mshichmanhodgesmosfet: wrong number of input parameters.\n");
  else
  {
  /* Retrieve input parameters */
  // Type of MOS-FET
  if (args(0).is_string())
    eltype = args(0).string_value();
  else
    error("Mshichmanhodgesmosfet: argument #1 is expected to be a string.\n");
  // Parameters and parameter names
  if (args(1).length() == args(2).length())
    {
      parameters     = args(1).column_vector_value();
      parameternames = args(2).all_strings();
    }
  else
    error("Mshichmanhodgesmosfet: parameters and parameternames are expected to have the same length.\n");
  // External pins
  if (args(3).length() == 5)
    extvar = args(3).column_vector_value();
  else
    error("Mshichmanhodgesmosfet: five external values expected.\n");
  // Internal variables
  if (args(4).is_empty())
    {}
  else if (args(4).length() == 5)
    intvar = args(4).column_vector_value();
  else
    error("Mshichmanhodgesmosfet: five internal values expected.\n");
  // Time point
  if (args(5).is_real_scalar())
    t = args(5).double_value();
  else
    error("Mshichmanhodgesmosfet: double type value expected as time instant.\n");
  }

  if (!error_state)
    {
      //FIXME: create enum of cases and use switch?
      if (eltype == "NMOS")
	{
	  //FIXME: change parameters to a single map or Octave_map
	  /* Default n-MOS parameters*/
	  rd     = 1e6;
	  W      = 1;
	  L      = 1;
	  mu0    = 1e-5;
	  Vth    = .5;
	  Cox    = 1e-9;
	  Cgb    = Cox;
	  Cgs    = .1*Cox;
	  Cgd    = .1*Cox;
	  Csb    = .1*Cox;
	  Cdb    = .1*Cox;
	  Tshift = 0;
      
	  /* Overwrite parameters */
	  set_parameters(parameters, parameternames, &rd, &W, &L, &mu0, &Vth, &Cox, &Cgs, &Cgd, &Cgb, &Csb, &Cdb, &Tshift);
	  //FIXME: debug
	  //print_parameters(rd, W, L, mu0, Vth, Cox, Cgs, Cgd, Cgb, Csb, Cdb, Tshift);
	  
	  /* Compute model conductance and capacitance */
	  nmos(extvar,mu0,Cox,W,L,Vth,rd,Tshift,&gm,&gd,&ids,&didT,&P,&dPdT,&dPdvgs,&dPdvds);
	  //FIXME: debug
	  //print_values(gm, gd, ids, didT, P, dPdT, dPdvgs, dPdvds);

	  /* Assemble output values*/
	  vg   = extvar(0); // V-gate
	  vs   = extvar(1); // V-source
	  vd   = extvar(2); // V-drain
	  vb   = extvar(3); // V-bulk
	  T    = extvar(4); // Temperature
      
	  Qgb  = intvar(0);
	  Qgs  = intvar(1);
	  Qgd  = intvar(2);
	  Qsb  = intvar(3);
	  Qdb  = intvar(4);
      
	  //FIXME: probably a better way to initialize Matrix exist!
	  /* Dynamic matrix (constant) */
	  a(0,5) =  1;
	  a(0,6) =  1;
	  a(0,7) =  1;
	  a(1,6) = -1;
	  a(1,8) =  1;
	  a(2,7) = -1;
	  a(2,9) =  1;
	  a(3,5) = -1;
	  a(3,8) = -1;
	  a(3,9) = -1;

	  /* Algebraic part (non-linear) */
	  b(1,0) = -gm;
	  b(1,1) = (gm+gd);
	  b(1,2) = -gd;
	  b(1,4) = -didT;
	  b(2,0) = gm;
	  b(2,1) = -(gm+gd);
	  b(2,2) = gd;
	  b(2,4) = didT;
	  b(4,0) = dPdvgs;
	  b(4,1) = -(dPdvgs+dPdvds);
	  b(4,2) = dPdvds;
	  b(4,4) = dPdT;

	  b(5,0) = Cgb;
	  b(5,3) = -Cgb;
	  b(6,0) = Cgs;
	  b(6,1) = -Cgs;
	  b(7,0) = Cgd;
	  b(7,2) = -Cgd;
	  b(8,1) = Csb;
	  b(8,3) = -Csb;
	  b(9,2) = Cdb;
	  b(9,3) = -Cdb;
      
	  b(5,5) = -1;
	  b(6,6) = -1;
	  b(7,7) = -1;
	  b(8,8) = -1;
	  b(9,9) = -1;
      
	  /* Residual */
	  c(0) = 0;
	  c(1) = -ids;
	  c(2) = ids;
	  c(4) = P;
	  c(5) = Cgb*(vg - vb) - Qgb;
	  c(6) = Cgs*(vg - vs) - Qgs;
	  c(7) = Cgd*(vg - vd) - Qgd;
	  c(8) = Csb*(vs - vb) - Qsb;			
	  c(9) = Cdb*(vd - vb) - Qdb; 

	  /* Return values */
	  retval(0) = octave_value(a);
	  retval(1) = octave_value(b);
	  retval(2) = octave_value(c);

	}
      else if (eltype == "PMOS")
	{
	  /* Default p-MOS parameters*/
	  rd     = 1e6;
	  W      = 1;
	  L      = 1;
	  mu0    = 1e-5;
	  Vth    = -.5;
	  Cox    = 1e-9;
	  Cgb    = Cox;
	  Cgs    = .1*Cox;
	  Cgd    = .1*Cox;
	  Csb    = .1*Cox;
	  Cdb    = .1*Cox;
	  Tshift = 0;
      
	  /* Overwrite parameters */
	  set_parameters(parameters, parameternames, &rd, &W, &L, &mu0, &Vth, &Cox, &Cgs, &Cgd, &Cgb, &Csb, &Cdb, &Tshift);

	  /* Compute model conductance and capacitance */
	  pmos(extvar,mu0,Cox,W,L,Vth,rd,Tshift,&gm,&gd,&ids,&didT,&P,&dPdT,&dPdvgs,&dPdvds);
	  
	  /* Assemble output values*/
	  vg   = extvar(0); // V-gate
	  vs   = extvar(1); // V-source
	  vd   = extvar(2); // V-drain
	  vb   = extvar(3); // V-bulk
	  T    = extvar(4); // Temperature
      
	  Qgb  = intvar(0);
	  Qgs  = intvar(1);
	  Qgd  = intvar(2);
	  Qsb  = intvar(3);
	  Qdb  = intvar(4);

	  /* Dynamic matrix (constant) */
	  a(0,5) =  1;
	  a(0,6) =  1;
	  a(0,7) =  1;
	  a(1,6) = -1;
	  a(1,8) =  1;
	  a(2,7) = -1;
	  a(2,9) =  1;
	  a(3,5) = -1;
	  a(3,8) = -1;
	  a(3,9) = -1;
      
	  /* Algebraic part (non-linear) */
	  b(1,0) = -gm;
	  b(1,1) = (gm+gd);
	  b(1,2) = -gd;
	  b(1,4) = -didT;
	  b(2,0) = gm;
	  b(2,1) = -(gm+gd);
	  b(2,2) = gd;
	  b(2,4) = didT;
	  b(4,0) = dPdvgs;
	  b(4,1) = -(dPdvgs+dPdvds);
	  b(4,2) = dPdvds;
	  b(4,4) = dPdT;

	  b(5,0) = Cgb;
	  b(5,3) = -Cgb;
	  b(6,0) = Cgs;
	  b(6,1) = -Cgs;
	  b(7,0) = Cgd;
	  b(7,2) = -Cgd;
	  b(8,1) = Csb;
	  b(8,3) = -Csb;
	  b(9,2) = Cdb;
	  b(9,3) = -Cdb;
      
	  b(5,5) = -1;
	  b(6,6) = -1;
	  b(7,7) = -1;
	  b(8,8) = -1;
	  b(9,9) = -1;

	  /* Residual */
	  c(0) = 0;
	  c(1) = -ids;
	  c(2) = ids;
	  c(4) = P;
	  c(5) = Cgb*(vg - vb) - Qgb;
	  c(6) = Cgs*(vg - vs) - Qgs;
	  c(7) = Cgd*(vg - vd) - Qgd;
	  c(8) = Csb*(vs - vb) - Qsb;			
	  c(9) = Cdb*(vd - vb) - Qdb;

	  /* Return values */
	  retval(0) = octave_value(a);
	  retval(1) = octave_value(b);
	  retval(2) = octave_value(c);
	}
      else
	error("Mshichmanhodgesmosfet: unknown element type.\n");

    }

  return retval;
}
