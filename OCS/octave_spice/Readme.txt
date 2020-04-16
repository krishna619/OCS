Extract the Tar file
How to use the function:
type "make" to simulate all cir files.
execute spice_readfile_test.m to view the simulation results

Note: the function does not handle complex plot structures as indicated 
by the ngspice sourcefile rawspice.c.
not supported functions are:
	* multidimension plots
        * plotting commands
	* grid settings, min/max, ..

after that
user@host:~$ ngspice -b rc2.cir

Now we can import the file ngspice produced ( named rawspice.raw) into octave with the following command:


user@host:~$ octave
GNU Octave, version 2.1.69 (i386-pc-linux-gnu).
Copyright (C) 2005 John W. Eaton.
This is free software; see the source code for copying conditions.
There is ABSOLUTELY NO WARRANTY; not even for MERCHANTIBILITY or
FITNESS FOR A PARTICULAR PURPOSE.  For details, type `warranty'.

Additional information about Octave is available at http://www.octave.org.

Please contribute if you find this software useful.
For more information, visit http://www.octave.org/help-wanted.html

Report bugs to <bug@octave.org> (but first, please read
http://www.octave.org/bugs.html to learn how to write a helpful report).

octave:1> k=spice_readfile("rawspice.raw");

The plot contained in file "rawspice.raw" is put into the k vector and the following output is generated:


s =
{
  commands = {}
date = Mon Sep 12 14:32:22  2005
  dimensions = 0
  no_points = 0
  no_variables = 0
  options = {}
plotname = Transient Analysis
  title = * simulation de RC2
}

octave:2> 

Now we can plot the results (as in nutmeg, but we are using Octave):

octane:2>plot (k)
