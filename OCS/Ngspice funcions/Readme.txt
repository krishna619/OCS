//How to use //

ngspice_init();
ngspice_command(["source simu.cir"; "tran 0.1us 6us uic"]);
[vn1, vn2] = ngspice_get_vector(["v(n1)"; "v(n2)"]);