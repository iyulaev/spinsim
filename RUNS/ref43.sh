./spinsim 0.000000e+00 0.000000e+00 1.000000e-02 6.000000e+05 6.000000e+05 6.000000e+02 0.000000e+00 4.000000e+05 6.000000e+04 1.000000e+07 0.000000e+00 0.000000e+00 6.000000e+00 1.200000e+01 0.000000e+00 26 -5.000000e+11 2.000000e+11 200 0.000000e+00 5.000000e+03 200 2.000000e-08 1 1 3.500000e+01 3.500000e+01 0.000000e+00 0.000000e+00 1.000000e-02 9.900000e-01 0.000000e+00 1.000000e-02 9.900000e-01 0.000000e+00 0.000000e+00 1.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 3.000000e+00 3.000000e+00 3.000000e+00 0.000000e+00 5.000000e-02 3.500000e-01 3.500000e-01 3.000000e-01 0 4 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 1.000000e-02 1.000000e-02 1.000000e-02 0.000000e+00

mkdir ref43_norm
mv mag* ref43_norm/.

mv spin_helpers.c spin_helpers.norm.c
mv spin_helpers.tunnel.c spin_helpers.c
make

./spinsim 0.000000e+00 0.000000e+00 1.000000e-02 6.000000e+05 6.000000e+05 6.000000e+02 0.000000e+00 4.000000e+05 6.000000e+04 1.000000e+07 0.000000e+00 0.000000e+00 6.000000e+00 1.200000e+01 0.000000e+00 26 -5.000000e+11 2.000000e+11 200 0.000000e+00 5.000000e+03 200 2.000000e-08 1 1 3.500000e+01 3.500000e+01 0.000000e+00 0.000000e+00 1.000000e-02 9.900000e-01 0.000000e+00 1.000000e-02 9.900000e-01 0.000000e+00 0.000000e+00 1.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 3.000000e+00 3.000000e+00 3.000000e+00 0.000000e+00 5.000000e-02 3.500000e-01 3.500000e-01 3.000000e-01 0 4 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 1.000000e-02 1.000000e-02 1.000000e-02 0.000000e+00

mkdir ref43_tunn
mv mag* ref43_tunn/.

mv spin_helpers.c spin_helpers.tunnel.c
mv spin_helpers.norm.c spin_helpers.c
make
