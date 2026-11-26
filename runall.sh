#!/bin/sh
# runall.sh -- reruns all the programs of the paper, writing one .log file for each, in about 15 minutes, most of it
# for autgroup.py, whose output aut.npy (about 180 MB) is used by groundstates.py, and is then removed, and for
# dynamics.py.  The enumerations of order 28 (gaethje.py, and transposes.py on matrices28, about 45 minutes each on
# 2 cores), 32 (golovkin.py and valentina.py, about 1 hour each on a laptop) and 36 (frazier.py, about 1 hour on a
# laptop), and the large simulations (stipe.py, in stipe/, about 1 hour on a laptop) are not rerun here, their outputs
# being in gaethje_output, matrices28T_raw, golovkin_output, valentina_output, frazier_output, and stipe/stipe_output,
# and neither are the permanents of order 28 (permanent.py on matrices28, about 25 minutes on 2 cores), which are in
# matrices28_perm.json.
set -x
python3 glowlaw.py > glowlaw.log 2>&1
python3 twoadic.py > twoadic.log 2>&1
python3 moments.py > moments.log 2>&1
python3 universal.py > universal.log 2>&1
python3 spectrum16.py > spectrum16.log 2>&1
python3 thermo.py > thermo.log 2>&1
python3 figures.py > figures.log 2>&1
python3 crystal.py > crystal.log 2>&1
python3 hypercube.py > hypercube.log 2>&1
cc -O3 -o dos dos.c
python3 hypercubedos.py > hypercubedos.log 2>&1
python3 make32.py > make32.log 2>&1
python3 autgroup.py matrices/had16_3.txt aut.npy > autgroup.log 2>&1
python3 groundstates.py > groundstates.log 2>&1
rm -f aut.npy
python3 multisets.py > multisets.log 2>&1
python3 mckay.py mckay/had24.txt 24 matrices24 > mckay.log 2>&1
python3 mckay.py mckay/had28.txt 28 matrices28 >> mckay.log 2>&1
cc -O3 -o dos dos.c
python3 dosall.py matrices24 dos24.json > dos24.log 2>&1
python3 structure.py dos24.json matrices24 > structure24.log 2>&1
python3 transposes.py matrices24 > transposes.log 2>&1
python3 transpose.py matrices24 dos24.json > transpose24.log 2>&1
cc -O3 -o perm perm.c
python3 permanent.py matrices24 dos24.json > permanent24.log 2>&1
python3 dosall.py matrices28 dos28.json gaethje_output > dos28.log 2>&1
python3 structure.py dos28.json matrices28 gaethje_output > structure28.log 2>&1
python3 transpose.py matrices28 dos28.json gaethje_output > transpose28.log 2>&1
python3 dosall.py matrices32 dos32.json golovkin_output valentina_output > dos32.log 2>&1
python3 structure.py dos32.json matrices32 > structure32.log 2>&1
python3 make36.py > make36.log 2>&1
python3 dosall.py matrices36 dos36.json frazier_output > dos36.log 2>&1
python3 structure.py dos36.json matrices36 frazier_output > structure36.log 2>&1
python3 crystal.py > crystal.log 2>&1
python3 residues.py > residues.log 2>&1
python3 pictures.py > pictures.log 2>&1
python3 dynlarge.py > dynlarge.log 2>&1
cc -O3 -o mc mc.c -lm
python3 dynamics.py > dynamics.log 2>&1
