###Run tleap and minimization with amber###
tleap -s -f tleap_imp.in > out.leap
parmed pd.parm7 parmed.in
pmemd.cuda -O -i min.in  -p pd_hmass.parm7 -o min.out -c pd.rst7 -r min_imp.rst7  -ref pd.rst7 -x min.crd -inf min.info
###Output processing###
ambpdb -p pd_hmass.parm7 -c min_imp.rst7 -bres > out-min.pdb
