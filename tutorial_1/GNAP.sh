
var2=(`ls COLVAR`)
if [ $var2 == 'COLVAR' ]; then rm COLVAR; fi



../install_plumed_with_gnap/plumed/utility1/bin/plumed  driver --mf_xtc 1000frames.xtc --plumed plumed.dat
