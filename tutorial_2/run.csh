#!/bin/csh -f
rm -rf DELETE

rm ref.gro
echo 0 | gmx trjconv -f start.pdb -s start.pdb -o ref.gro
sed -i 1,2d ref.gro

rm ref.pdb
echo 0 | gmx trjconv -f start.pdb -s start.pdb -o ref.pdb

rm xyz.dat
python xyz.py >> xyz.dat

rm Romega.dat
rm a.out
g++ romega.cpp -lm
./a.out

rm a.out
rm neighbourhood.dat
rm for_deri.dat
rm population_neigh.dat
rm PMATRIX.dat
gcc pmat.c -Ddsyev=dsyev_ -lm -llapack -lblas
./a.out

rm COLVAR
mv neighbourhood.dat  neighbourhood1.dat
mv for_deri.dat  for_deri1.dat
mv population_neigh.dat  population_neigh1.dat
mv PMATRIX.dat  PMATRIX1.dat
mv ref.gro ref1.gro
time ../install_plumed_with_gnap/plumed/utility1/bin/plumed  driver --mf_xtc 1000frames.xtc --plumed plumed.dat



mkdir DELETE
 mv ref1.gro DELETE/ref1.gro
 mv PMATRIX1.dat DELETE/PMATRIX1.dat
 mv population_neigh1.dat  DELETE/population_neigh1.dat
 mv for_deri1.dat  DELETE/for_deri1.dat
 mv neighbourhood1.dat  DELETE/neighbourhood1.dat
 mv ref.pdb DELETE/ref.pdb
 mv Romega.dat DELETE/Romega.dat
 mv xyz.dat  DELETE/xyz.dat

