#!/bin/zsh

# select grad or pote
KANSOKU=$1

# K
K=$2

# param
SCATTER_MAX=$3

# move mass on lattice
./onlattice $((K)) > $KANSOKU/log/onlattice.txt < $KANSOKU/dp/log/xm.txt

# initialize scattering
# rm $KANSOKU/$/scatter/*
rm $KANSOKU/PMS/*
cp $KANSOKU/log/onlattice.txt $KANSOKU/PMS/scatter0.txt
cp $KANSOKU/log/onlattice.txt $KANSOKU/PMS/old.txt

# scattering
for ((k = 1; k <= $SCATTER_MAX; ++k)); do
  if [ $((k % 10)) -eq 0 ] 
  then
    ./scatter $((k)) < $KANSOKU/PMS/old.txt | tee $KANSOKU/PMS/new.txt > $KANSOKU/PMS/scatter$((k)).txt
  else 
    ./scatter $((k)) < $KANSOKU/PMS/old.txt 1> $KANSOKU/PMS/new.txt 
    # ./scatter $((k)) < $KANSOKU/scatter/old.txt | tee $KANSOKU/scatter/new.txt > $KANSOKU/scatter/scatter$((k)).txt
  fi
  # ./scatter $((k)) < $KANSOKU/PMS/old.txt | tee $KANSOKU/PMS/new.txt > $KANSOKU/PMS/scatter$((k)).txt

  if [ $? -eq 1 ]
  then
    echo end scattering
    exit 1
  else
    cp $KANSOKU/PMS/new.txt  $KANSOKU/PMS/old.txt
  fi

done

echo need more scattering


