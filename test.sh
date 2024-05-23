#!/bin/zsh

# ./test kansoku N K R
# ./test $KANSOKU 9 3 50 >> $KANSOKU/log/xm.txt 

rm grav/log/xm.txt

for (( i=1; i < 1000; ++i))
do
  echo $((i))
  ./test g l 100 6 10 > grav/log/xm.txt
done


