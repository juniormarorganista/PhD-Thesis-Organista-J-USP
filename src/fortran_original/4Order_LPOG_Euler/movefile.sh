#!/bin/bash
rm -rf Results_2
sleep 1
mkdir -p Results_2 || true 
sleep 1
cd Results_2/

PP=1
tt=5
#for ptsxy in "'17'" "'33'" "'49'" "'65'" "'81'" "'97'" "'113'" "'129'" "'145'" "'161'" "'177'" "'193'" "'209'" "'225'" "'241'" "'257'" ; do
for ptsxy in "'17'" "'33'" "'49'" "'65'" "'81'" "'97'" "'113'" "'129'" ; do
   BB=1
   for beta in "'0.1d0'" "'0.5d0'" "'0.9d0'" ; do
      #for caso in 1 2 3 4 5 6 7 8 9 10 11 12 ; do
      for caso in 12 ; do
      WW=1
      for wi in "'1.0d0'" ; do
         RR=1
         for re in "'100.0d0'" "'400.0d0'" "'1000.0d0'" ; do
            GG=1
            for alphag in "'0.0d0'" "'0.5d0'" ; do
               XX=1
               for xi in "'0.0d0'" ; do
                
                   cp ../'Case_0'$caso'_PTS'$PP'_Wi'$WW'_Re'$RR'_AlphaG'$GG'_Xi'$XX/*.dat .
                   cp ../'Case_0'$caso'_PTS'$PP'_Wi'$WW'_Re'$RR'_AlphaG'$GG'_Xi'$XX/*.vtk .

                   echo $ptsxy:'Case_0'$caso'_PTS'$PP'_Wi'$WW'_Re'$RR'_AlphaG'$GG'_Xi'$XX/*.dat
               
                   #sleep $(($tt))
                   sleep 0.1
               XX=$(($XX+1))
               done
            GG=$(($GG+1))
            done
         RR=$(($RR+1))
         done
      WW=$(($WW+1))
      done
   done
   BB=$(($BB+1))
   done
PP=$(($PP+1))
tt=$(($tt*$PP))
done


