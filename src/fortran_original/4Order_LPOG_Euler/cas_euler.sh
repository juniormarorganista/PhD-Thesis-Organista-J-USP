#!/bin/bash
PP=41
tt=5
for ptsxy in "'33'" "'65'" "'129'" "'257'" "'513'" ; do
   BB=1
   for beta in "'0.1d0' '0.5d0' '1.0d0'" ; do
      for caso in 3 ; do
      WW=1
      for wi in "'2.0d0'" ; do
         RR=1
         for re in "'10.0d0'" "'100.0d0'" "'1000.0d0'" ; do
            GG=1
            for alphag in "'0.0d0'" ; do
               XX=1
               for xi in "'0.5d0'" ; do
               
                   cat cas.sh.template | sed -e "s/__PPPP__/$ptsxy/" | sed -e "s/__CC__/$caso/" | sed -e "s/__BBBB__/$beta/" | sed -e "s/__WWWW__/$wi/" | sed -e "s/__RRRR__/$re/" | sed -e "s/__AAAA__/$alphag/" | sed -e "s/__XXXX__/$xi/" > cas_0$(($caso))_PTS$(($PP))_Wi$(($WW))_Re$(($RR))_AlphaG$(($GG))_Xi$(($XX)).sh
               
                   cat cas.job.template | sed -e "s/__CC__/$caso/" | sed -e "s/__PP__/$PP/" | sed -e "s/__RR__/$RR/" | sed -e "s/__GG__/$GG/" | sed -e "s/__WW__/$WW/" | sed -e "s/__XX__/$XX/"  > cas_0$(($caso))_PTS$(($PP))_Wi$(($WW))_Re$(($RR))_AlphaG$(($GG))_Xi$(($XX)).job
               
                   mkdir -p 'Case_0'$caso'_PTS'$PP'_Wi'$WW'_Re'$RR'_AlphaG'$GG'_Xi'$XX || true 
                   cd 'Case_0'$caso'_PTS'$PP'_Wi'$WW'_Re'$RR'_AlphaG'$GG'_Xi'$XX
               
                   cp ../*.f . 
                   cp ../makefile . 
                   cp ../par.nn* . 
                   cp ../comm.var* .
                   mv ../cas_0$(($caso))_PTS$(($PP))_Wi$(($WW))_Re$(($RR))_AlphaG$(($GG))_Xi$(($XX)).job .
                   mv ../cas_0$(($caso))_PTS$(($PP))_Wi$(($WW))_Re$(($RR))_AlphaG$(($GG))_Xi$(($XX)).sh .
               
                   qsub -P project1 cas_0$(($caso))_PTS$(($PP))_Wi$(($WW))_Re$(($RR))_AlphaG$(($GG))_Xi$(($XX)).job
                   qstat -u juniormr
                   totaljob
                   
                   #sleep $(($tt))
                   #sleep 1
                   cd ..
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
