#!/bin/bash
PP=1
tt=5
for ptsxy in "'17'" "'33'" "'49'" "'65'" "'81'" "'97'" "'113'" "'129'" "'145'" "'161'" "'177'" "'193'" "'209'" "'225'" "'241'" "'257'" ; do
   BB=1
   for beta in "'0.d0'" "'0.1d0'" "'0.5d0'" "'0.9d0'" "'1.0d0-1.d-10'" ; do
      CC=1
      for caso in 1 2 3 4 5 6 7 8 9 10 11 12 ; do
      TS=1
      for tsim in 1 ; do
      RR=1
      for re in "'1.0d0'" "'10.0d0'" "'100.0d0'" "'400.0d0'" "'1000.0d0'" "'10000.0d0'" ; do
         WW=1
         for wi in "'1.0d0'" "'5.0d0'" "'100.0d0'" ; do
            GG=1
            for alphag in "'0.0d0'" "'0.5d0'" "'0.85d0'" ; do
               XX=1
               for xi in "'0.0d0'" "'0.1d0'" "'0.5d0'" ; do
               EP=1
               for epsilon in "'0.0d0'" "'0.5d0'" ; do
               AT=1
               for at in "'-0.5d0'" "'0.0d0'" "'0.5d0'" ; do
               DT=1
               for dt in "'1.0d-05'" ; do
               
                   cat cas.sh.template | sed -e "s/__TSTS__/$tsim/" | sed -e "s/__DTDT__/$dt/" | sed -e "s/__EPEP__/$epsilon/" | sed -e "s/__ATAT__/$at/" | sed -e "s/__PPPP__/$ptsxy/" | sed -e "s/__CC__/$caso/" | sed -e "s/__BBBB__/$beta/" | sed -e "s/__WWWW__/$wi/" | sed -e "s/__RRRR__/$re/" | sed -e "s/__AAAA__/$alphag/" | sed -e "s/__XXXX__/$xi/" > cas_0$(($caso))_PTS$(($PP))_Wi$(($WW))_Re$(($RR))_AlphaG$(($GG))_Xi$(($XX)).sh
               
                   cat cas.job.template | sed -e "s/__CC__/$caso/" | sed -e "s/__PP__/$PP/" | sed -e "s/__RR__/$RR/" | sed -e "s/__GG__/$GG/" | sed -e "s/__WW__/$WW/" > cas_0$(($caso))_PTS$(($PP))_Wi$(($WW))_Re$(($RR))_AlphaG$(($GG)).job

                   mkdir -p 'Case_0'$caso'_PTS'$PP'_Wi'$WW'_Re'$RR'_AlphaG'$GG'_Xi'$XX || true 
                   cd 'Case_0'$caso'_PTS'$PP'_Wi'$WW'_Re'$RR'_AlphaG'$GG'_Xi'$XX
               
                   cp ../*.f . 
                   cp ../makefile . 
                   cp ../par.nn* . 
                   cp ../comm.var* .
                   mv ../cas_0$(($caso))_PTS$(($PP))_Wi$(($WW))_Re$(($RR))_AlphaG$(($GG)).job .
                   mv ../cas_0$(($caso))_PTS$(($PP))_Wi$(($WW))_Re$(($RR))_AlphaG$(($GG))_Xi$(($XX)).sh .
                   qsub -P project1 cas_0$(($caso))_PTS$(($PP))_Wi$(($WW))_Re$(($RR))_AlphaG$(($GG)).job
                   qstat -u juniormr
                   totaljob
                   pwd
                   sleep 10
                   cd ..

               DT=$(($DT+1))
               done
               AT=$(($AT+1))
               done
               EP=$(($EP+1))
               done
               XX=$(($XX+1))
               done
            GG=$(($GG+1))
            done
         WW=$(($WW+1))
         done
      RR=$(($RR+1))
      done
      TS=$(($TS+1))
      done
   CC=$(($CC+1))
   done
   BB=$(($BB+1))
   done
PP=$(($PP+1))
tt=$(($tt*$PP))
done
