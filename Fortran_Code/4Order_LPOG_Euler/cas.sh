#!/bin/bash
for caso in 10 11 ; do
    PP=1
for ptsxy in "'33'" "'65'" "'129'" "'257'" "'513'" "'1025'" ; do
    WW=1
for wi in "'0.05d0'" ; do
    RR=1
for re in "'0.28d0'" ; do
    GG=1
for alphag in "'0.0d0'" ; do

    cat cas.sh.template | sed -e "s/__CC__/$caso/" | sed -e "s/__PPPP__/$ptsxy/" | sed -e "s/__RRRR__/$re/" | sed -e "s/__AAAA__/$alphag/" | sed -e "s/__WWWW__/$wi/" > cas_0$(($caso))_PTS$(($PP))_Wi$(($WW))_Re$(($RR))_AlphaG$(($GG)).sh

    cat cas.job.template | sed -e "s/__CC__/$caso/" | sed -e "s/__PP__/$PP/" | sed -e "s/__RR__/$RR/" | sed -e "s/__GG__/$GG/" | sed -e "s/__WW__/$WW/" > cas_0$(($caso))_PTS$(($PP))_Wi$(($WW))_Re$(($RR))_AlphaG$(($GG)).job

    mkdir -p 'Case_0'$caso'_PTS'$PP'_Wi'$WW'_Re'$RR'_AlphaG'$GG || true 
    cd 'Case_0'$caso'_PTS'$PP'_Wi'$WW'_Re'$RR'_AlphaG'$GG
    cp ../*.f . 
    cp ../makefile . 
    cp ../par.nn* . 
    cp ../comm.var* .
    mv ../cas_0$(($caso))_PTS$(($PP))_Wi$(($WW))_Re$(($RR))_AlphaG$(($GG)).job .
    mv ../cas_0$(($caso))_PTS$(($PP))_Wi$(($WW))_Re$(($RR))_AlphaG$(($GG)).sh .

    qsub -P project1 cas_0$(($caso))_PTS$(($PP))_Wi$(($WW))_Re$(($RR))_AlphaG$(($GG)).job

    qstat -u juniormr

    totaljob
    
    sleep 1
    cd ..
    GG=$(($GG+1))
done
    RR=$(($RR+1))
done
    WW=$(($WW+1))
done
    PP=$(($PP+1))
done
done
