#! /bin/bash


#N=10
#num_of_Is="40 80 120 160 200 240 280 320"
#Bs_ratio="0 0.1 0.4 0.8"
#GRN_name="mammalian fission budding arabidopsis thelper tcrNet"


# N=1
# num_of_Is="440 480 520"
# GRN_name="mammalian"


# for E in $num_of_Is
# do
#   python3.6 gen-E-compare-multiB.py mammalian $E 0.1 compare 4 multiB.log
# done

num_of_Is="200 240 280 320 360 400 440 480"
for E in $num_of_Is
do
  for (( K=0; K<10; K++ ))
  do 
    #echo $(date +%R) 
    python3.6 gen-E-compare.py budding $E 0.1 compare 5 fact.log
  done
done

# for (( K=0; K<10; K++ ))
# do 
#   echo $(date +%R) 
#   clingo MHSC_v0.lp auxiliary/mammalian.lp examples_asp/mammalian_v0_$K.lp --outf=0 --verbose=1  --time-limit=3600 --models 0 --opt-mode=opt --quiet=2 
# done


# N=1
# num_of_Is="" 
# Bs_ratio="0.1"
# GRN_name="arabidopsis thelper tcrNet"

# for (( K=1; K<=N; K++ ))
#   do
#   for G in $GRN_name
#   do

#     if [[ $G == "budding" || $G == "arabidopsis"  || $G == "thelper"  ]]
#     then
#       num_of_Is="400 440 480"
#     elif [ $G == "tcrNet"  ]
#     then
#       num_of_Is="360 400 440 480"
#     else #default for mammalian and fission
#       num_of_Is="40 80 120 160 200 240 280 320 360 400 440 480"
#     fi

#     for E in $num_of_Is
#     do
#       python3.6 gen-E-compare.py $G $E $Bs_ratio compare 4 margin.log
#     done
#   done
# done
