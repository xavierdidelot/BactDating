rm -rf ~/simuCreDating
mkdir ~/simuCreDating
nohup seq 1 100 | parallel -j 30 ~/CreDating/R/simu/simu.sh {} &
