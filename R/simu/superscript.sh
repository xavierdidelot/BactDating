rm -rf ~/simuCreDating
mkdir ~/simuCreDating
nohup seq 1 200 | parallel -j 30 ~/CreDating/R/simu/simu.sh {} &
