rm -rf ~/simuBactDating
mkdir ~/simuBactDating
nohup seq 1 200 | parallel -j 30 ~/BactDating/R/simu/simu.sh {} &
