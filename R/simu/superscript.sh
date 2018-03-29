rm -rf ~/simuBactDating
mkdir ~/simuBactDating
nohup seq 1 300 | parallel -j 20 ~/BactDating/R/simu/simu.sh {} &
