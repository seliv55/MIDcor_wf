# standing on the directory where the package is (midcor)

library(devtools)
build()
install()
******************-------------------------******************
# run script:
vs@vit:~/phenomen/docker-midcor/scripts$ ./runMidcor.R -i ~/phenomen/outin.csv -o ~/phenomen/tes_out.csv

******************-------------------------******************

# crear Docker container:

# go to the directory where the dockerfile is
cd /home/vs/phenomen/docker-midcor
# create container from dockerfile
sudo docker build -t midcor:0.1 .

# check existing images:
sudo docker images

# try to run the tool from docker
sudo docker run --name=midcor-d -i -t midcor:0.1 -i outin.csv -o my_out.csv 

sudo docker run -i -t -v $PWD:/data midcor:0.1 -i /data/outin.csv -o /data/my_out.csv
