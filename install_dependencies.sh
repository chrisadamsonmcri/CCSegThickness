#!/bin/bash

if [ $UID != 0 ] ;
then
printf "\nNeeds to be root | sudo user to install dependencies\nTry running as \n\tsudo bash install.sh\n " ;
exit
fi

printf "\n Updating system and Installing required packages ......\n\n\n  "

sudo apt install python3-dev python3-pip python3-numpy python3-scipy python3-nibabel python3-opencv python3-matplotlib libsuitesparse-dev -y

sudo -H pip3 install -U scikits.sparse
