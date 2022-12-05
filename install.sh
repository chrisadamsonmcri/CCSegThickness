#!/bin/bash


if [ $UID != 0 ] ;
then
printf "\nNeeds to be root | sudo user to install dependencies\nTry running as \n\tsudo bash install.sh\n " ;
exit
fi

printf "\n Updating system and Installing required packages ......\n\n\n  "

sudo apt install python python-dev python-pip python-numpy python-scipy python-nibabel python-opencv  python-matplotlib libsuitesparse-dev cython3 -y

sudo -H pip3 install -U scikits.sparse
