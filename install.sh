#!/bin/bash


if [ $UID != 0 ] ;
then
printf "\nNeeds to be root | sudo user to install dependencies\nTry running as \n\tsudo bash install.sh\n " ;
exit
fi

printf "\n Updating system and Installing required packages ......\n\n\n  "


sudo apt update
sudo apt upgrade -y
sudo apt install python python-dev python-pip python-numpy python-scipy python-nibabel python-opencv  python-matplotlib -y
sudo apt install libsuitesparse-dev -y

sudo -H pip install -U scikits.sparse
