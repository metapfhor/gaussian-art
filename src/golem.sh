#!/bin/bash
scp -P 2222 *.f90 mmotif.concordia.ca:~/art/src
ssh -p 2222 mmotif.concordia.ca 'cd ~/art/src; make clean; make all'
scp -P 2222 mmotif.concordia.ca:~/art/arttest.exe laurent@golem.concordia.ca:/home/laurent/art
