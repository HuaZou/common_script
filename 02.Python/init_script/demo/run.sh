# -m python:module R:packages perl:module
python ../init.py -t python -o test -m sys,os,re
python ../init.py -t R -o test -m dplyr,ggplot2
python ../init.py -t perl -o test -m test
