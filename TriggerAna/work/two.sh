#!/bin/bash

root -b -q 'Fun4SimDst.C(0, "list_e906.txt" ,"uDST_e906.root")'
root -b -q 'Fun4SimDst.C(0, "list_e1039.txt","uDST_e1039.root")'

root -b -q 'Fun4SimMicroDst.C("uDST_e906.root","sim_e906.root")'
root -b -q 'Fun4SimMicroDst.C("uDST_e1039.root","sim_e1039.root")'
