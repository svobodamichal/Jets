#!/bin/sh
export PARTICLE="Pi"
export FILELIST="filelists/Pi_file.list"
export CENTRAL=0
export GLOBAL=0
export NFIT=15

starver SL11d_embed
pwd
#ls -latr
root4star -l <<EOF
.L StMiniMcTree.C
StMiniMcTree uu
uu.Loop()
EOF
#.q
