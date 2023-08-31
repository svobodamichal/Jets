#!/bin/bash

PLOTEFFI=1
PLOTMOMRES=0

for CUTSET in "cutset1" #"nfit12" "nfit15" "nfit20"
do
if [ $PLOTEFFI -eq 1 ]; then
for SYSTEM in "pp" "AuAu"
do
for CENTRAL in 0 1
do
root -l -b <<EOF
.L momentum_resolution_and_efficiency.C
efficiency($CENTRAL,"$CUTSET","$SYSTEM")
.q
EOF
done
done
fi

if [ $PLOTMOMRES -eq 1 ]; then
for PARTICLE in "Pi0" "Pi" "K" "P" "jet"
do
if [ $PARTICLE == "Pi0" -o $PARTICLE == "jet" ]; then
CENTRALS="2"
else 
CENTRALS="0 1"
fi
for CENTRAL in $CENTRALS
do
root -l -b <<EOF
.L momentum_resolution_and_efficiency.C
momentum_resolution("$PARTICLE",$CENTRAL, "$CUTSET")
.q
EOF
done
done
fi

done #cutset
