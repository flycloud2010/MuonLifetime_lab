
for i in {4,5,6,7,8,9}; do

echo " Number of weeks $1"
echo "----- Analysing n_tau: $i -----"

root -b -l << EOF
.L double_exp.C
Fit(1000,true,false,$1,$i)
.q

EOF

done
