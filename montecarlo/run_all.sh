for i in {4,4.5,5,5.5,6,6.5,7,7.5,8,8.5,9}; do

echo "Analysing n_tau: $i"

root -b -l << EOF
.L exp_baseline.C
FitBaseline(5000,true,true,1,$i)
.q

EOF

done