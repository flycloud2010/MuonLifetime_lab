for i in {4,5,6,7,8,9}; do

echo "Analysing n_tau: $i"

root -b -l << EOF
.L exp_baseline.C
ExpBase(5000,true,true,4,$i)
.q

EOF

done