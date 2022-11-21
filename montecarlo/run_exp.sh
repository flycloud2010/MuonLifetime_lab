for i in {1,2,3}; do

echo "Analysing week: $i"

root -b -l << EOF
.L exponentialMC.C
MC_CDF(10000, true, true, $i)
.q

EOF

done