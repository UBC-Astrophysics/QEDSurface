(NR==1) { mu=$3; nu=$4; mass=$5; radius=$6; alpha=$7; qsum=0; nq=0; }
($1!="#" && $10<4.1) { qsum+=$15; nq++; }
END { print mu, nu, mass, radius, alpha, qsum, nq; }
