for ang in 15 30 45 60 75
do
	   awk "(\$5==${ang} && \$4<1.15e6) { ii=\$2/0.25e18; ss[ii]+=\$6; nn[ii]+=\$7; } END { for (ii=0;ii<11;ii++) { print 1e30, ii*0.25e18, 2.06e5, 1.15e6, ${ang}, ss[ii], nn[ii];}}" cone_summary_6.1 > cone_6.1_${ang}_1000m
	   awk "(\$5==${ang} && \$4<1.25e6) { ii=\$2/0.25e18; ss[ii]+=\$6; nn[ii]+=\$7; } END { for (ii=0;ii<11;ii++) { print 1e30, ii*0.25e18, 2.06e5, 1.25e6, ${ang}, ss[ii], nn[ii];}}" cone_summary_6.1 > cone_6.1_${ang}_2000m
	   awk "(\$5==${ang} && \$4<1.35e6) { ii=\$2/0.25e18; ss[ii]+=\$6; nn[ii]+=\$7; } END { for (ii=0;ii<11;ii++) { print 1e30, ii*0.25e18, 2.06e5, 1.35e6, ${ang}, ss[ii], nn[ii];}}" cone_summary_6.1 > cone_6.1_${ang}_3000m
	   awk "(\$5==${ang} && \$4<1.65e6) { ii=\$2/0.25e18; ss[ii]+=\$6; nn[ii]+=\$7; } END { for (ii=0;ii<11;ii++) { print 1e30, ii*0.25e18, 2.06e5, 1.65e6, ${ang}, ss[ii], nn[ii];}}" cone_summary_6.1 > cone_6.1_${ang}_6000m
ctioga2 --name ns_height_${ang} \
	--xfact 4.14e-18 \
	--legend '$h_c=0$m' \
	cone_summary_6.1@"(if \$5==${ang} && \$4==1e6 then \$2 else 0.0/0.0 end):-\$6/\$7" \
	--legend '$h_c=1000$m' \
	cone_6.1_${ang}_1000m@'$2:-$6/$7' \
	--legend '$h_c=2000$m' \
	cone_6.1_${ang}_2000m@'$2:-$6/$7' \
	--legend '$h_c=3000$m' \
	cone_6.1_${ang}_3000m@'$2:-$6/$7' \
	--legend '$h_c=6000$m' \
	cone_6.1_${ang}_6000m@'$2:-$6/$7' \
	-x 'Energy [keV]' -y "Polarized Fraction at ${ang} degrees"
done


