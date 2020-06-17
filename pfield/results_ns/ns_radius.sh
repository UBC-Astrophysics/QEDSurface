for ang in 15 30 45 60 75
do
    ctioga2 --name ns_radius_${ang} \
	    --xfact 4.14e-18 \
	    --legend '10 km' \
	    cone_summary_6.1@"(if \$5==${ang} && \$4==1e6 then \$2 else 0.0/0.0 end):-\$6/\$7" \
	    --legend '12 km' \
	    @"(if \$5==${ang} && \$4==1.2e6 then \$2 else 0.0/0.0 end):-\$6/\$7" \
	    --legend '14 km' \
	    @"(if \$5==30 && \$4==1.4e6 then \$2 else 0.0/0.0 end):-\$6/\$7" \
	    -x 'Energy [keV]' -y "Polarized Fraction at ${ang} degrees"
done

