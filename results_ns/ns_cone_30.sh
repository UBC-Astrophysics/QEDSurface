ctioga2 --name ns_cone_30 \
	--xfact 4.14e-18 \
	--legend '$r_c=700$m' \
	cone_summary_4.1@'(if $5==30 && $4==1e6 then $2 else 0.0/0.0 end):-$6/$7' \
	--legend '$r_c=1000$m' \
	cone_summary_6.1@'(if $5==30 && $4==1e6 then $2 else 0.0/0.0 end):-$6/$7' \
	--legend '$r_c=1800$m' \
	cone_summary_10.1@'(if $5==30 && $4==1e6 then $2 else 0.0/0.0 end):-$6/$7' \
	-x 'Energy [keV]' -y 'QED Decrement at 30 degrees'
