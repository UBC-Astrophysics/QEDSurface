ctioga2 --name ns_nu_60 \
	--xfact 4.14e-18 \
	--legend '10 km' \
	cone_summary_6.1@'(if $5==60 && $4==1e6 then $2 else 0.0/0.0 end):-$6/$7' \
	--legend '12 km' \
	@'(if $5==60 && $4==1.2e6 then $2 else 0.0/0.0 end):-$6/$7' \
	--legend '14 km' \
	@'(if $5==60 && $4==1.4e6 then $2 else 0.0/0.0 end):-$6/$7' \
	-x 'Energy [keV]' -y 'QED Decrement at 60 degrees'
