ctioga2 --name ns_radius \
	cone_summary_6.1@'(if $5>0 && $2==0.5e18 && $4==1e6 then $5 else 0.0/0.0 end):-$6/$7' \
	@'(if $5>0 && $2==1e18 && $4==1e6 then $5 else 0.0/0.0 end):-$6/$7' \
	@'(if $5>0 && $2==2e18 && $4==1e6 then $5 else 0.0/0.0 end):-$6/$7' 
