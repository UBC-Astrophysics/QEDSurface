BEGIN {
    print("#   Radius[cm]     Phi[deg]   Energy[keV]     I          Q/I        mu[Gcm3]      Mass[cm]        Q      Phi[rad]")
}
($1!~"#") {
    printf("%s %8g\n",$0,$2*atan2(1.0,1.0)/45.0);
}
