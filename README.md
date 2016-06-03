# QEDSurface

pfield will calculate the outgoing intensity in Stokes's Q using
projection of the star's the magnetic moment into the plane of
the sky as the reference direction under the assumption that at
the point of emission the intensity is entirely in the X mode
with unit value.

It spaces the image evenly so that one can calculate the total
polarized flux by simply summing the product of resulting Q values
with the (X-O) intensity at the point of emission.

It also calculates the location in magnetic colatitude of where
the photon is emitted at the given radius as well as the initial
zenith angle and azimuth with respect to the magnetic pole.

It accounts for light bending and gravitational redshift
in the Schwarzschild metric even if the radius is less than the
photon orbit.

  Format:

     pfield _mu_ _nu_ _mass_ _radius_ _alpha_ [optional parameters]
 
   _mu_     magnetic dipole moment in G cm^3
   _nu_     photon frequency in Hz [observer's frame]
   _mass_   mass in cm, i.e. GM/c^2
   _radius_ radius in cm
   _alpha_  angle of magnetic moment with line of sight in degrees
 
   [optional parameters]
 
            --doall   calculate for the top and bottom of the image
                      (for non-symmetric fields)
            --nb      number of impact parameters to calculate (10)
            optional spectral models to use

