# Using Integrate.py

python integrate.py herx1_pulsar10on

will integrate the model defined in herx1_pulsar10on over phase and output the data.  I have written herx1_pulsar10on.py so that it will look for ximpol, but if it fails to find it, it will still work with integrate.py.

It outputs the following data in columns:

Energy in keV
Mean Flux
Mean Q
Mean Q/I

We use the rotation axis as the point of reference of the polarization.
