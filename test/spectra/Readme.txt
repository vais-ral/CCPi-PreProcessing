These directories store spectra for X-ray target materials.
 - spectra/W contains the spetra for tungsten. In this case the source is spekcalc as run at QMUL.
   data is given for energies in range 1KeV to 150KeV in steps of 1KeV and take-off angles between
   0 and 89 degress. Files are named "eeeaaa.spc" where the "eee" values are the KeVs, e.g. 091 for 91KeV
   and the "aaa" is angle*10, e.g. 420 for 42 degrees.
   This data is more detailed than needed.
 - spectra/Mo contains a few Molybdenum spectra. These are from a Siemens's web site and are not take-off
   angle corrected. To be consistent with the W naming we assume an angle of 45 degrees. Any take off angle
   dependence needs to be approximated by the fitting process, as is done in practice for W anyway.
   The normalisation of Mo and W data are very different.
