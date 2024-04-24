# Interpolation of the complex refractive index w.r.t. the bandgap energy
This repository contains an implementation of the interpolation method developed by Y. Blom et al. 
It allows for predicting the complex refractive index of perovskite across a continuous range of bandgap energies, based on refractive index measurements from a few samples. 
For detailed information, refer to [1].

To utilize this implementation, you can use either the MATLAB or Python script provided. 
Both scripts feature a section labeled "user inputs," enabling users to adjust data files, degree of parameter fitting, number of considered oscillators, and desired bandgap energies for refractive index prediction. 
Measured refractive index data for various bandgap energies should be stored as CSV files in the data folder. 
For instance, you can use refractive index data from three perovskite samples sourced from Manzoor et al. [2] as examples.

When using these scripts, please cite article [1].

For questions about the method or implementation, feel free to reach out to y.blom@tudelft.nl.

We hope this method is helpful in your research!

Best regards,

Youri Blom





[1] Youri Blom, Malte Ruben Vogt, Olindo Isabella, and Rudi Santbergen, "Method for bandgap interpolation of perovskite’s spectral complex refractive index," Opt. Express 32, 4365-4375 (2024), https://doi.org/10.1364/OE.509982

[2] Salman Manzoor, Jakob Häusele, Kevin A. Bush, Axel F. Palmstrom, Joe Carpenter, Zhengshan J. Yu, Stacey F. Bent, Michael D. Mcgehee, and Zachary C. Holman, "Optical modeling of wide-bandgap perovskite and perovskite/silicon tandem solar cells using complex refractive indices for arbitrary-bandgap perovskite absorbers," Opt. Express 26, 27441-27460 (2018), https://doi.org/10.1364/OE.26.027441
