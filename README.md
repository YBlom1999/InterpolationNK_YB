# InterpolationNK_YB
This repository contains the implementation the interpolation method developed by Y. Blom et al. 
This implementation allows you to predict the complex refractive index of perovskite for a continious range of bandgap energies, based on the measured refractive index of a few samples.
The full details of this implementation are described in [1].

To use this implementation, either the matlab or the python script can be used. 
Both scripts contain a block that is called user inputs. This sections allows the user to change the used data files, the degree of the parameter fitting, the number of considered oscillators, and the desired bandgap energies for which the refractive index will be predicted.
The measured refractive index for a different bandgap energies should be stored in the data-folder as csv-file.
As example, the measured complex refractive index of three perovskite samples can be used, which are taken from Manzoor et al. [2].

If you use either one of this scripts, please refer to article [1].

If there are any questions regarding the method or the implementation, you can send them to y.blom@tudelft.nl.

Hopefully, this method can help you further in your research!

With kind regards,

Youri Blom





[1] Youri Blom, Malte Ruben Vogt, Olindo Isabella, and Rudi Santbergen, "Method for bandgap interpolation of perovskite’s spectral complex refractive index," Opt. Express 32, 4365-4375 (2024), https://doi.org/10.1364/OE.509982

[2] Salman Manzoor, Jakob Häusele, Kevin A. Bush, Axel F. Palmstrom, Joe Carpenter, Zhengshan J. Yu, Stacey F. Bent, Michael D. Mcgehee, and Zachary C. Holman, "Optical modeling of wide-bandgap perovskite and perovskite/silicon tandem solar cells using complex refractive indices for arbitrary-bandgap perovskite absorbers," Opt. Express 26, 27441-27460 (2018), https://doi.org/10.1364/OE.26.027441
