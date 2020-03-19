# Forbidden_Transition
**[Kieran F. Thomas](https://github.com/KF-Thomas),  [Jacob A. Ross](https://github.com/GroundhogState),  [Bryce M. Henson](https://github.com/brycehenson)**

Analysis code for an experimental measurment of the 427 nm 2<sup>3</sup>S<sub>1</sub> -- 3<sup>3</sup>S<sub>1</sub> highly forbidden transition in metastable helium.


| ![Level Diagram of He*](/figs/level_scheme_v5.png "Fig1") | 
|:--:| 
| **Figure 1**- Level diagram of He\* showing the forbidden 427nm transition.  |


A write up of the results of this experiment can be found at https://arxiv.org/abs/2002.04811.

| ![Direct detection scan*](/figs/direct_scan.png "Fig2") | 
|:--:| 
| **Figure 2**- Normalised detected number of scattered atoms against probe frequency, Fig. 3 in https://arxiv.org/abs/2002.04811.  |

The results of the direct detection method showin in Figure 2 are produced with the file '*fig_direct_scan.m*', while the center frequency can be calculated with '*calc_freq.m*'. The Einstein A coefficent fro the 2<sup>3</sup>S<sub>1</sub> -- 3<sup>3</sup>S<sub>1</sub> can be calculated from our data using the file '*calc_A.m*', and the excited state lifetime of the 3<sup>3</sup>S<sub>1</sub> state can be derived from our data using the file '*calc_tau.m*'.

## Install
``` 
git clone --recurse-submodules -j8 https://github.com/KF-Thomas/Forbidden_Transition.git
```
then to update 
```
git submodule update --init --recursive --remote --merge
```

## Running
The code only requires MATLAB to run. Note that it does need all subdirectories to be added to the current path in order to exectue properly.

## Contributions  
This project would not have been possible without the many open source tools that it is based on. In no particular order: 

* ***James Conder*** [gaussfilt](https://au.mathworks.com/matlabcentral/fileexchange/43182-gaussfilt-t-z-sigma)
* ***Ander Biguri*** [Perceptually uniform colormaps](https://au.mathworks.com/matlabcentral/fileexchange/51986-perceptually-uniform-colormaps)
* ***Jan*** [FileTime](https://au.mathworks.com/matlabcentral/fileexchange/24671-filetime)
* ***Benjamin Kraus*** [nanconv](https://au.mathworks.com/matlabcentral/fileexchange/41961-nanconv)
* ***M. A. Hopcroft**** [allan](https://au.mathworks.com/matlabcentral/fileexchange/13246-allan)
* ***Daniel Eaton***  [sfigure](https://au.mathworks.com/matlabcentral/fileexchange/8919-smart-silent-figure)
* ***Denis Gilbert***  [M-file Header Template](https://au.mathworks.com/matlabcentral/fileexchange/4908-m-file-header-template)
* ***DrosteEffect***  [CIECAM02](https://github.com/DrosteEffect/CIECAM02)
* ***Holger Hoffmann*** [violin](https://au.mathworks.com/matlabcentral/fileexchange/45134-violin-plot)
