# Forbidden_Transition
**[Kieran F. Thomas](https://github.com/KF-Thomas), [Jacob A. Ross](https://github.com/GroundhogState),[Bryce M. Henson](https://github.com/brycehenson)**

Analysis code for an experimental measurment of the 427 2<sup>3</sup>S<sub>1</sub> -> 3<sup>3</sup>S<sub>1</sub> highly forbidden transition in metastable helium.


| ![Level Diagram of He*](/figs/level_scheme_v4.svg "Fig2") | 
|:--:| 
| **Figure1**- Level diagram of He\* showing the forbidden 427nm transition, modified from [J. Simonet, Optical traps for Ultracold Metastable Helium atoms. PhD thesis](https://tel.archives-ouvertes.fr/tel-00651592/file/Simonet_PhD_Thesis.pdf)   |



## TO Do
- [ ] include core BEC as submodule for easy install
- [ ] thing 2


## Install
``` 
git clone --recurse-submodules -j8 https://github.com/KF-Thomas/Forbidden_Transition.git
```
then to update 
```
git submodule update --init --recursive --remote --merge
```

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
