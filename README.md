This repository is used to parameterize a second-order equivalent circuit model.

`idtfRCParam_IGWO.m`: the main script using improved gray wolf optimizer (IGWO) to find optimal ECM parameters.

`simECM.m`: a function to simulate the second-order ECM.

`buildModel.m`: a script to establish the structure named "model" in `idtfRCParam_IGWO.m`

Folder `IGWO` contains source code of the IGWO algorithm.

Folder `modelData` contains the structure "model" of each kind of battery.

Folder `optParams`: to save the optimized results after running `idtfRCParam_IGWO.m`

Folder `SOC-OCV_function` contains functions to calculate a OCV from a SOC, or vice versa.

Folder `tools` contains other tools for data process.