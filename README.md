# anisotropicScattering

Source code for paper Gerber et al., (2024): *Anisotropic Scattering in Radio-Echo Sounding: Insights from Northeast Greenland*

data files can be found at Zenodo XXX repository. 

### How to cite:


### script overview

#### Simulate azimuthal power response from EastGRIP COF with birefringence, scattering, and both: `EgripCOF_model.m`
* uses EGRIP eigenvalues from Zeising et al, 2022 (in EGRIP_cAxis.csv) and the fujita model (Fujita et al,. 2006) to simulate azimuthal response at EGRIP
* dependencies: `fujitaModel.m`, `prepCofInput.m`, `calculatePowerRatio.m`, `computePowerAnomalies.m`, `AverageDepth.m`, `EGRIP_cAxis.csv`
* generates Fig.3

#### Compare the modeled results with turning circle and synthesized response from quad-pol measurement near the turning circle: `obsSyntModel_comp.m` 
* comparison between azimuthal response from turning circle, synthesized and Fujita model; 
* dependencies: `computePowerAnomalies`, `computePhaseDerivative.m`, `prepCofInput.m`, `calculatePowerRatio.m`, `fujitaModel.m`
* generates Fig. 4 & Fig. 5

#### Calculate synthetic azimuthal response in 5km intervals along radar lines 
**Step 1: `prepareProfiles.m`**
* determining synthetic response at 5km intervals for lines with *linenum* as input for name of radar profile.
* dependencies: `readRadar.m`, `computePowerAnomalies.m`, `AverageDepth.m`, `computePhaseDerivative.m`
* reads .h5 radar files and saves synthetic response of analysis points from each line segment as .mat file.

**Step 2: `profilesCurveFitting.m`**
* calculate orientation and strength of scattering/birefringence vs depth
* use mode = 0 for 5 depth intervals, use mode = 1 for high-resolution depth calculation
* reads .mat files with analysis points and saves fitting parameters in .mat file.
* dependencies: input files, e.g. `profile20220701__091615__20trace_average.mat`

**plotProfileFig.m**
* plot Fig. 6 and Fig. 7, and Figures in appendix.
* mode = 0 for plotting the azimuthal response for every other analysis points (e.g. Fig. 6, Fig. 7)
* mode = 1 for plotting amplitudes and r^2 for high-res curve fitting
* dependencies: `saveCSVfile.m`

