# anisotropicScattering

Source code for paper Gerber et al., (2024): *Anisotropic Scattering in Radio-Echo Sounding: Insights from Northeast Greenland*

data files can be found at Zenodo XXX repository. 

### How to cite:


### script overview:

**Simulate azimuthal power response from EastGRIP COF with birefringence, scattering, and both: `EgripCOF_model.m`** 
* uses EGRIP eigenvalues from Zeising et al, 2022 (in EGRIP_cAxis.csv) and the fujita model (Fujita et al,. 2006) to simulate azimuthal response at EGRIP
* dependencies: `fujitaModel.m`, `prepCofInput.m`, `calculatePowerRatio.m`, `computePowerAnomalies.m`, `AverageDepth.m`, `EGRIP_cAxis.csv`
* generates Fig.3

**Compare the modeled results with turning circle and synthesized response from quad-pol measurement near the turning circle: `obsSyntModel_comp.m`** 
* comparison between azimuthal response from turning circle, synthesized and Fujita model; 
* dependencies: `computePowerAnomalies`, `computePhaseDerivative.m`, `prepCofInput.m`, `calculatePowerRatio.m`, `fujitaModel.m`
* generates Fig. 4 & Fig. 5

**Calculate synthetic azimuthal response in 5km intervals along radar lines: **
**Step 1: `prepareProfiles.m`**
* determining synthetic response at 5km intervals for lines
* dependencies:
* 
  
**profilesCurveFitting.m**
* calculate orientation and strength of scattering/birefringence
* dependencies:

**plotProfileFig.m**
* 
* 

