# anisotropicScattering

Source code for paper Gerber et al., (2024): *Anisotropic Scattering in Radio-Echo Sounding: Insights from Northeast Greenland*

data files can be found at Zenodo XXX repository. 

### How to cite:


### script overview:

**EgripCOF_model.m**: 
* uses EGRIP eigenvalues from Zeising et al, 2022 (in EGRIP_cAxis.csv) and the fujita model (Fujita et al,. 2006) to simulate azimuthal response at EGRIP
* dependencies: `fujitaModel.m`, `prepCofInput.m`, `calculatePowerRatio.m`, `computePowerAnomalies.m`, `AverageDepth.m`, `EGRIP_cAxis.csv`
* generates Fig.3

**obsSyntModel_comp.m** --> reproduce Fig. 4 & Fig. 5
* comparison between azimuthal response from turning circle, synthesized and Fujita model
* dependencies: 
* generates Fig. 4 & Fig. 5

**prepareProfiles.m**
* determining synthetic response at 5km intervals for lines
* dependencies:
* 
  
**profilesCurveFitting.m**
* calculate orientation and strength of scattering/birefringence
* dependencies:

**plotProfileFig.m**
* 
* 

