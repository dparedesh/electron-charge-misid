# Introduction

This repository contains the macros used to estimate the background from misidentifying the electron charge in the LHC data with the ATLAS detector. 

While analyzing the data taken by the LHC where the final state of interest is defined by two leptons with the same electric charge, the background coming from [Standard Model (SM)](https://en.wikipedia.org/wiki/Standard_Model) processes is very small. Therefore, it is crucial to consider the detector-related backgrounds: *objects that are misidentified or misreconstructed such that they appear to have a same-sign dilepton final state.* 

Opposite-sign leptons from SM processes as Drell-Yan, $W^+ W^-$, and mainly $t\bar t$, could contribute to the same-sign charge dilepton background if the charge of one of the leptons coming from the dileptonic decay of these processes  is mismeasured. This process is called *charge misidentification*.  The contribution of this background to the same-sign dilepton signature is estimated by measuring the probability that the lepton charge is misreconstructed using a data-driven technique.  

Here are two main sources of electron charge misidentification: 

- Hard Bremsstrahlung producing trident electrons ($e^\pm \rightarrow e^\pm \gamma^* \rightarrow e^\pm e^+e^-$) whose electromagnetic cluster is identified with the wrong electron's track, leading to a misidentification of the charge. This source represents the main contribution to the background. The fraction of trident electrons depends on the amount of material that the electrons traverse. In the ATLAS detector, the distribution of the material depends on $|\eta|$. Therefore, a strong dependence on $|\eta|$ is expected in the misidentification rates.
- A slightly curved track that induces a measurement error. This effect is important at high transverse momentum ($p_\text{T}$).  Thus, a small dependence on electron $p_\text{T}$ is also expected in the misidentification rates.
    

# Description

The full process is divided into several steps:

1. **Sample selection:** The rates are estimated in a  $Z\rightarrow e^+e^-+jets$  sample for a specific electron selection.
2. **Estimate the charge misidentification rates**:  They are estimated using the likelihood method as a function of $|\eta|$ and $p_\text{T}$ of the electron (to be more precise, the full $|\eta|$ and $p_\text{T}$ ranges have been divided into regions, named *bins*, so that the rates are determined as a function of $|\eta|$ and $p_\text{T}$ bins).
5. **Estimate the background coming from charge misidentification:** These rates are then applied to scale the data but having opposite-sign electrons, which provide the expected background contribution for the same-signal final state. 

Details about how those steps are performed are shown below.

# Sample selection


Skimming.C

# Estimate the charge misidentification rates

Likelihood.C


# Estimate the background from charge misidentification

GetQMisIDNtuple.C

# Additional studies

## Compute Dependencies

ComputeDependencies.C


## Batch processing

ChargeFlipBatchNoCF.C
SubmitChargeFlipRates.py
SubmitValidationPlots.py

# Official results obtained with this tool




