# Introduction

This repository contains the macros used to estimate the background coming from the misidentification  of the electron charge in the LHC data with the ATLAS detector. 

While analyzing the data taken by the LHC, where the final state of interest is defined by two leptons with the same electric charge, the background coming from [Standard Model (SM)](https://en.wikipedia.org/wiki/Standard_Model) processes is very small. Therefore, it is crucial to consider the detector-related backgrounds: *objects that are misidentified or misreconstructed such that they appear to have a same-sign dilepton final state.* 

Opposite-sign leptons from SM processes as Drell-Yan, $W^+ W^-$, and mainly $t\bar t$, could contribute to the same-sign dilepton background if the charge of one of the leptons coming from the dileptonic decay of these processes  is mismeasured. This process is called *charge misidentification*.  The contribution of this background to the same-sign dilepton signature is estimated by measuring the probability that the lepton charge is misreconstructed using a data-driven technique.  

There are two main sources of electron charge misidentification: 

- Hard Bremsstrahlung producing trident electrons ($e^\pm \rightarrow e^\pm \gamma^* \rightarrow e^\pm e^+e^-$) whose electromagnetic cluster is identified with the wrong electron's track, leading to a misidentification of the charge. This source represents the main contribution to the background. The fraction of trident electrons depends on the amount of material that the electrons traverse. In the ATLAS detector, the distribution of the material depends on $|\eta|$. Therefore, a strong dependence on $|\eta|$ is expected in the misidentification rates.
- A slightly curved track that induces a measurement error. This effect is important at high transverse momentum ($p_\text{T}$).  Thus, a small dependence on electron $p_\text{T}$ is also expected in the misidentification rates.


### [Preliminary concepts](#anchor-rates) 

When a true opposite-sign event (for example $t\bar{t} \rightarrow bW^+\bar{b}W^- \rightarrow
b\bar{b}e^+e^-\nu\bar{\nu}$) is produced, and assuming that $\epsilon$ is the rate of charge misidentification for a single electron, there are three possibilities for this event to be reconstructed:

- $e^+e^- +X$ without any charge misidentification, with a probability of $(1-\epsilon)^2$,
- $e^+e^- +X$ with the two electrons having a charge flip, with a probability of $\epsilon^2$,
- $e^{\pm}e^{\pm} +X$ when only one of the two electrons is misidentified, with a probability of $2\epsilon(1-\epsilon)$.


Therefore, if there are $N$ true opposite-sign events, the reconstructed events will be:

- $N^{os} = (1-2\epsilon+2\epsilon^2) N$ opposite-sign events,
- $N^{ss} = 2\epsilon(1-\epsilon) N \simeq 2\epsilon N$ same-sign events,

where the last approximation for $N^{ss}$ corresponds to the assumption that $\epsilon$ is
very small.

Knowing this charge misidentification rate $\epsilon$, it is therefore possible to compute the
estimated number of same-sign events $N^{ss}$ from the measured number of opposite-sign
events $N^{os}$, using the following expressions:


 $$
 N^{ss} = \frac{\epsilon_i +\epsilon_j -2\epsilon_i \epsilon_j}{1-\epsilon_i -\epsilon_j +2\epsilon_i \epsilon_j} N^{os}  \tag{1}
 $$

for the $ee$ channel, and

 $$
 N^{ss} = \frac{\epsilon}{1-\epsilon} N^{os}  \tag{2}
 $$ 
 
for the $e\mu$ channel, where $\epsilon_i$ and $\epsilon_j$ are the charge misidentification rates for the two
different electrons.
    

# Description

The full process is divided into several steps:

1. **Sample selection:** The rates are estimated in a  $Z\rightarrow e^+e^-+jets$  sample for a specific electron selection.
2. **Estimate the charge misidentification rates**:  They are estimated using the likelihood method as a function of $|\eta|$ and $p_\text{T}$ of the electron (to be more precise, the full $|\eta|$ and $p_\text{T}$ ranges have been divided into regions, named *bins*, so that the rates are determined as a function of $|\eta|$ and $p_\text{T}$ bins).
5. **Estimate the background coming from charge misidentification:** These rates are then applied to scale the data but having opposite-sign electrons, which provide the expected background contribution for the same-signal final state. 

Details about how those steps are performed are shown below.

## Sample selection

The sample is selected by requiring the data to contain only events coming from the decay $Z\rightarrow e^+e^-+jets$, i.e. only events with two electrons and no muons are selected. In addition,  electrons must satisfy a set of minimum data quality criteria.  

This is done via the macro:

    Skimming.C

This macro also requires input the variables that will be saved in the final sample.  The macro will output the *sample* saved in  a .root file containing all the variables requested after applying the selection described above. 


## Estimate the charge misidentification rates

The misidentification rates of the electron charge are estimated using the likelihood method. 

The likelihood method  assumes that  the misidentification rates of the electron charge are independent for different $|\eta|$ regions. Therefore, the probability of having a number of same-sign events ($N^{ij}_{ss}$) with electrons in $|\eta|$ region $i$ and $j$ can be written as a function of the number of events $N^{ij}$ as follows:


$$
N^{ij}_{ss}=N^{ij}(\epsilon_i+\epsilon_j).
$$


If all the same-sign events in the $Z$ peak are produced by charge flip, then $N^{ij}\_{ss}$ is described by a Poisson distribution:

$$
f(k,\lambda)=\frac{\lambda^k e^{-\lambda}}{k!},
$$


where $k$ is the observed number of occurrences of the event, i.e. $k=N^{ij}\_{ss}$, and $\lambda$ is the expected number,  i.e. $\lambda=N^{ij}(\epsilon_i+\epsilon_j)$. Thus, the probability for both electrons to produce a charge flip is expressed by:


$$
P(\epsilon_i,\epsilon_j| N^{ij}\_{ss},N^{ij}) = \frac{[N^{ij}(\epsilon_i+\epsilon_j)]^{N_{ss}^{ij}}e^{-N^{ij}(\epsilon_i+\epsilon_j)}}{N^{ij}\_{ss}!}.
$$


The likelihood $L$ for all the events is obtained by evaluating all the $|\eta|$ combinations:


$$
L(\epsilon|N_{ss},N)=\prod_{i,j}\frac{[N^{ij}(\epsilon_i+\epsilon_j)]^{N_{ss}^{ij}}e^{-N^{ij}(\epsilon_i+\epsilon_j)}}{N^{ij}\_{ss}!},
$$


where the rates $\epsilon_i$ and $\epsilon_j$ can be obtained by minimizing the likelihood function. In this process, the $-\ln L$ is used to simplify and make easier the minimization. Terms which do not depend on the rates $\epsilon_i$ and $\epsilon_j$ are removed in this step. This way, the final function to minimize is given by the following expression:


$$
-\ln L(\epsilon|N_{ss},N)\approx \sum_{i,j}\ln[N^{ij}(\epsilon_i+\epsilon_j)]N^{ij}_{ss}-N^{ij}(\epsilon_i+\epsilon_j).
$$

The events are selected requiring the invariant mass of the two electrons to be compatible with the one of the  $Z$-boson and stored --with the electron order by $|\eta|$-- in two triangular matrices: one for the same-sign events $N^{ij}_{ss}$,  and the other one for all events $N^{ij}$. The likelihood method takes into account
electron pairs with all $|\eta|$ combinations, which allows to use of the full available statistics  getting therefore lower statistical uncertainties. Moreover, it does not bias the kinematical properties of the electrons, compared to other methods like tag-and-probe.

The likelihood  method can be easily extended to measure the charge misidentification rates as a function of  two parameters. In this study, the interest lies not only on the measurement of the rates   as a function of the $|\eta|$, but also on $p_\text{T}$. Thus, the probability of finding a same-sign event given the rates for each electron is ($\epsilon_{i,k}+\epsilon_{j,l}$), where the two indices represent binned $|\eta|$- and $p_\text{T}$-dependence. Thus, the last equation transforms into

$$
-\ln L(\epsilon|N_{ss},N)\approx \sum_{i,j,k,l}\ln[N^{ij,kl}(\epsilon_{i,k}+\epsilon_{j,l})]N^{ij,kl}\_{ss}-N^{ij,kl}(\epsilon_{i,k}+\epsilon_{j,l}) \tag{3} 
$$

The likelihood method uses only $Z$ *signal* events. Therefore, background coming from other processes where the dilepton invariant mass corresponds to the one of the $Z$ boson needs to be subtracted. The background subtraction is done using a simple *side-band method*.   This method consists in dividing the $Z$ invariant mass in three regions, i.e. $A$, $B$ and $C$, where $B$ is the central region corresponding to the $Z$ peak. The number of events is counted in the regions on the sides of the peak, i.e. $n_A$ and $n_C$, and removed  from the total number of events in the peak region $B$, $n_B$. This way, the number of signal events $N_Z$ is given by

$$
N_Z=n_B-\frac{n_A+n_C}{2}.
$$

Once the background has been subtracted, the likelihood method can be applied. 

All the procedure described in this section, including the minimization of the likelihood function defined by Equation (3) is performed by the macro:

    Likelihood.C

The macro will output the *charge misidentification rates that minimize the likelihood* as a function of the electron $|\eta|$ and $p_\text{T}$. Those rates are stored in 2D histograms, and saved into a file in a .root format. 


## Estimate the background from charge misidentification

Once the charge misidentification rates are estimated, every event in the opposite-sign sample is scaled by applying a weight computed using Equations (1) and (2). 

This is done with the following macro:

    GetQMisIDNtuple.C

For the input sample, the macro will output the *same sample with an extra variable* called *weight_qFlip*, corresponding to the charge flip weight. 


# Additional studies

## Compute Dependencies

ComputeDependencies.C


## Batch processing

ChargeFlipBatchNoCF.C

SubmitChargeFlipRates.py

SubmitValidationPlots.py

# Official results obtained with this tool




