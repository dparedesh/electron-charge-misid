# Introduction

This repository contains the macros used to estimate the background coming from the misidentification  of the electron charge in the LHC data with the ATLAS detector. 

While analyzing the data taken by the LHC, where the final state of interest is defined by two leptons with the same electric charge, the background coming from [Standard Model (SM)](https://en.wikipedia.org/wiki/Standard_Model) processes is very small. Therefore, it is crucial to consider the detector-related backgrounds: *objects that are misidentified or misreconstructed such that they appear to have a same-sign dilepton final state.* 

Opposite-sign leptons from SM processes as Drell-Yan, $W^+ W^-$, and mainly $t\bar t$, could contribute to the same-sign dilepton background if the charge of one of the leptons coming from the dileptonic decay of these processes  is mismeasured. This process is called *charge misidentification*.  The contribution of this background to the same-sign dilepton signature is estimated by measuring the probability that the lepton charge is misreconstructed using a data-driven technique.  

There are two main sources of electron charge misidentification: 

- Hard Bremsstrahlung producing trident electrons ($e^\pm \rightarrow e^\pm \gamma^* \rightarrow e^\pm e^+e^-$) whose electromagnetic cluster is identified with the wrong electron's track, leading to a misidentification of the charge. This source represents the main contribution to the background. The fraction of trident electrons depends on the amount of material that the electrons traverse. In the ATLAS detector, the distribution of the material depends on $|\eta|$. Therefore, a strong dependence on $|\eta|$ is expected in the misidentification rates.
- A slightly curved track that induces a measurement error. This effect is important at high transverse momentum ($p_\text{T}$).  Thus, a small dependence on electron $p_\text{T}$ is also expected in the misidentification rates.


### Preliminary concepts

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
    

# Methodology

The full process is divided into several steps:

1. **Sample selection:** The rates are estimated in a  $Z\rightarrow e^+e^-+jets$  sample for a specific electron selection.
2. **Estimation of the charge misidentification rates**:  They are estimated using the *likelihood method* as a function of $|\eta|$ and $p_\text{T}$ of the electron (to be more precise, the full $|\eta|$ and $p_\text{T}$ ranges have been divided into regions, named *bins*, so that the rates are determined as a function of $|\eta|$ and $p_\text{T}$ bins).
5. **Estimation of the background coming from charge misidentification:** These rates are then applied to scale the data but having opposite-sign electrons, which provide the expected background contribution for the same-signal final state. 

Details about how those steps are performed are shown below.

## Sample selection

The sample is selected by requiring the data to contain only events coming from the decay $Z\rightarrow e^+e^-+jets$, i.e. only events with two electrons and no muons are selected. In addition,  electrons must satisfy a set of minimum data quality criteria.  

This is done via the macro:

    Skimming.C

This macro also requires input the variables that will be saved in the final sample.  The macro will output the *sample* saved in  a .root file containing all the variables requested after applying the selection described above. 


## Estimation of the charge misidentification rates

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


## Estimation of the background from charge misidentification

Once the charge misidentification rates are estimated, every event in the opposite-sign sample is scaled by applying a weight computed using Equations (1) and (2). 

This is done with the following macro:

    GetQMisIDNtuple.C

For the input sample, the macro will output the *same sample with an extra variable* called *weight_qFlip*, corresponding to the charge flip weight. 


# Additional studies...

The studies described above are done with *real data* from the LHC and they correspond to the standard methodology to compute the charge misidentification rates. However, samples of simulated events can be used to cross-check the results obtained with real data and to validate the method used to compute the rates. 

 Moreover, additional studies were also performed to check:

 1. If there are other dependencies apart from the ones on $|\eta|$ and $p_\text{T}$.
 2. If the rates for electrons are compatible with those from positrons.
 3. If the rates for electrons from photon conversion are compatible with the ones from wrong electron's track.
 4. If the rates are compatible for the different physics processes. Our main interest is to check the compatibility of the rates for  $Z+jets$  and $t\bar t+jets$.


## Truth matching method

In samples of simulated events, all the information related to the *true* generated electron is known. The charge of the *true* electron can be compared to the one of the *reconstructed* electron. Therefore, the misidentification rates can be computed as the ratio of the number of reconstructed electrons with their charge misidentified with respect to the total number of reconstructed electrons in the sample. This method is called  *truth-matching*. 

In addition, samples of simulated events also save the information related to the *origin* of the electron, i.e. is the electron coming from hard Bremsstrahlung? or was the curvature of the electron's track mismeasured? the misidentification rates can be computed for every one of them.  


### Estimation of the charge flip rates in parallel: batch jobs

The run time for estimating the charge misidentification rates can be long. If the charge flip rates are computed for different dependencies the most efficient way to do that is by running parallel jobs.  

The submission of the batch jobs to compute the charge flip rates using the truth-matching method for all the desired dependencies (in addition to $|\eta|$ and $p_\text{T}$) can be done using the script:

    SubmitChargeFlipRates.py

This script needs an input dictionary `conf[<key>] = <value>` with the ID assigned to the desired configuration as `key` and the selection to be applied to the data as `value`. For example,

   ```python
   conf["_nJets_ge5"]="nJets>=5"
   ```

The script calls internally the macro 

    TruthMatchingBatch.C

for every `key` in the `conf` dictionary, and will output the charge flip rates for the `value` provided. 

In addition, the script will output the rates computed for electrons and positrons, and also the charge flip rates depending on the *origin* of the electron. 

The procedure can be done for every physical process, i.e. $Z+jets$  and $t\bar t+jets$.  


### Charge flip rates for different dependencies

The ratios of the charge flip rates as a function of $|\eta|$, $p_\text{T}$, and parametrized on the different dependencies can be obtained with the macro:

    ComputeDependencies.C

This macro will output the plot for every variable for which the dependency wants to be study. It needs as input the rates previosly calculated. 


### Charge flip rates from different physics processes

The ratios of the charge flip rates for $Z+jets$  and $t\bar t+jets$ as a function of $|\eta|$ and $p_\text{T}$ are obtained with the following macro:

    CompareProcesses.C

It needs as input the rates previosly calculated for every process.


## Validation of the charge misidentification rates

The charge flip rates obtained previously are validated by comparing the number of same-sign events with the number of *weighted* opposite-sign events. Ideally, this is done by comparing the distributions of the different variables used in the analysis. As the samples size can be huge, the most effective way to do this is with parallel jobs.

The batch submission of the jobs to get the validation plots for every desired variable is with the script: 

    SubmitValidationPlots.py


This script needs as input a list of the `variables` to be plotted. For example:

```python
variables=["pt1","pt2","eta1","met","Ht","Htlep","Htjets","nJets","BJets","mu","pv","Mll"]
```

For every element in the list `variables`, the script will call the macro:
    
    ValidationPlots.C

which will output the 1D distribution of the same-sign and weighted opposite events overlaid in the same plot. 



# Official results obtained with this tool

This tool has been used to estimate the background coming from the charge misidentification of the following paper publications:


- [[1] JHEP 07 (2023) 203](https://inspirehep.net/literature/2175533): Search for $t\bar tH/A \rightarrow t\bar tt\bar t$ production in the multilepton final state in proton–proton collisions at $\sqrt{s}=13$ TeV with the ATLAS detector.
- [[2] JHEP 11 (2021) 118](https://inspirehep.net/literature/1869695): Measurement of the $t\bar tt\bar t$ production cross section in $pp$ collisions at $\sqrt{s}​=13$ TeV with the ATLAS detector.
- [[3] Eur. Phys. J. C 80 (2020) 1085](https://inspirehep.net/literature/1809244): Evidence for $t\bar tt\bar t$ production in the multilepton final state in proton–proton collisions at $\sqrt{s}
​=13$ TeV with the ATLAS detector.

- [[4] Phys.Lett.B 749 (2015) 519-541](https://inspirehep.net/literature/1377202): Search for the associated production of the Higgs boson with a top quark pair in multilepton final states with the ATLAS detector.

- [[5] JHEP 10 (2015) 150](https://inspirehep.net/literature/1361912): Analysis of events with $b$-jets and a pair of leptons of the same charge in $pp$ collisions at $\sqrt{s}=8$ TeV with the ATLAS detector.


<p align="center">
<a href="url"><img src="https://github.com/dparedesh/electron-charge-misid/assets/13987503/82b06d02-b3e3-4f6c-92dd-2090312b0c23" align="center" height="250"  ></a>
<a href="url"><img src="https://github.com/dparedesh/electron-charge-misid/assets/13987503/b868eac9-215c-4b48-a931-b126b30f30dc" align="center" height="265"  ></a>
</p>
<p align="center"><sub>Figure 1: Comparison of data and total background in the main discriminant variable in this analysis (left). The contribution coming from charge misidentification is labeled as *QmisID*.  Those results were used to compute upper limits on the production rate of the $t\bar tH \rightarrow t\bar tt\bar t$ signal  as a function of the mass of the $H/A$ particle (right)[1]. </sup></p>  



<p align="center">
<a href="url"><img src="https://github.com/dparedesh/electron-charge-misid/assets/13987503/0af12b9c-3c50-40c3-a0b9-5db07ffe9d27" align="center" height="250"  ></a>
<a href="url"><img src="https://github.com/dparedesh/electron-charge-misid/assets/13987503/5d31dd18-8c01-435e-8157-7ab2129cd80a" align="center" height="270"  ></a>
</p>
<p align="center"><sub>Figure 2: Comparison of data and total background in the main discriminant variable in this analysis (left). The contribution coming from charge misidentification is labeled as *Q mis-id*. Those results were used to compute the production rate of the $t\bar tt\bar t$  process in the same-sign dilepton channel [2], and its statistical combination with the 1L/2LOS channel (right)[3]. </sup></p>  



<p align="center">
<a href="url"><img src="https://github.com/dparedesh/electron-charge-misid/assets/13987503/43086def-dcb4-4ea9-a28f-90cc12016ba9" align="center" height="230"  ></a>
<a href="url"><img src="https://github.com/dparedesh/electron-charge-misid/assets/13987503/91980197-9f01-4c8c-888e-56cbf6343caa" align="center" height="215"  ></a>
</p>
<p align="center"><sub>Figure 3: Comparison of the charge misidentification background with the $t\bar tH$ signal (left). Those results were used to compute the signal strength $\mu$ of the $t\bar tH$ signal (right)[4]. </sup></p>  

<p align="center">
<a href="url"><img src="https://github.com/dparedesh/electron-charge-misid/assets/13987503/f2030758-2202-447f-857a-745b4db79500" align="center" height="270"  ></a>
<a href="url"><img src="https://github.com/dparedesh/electron-charge-misid/assets/13987503/eea88a3f-6234-4fed-8f33-aec9dc3e575d" align="center" height="270"  ></a>
</p>
<p align="center"><sub>Figure 4:  Comparison of data and total background in this analysis (left). The background coming from charge misidentification is labeled *Q Mis-id*. Those results were used to perform the statistical interpretation of the data for different New Physics models (right) of Ref. [5]. </sup></p>  



