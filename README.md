# Description of the main code
These are code and figures for several different things including random walk, GLM, tracking and pooling. Some functions are be used for different purpose. 
Some code is copied with modification to meet different purposes. 
So there are multiple code with similar names (e.g. DistHaz_XXX.m). I will introduce these "family" code together. 
There are some files I labeled as "unuseful". I believe they are some previous versions of code but I am not sure if some parameters are still valuable.
(to draw a specific case, etc) So I keep them in the folder as well. 

## Random walk neuron model simulations and analysis
### DistHaz_Ult.m
It simultes spike trains from random walk neuron model using function "GetISI.m" (means get interspike interval). \
Then it plots ISI histogram with a hazard function calculated by time bins. \
In the bottom of the same figure, it plots the random inputs, signal inputs, voltage and spike trains. \
After that, spike trains are fitted to three GLM models: 
- Fixed Length Filters (FLF) GLM (Pillow)\
Using "fit_glm.m"
- Hazard function (aka only the first spike)\
Using "fit_hazard.m"
- Fixed Number Filter (FNF) GLM (default is three spikes)\
Using "fit_FNL.m"

### DistHaz_CV.m
It returns the coefficient of variance of interspike interval (CV) as membrane time constant (\tau_m) and synaptic time constant (\tau_s) change. 

### DistHaz_MapPara.m
It plots the spike trains as membrane time constant (\tau_m) and synaptic time constant (\tau_s) change. 
The four extreme cases in the two dimensional parameter space: Poisson spiking, inverse-Gaussian spiking, isochronal spiking and "burst". 

### DistHaz_Mix.m
It is like "DistHaz_Ult.m", but there are four kinds of synapses: AMPA, NMDA for excitatory, GABA-A and GABA-B for inhibitory. 
The corresponding decay time constants are: 5ms,67ms,5ms,150ms. 

### DistHaz_PopInput.m
Unuseful. It's probably a early version for population input. The code for population input should be in the folder "\Code coupling filter\CouplingFilter.m". 

### DistHaz.m 
Unuseful. I think it should be an early version of "DistHaz_Ult.m". But maybe I left it there without deleting it because I still want to keep some parameters. 

### Copy_of_DistHaz_Ult.m
It's the same as DistHaz_Ult.m, but it should be useful. I believe I copied it to keep some parameters for a final figure. 

### Conductence_DistHaz_Ult.m
In the random walk model, voltage step sizes V_E and V_I are fixed. But relating the step size to reversal potential, random walk neuron model becomes just LIF model. 
I call it conductance model. Also, I showed that RW model can generate similar spike trains as conductance models if we choose a proper lower boundary. 

##  Code coupling filter
Instead of a specific current or a specific stimulus, the inputs here are a neuron population. Default number of neurons in the population is one hundred. 
All neuron has the same underlying firing rate. GLM with coupling filter is fitted with only a proportion of the neurons. It turns out that pooling improves the fitting. 
### Code coupling filter\CouplingFilter.m
The main code that illustate the results with a single trial. 
### Code coupling filter\plot_IMSE_nlogL.m
It shows fitted resuls (using IMSE and negative log likelihood to compare goodness of fit)
- with different number of trials
- with or without pooling. 
### Code coupling filter\plot_differentKS.m
Unuseful. It's the same as plot_IMSE_nlogL.m but compares the p values using KS test.  I meant to find a way that tells us how many trials we need. 

## Tracking signal
### TrackingSignal.m
It a lot like the previous one (Random walk neuron model simulations and analysis), spike trains are also simulated from "GetISI.m", 
but the signal input is either square waves or Gamma functions. 
After plotting a raster, it plots the firing rate from simulations along with input signal. MISE is used to compare output firing rate and input signal. 
Finally, it fits FLF-GLM and do KS test. (Optionally it can estimate hazard function for baseline input spike trains.)

### TrackingSignal_compare.m
Compare different definitions of SNR (and MISE) in Poisson process. Spike trains are under different number of inputs. 
There are five SNRs: 
- Rob's SNR $\frac{\int (\lambda-\bar{\lambda})^2 dt}{\bar{\lambda}}$
- Rob's SNR but with baseline firing rather than mean firing rate as the denominator
- An SNR defined by [Lesica et al](https://doi.org/10.3389/fncom.2010.00144) ， which actually converges to the first SNR
- A model based SNR (difference in the likelihood with and without stimulus effect)

### StudySNR_Refractory.m
Variance-based SNRs decreases as absolute refractory period gets longer. Also I plot the numerator and denominator in the SNR definition： 
Both numerator and denominator goes down, but numerator goes down faster. 

### StudySNR_TrialNumber.m 
The code stduies how number of trials influences different SNRs. With more trials, variance-based SNRs decrease while GLM-based SNRs remain the same. 

## LIFtoGLM.m 
It draw a figure to illustate balanced position in LIF can be viewed as log firing rate in GLM. In other words, firing rate ~ exp(mean voltage) 

## In Folder \Normalization
I tried to illustate divisive normalization with HH model, GLM and random walk models. 
But the result seems to be something between divisive normalization and subtractive normalization. 

## In Folder \Code mix receptors 
Unuseful. It's a previous version of "DistHaz_Mix.m". Only two types of receptors in considered in the code. 

## In Folder \Code lower threshold
DistHaz.m and DistHaz_CV.m are like the ones I described. But here it focus on studying the effects of lower boundary. 

# Description of the functions
Detailed descriptions are in the corresponding file
- sameconv.m \
get the convolution of post spike train and post-spike filter in FLF GLM. 
- sameconvSti.m \
get the convolution of stimulus and stimulus filter. 
- sameconv_FNL.m \
get the convolution of post spike train and post-spike filter in FNF GLM. 
- sameconv_OnlyLastOne.m
get the convolution of last spike train and post-spike filter to estimate hazard function. 
- plotraster.m
draw a raster of spike trains. 
- normalizecols.m
normalizes the columns of a matrix, so each is a unit vector. 
- makeBasis_StimKernel.m
get Pillow's basis for stimulus filters. 
- makeBasis_PostSpike.m
get Pillow's basis for post-spike filters. 
- KStest.m
conduct KS test on a fitted GLM. (return a figure and a p value of KS test with spike trains, stimulus filter, post-spike filter, input and biased rate)
- GetISI.m
simulate RW model with given parameters. 
- GetMeanISI_J.m
use this function to know if how far the mean ISI is away from the target mean. 
- GetISI_Mix.m
simulate RW model with multiple receptors. 
- get_features.m
returns the value and position of dip and peak of a filter. 

