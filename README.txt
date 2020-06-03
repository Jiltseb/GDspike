This folder contains scripts to execute the GDspike algorithm for spike estimation from calcium fluorescence signals. This algorithm exploits the high-resoluytion property of group delay to enhance the peaks of the signal and is appropriately post-processed using a triangulation step. The details of the algorithm can be found at:

Jilt Sebastian,  Mari Ganesh Kumar M, Y. S. Sreekar, Rajeev Vijay Rikhye, Mriganka Sur and Hema A. Murthy, "GDspike: An Accurate Spike estimation Algorithm from Noisy Calcium Flourescence Signals", in ICASSP 2017, DOI: 10.1109/ICASSP.2017.7952315

and : 
https://ieeexplore.ieee.org/document/8680673

J. Sebastian, M. G. Kumar, V. S. Viraraghavan, M. Sur and H. A. Murthy, "Spike Estimation From Fluorescence Signals Using High-Resolution Property of Group Delay," in IEEE Transactions on Signal Processing, vol. 67, no. 11, pp. 2923-2936, 1 June1, 2019, doi: 10.1109/TSP.2019.2908913.

Usage:
Download the folder "GDspike". You would see some binaries and some matlab scripts. "GDspike.m" is the main function which calls "WordSegmentWithSilenceRemoval" binary for computation of group delay. Usage of the GDspike function is as follows:

[group_delay_output, triangulation_output, discrete_spike_estimates]= GDspike(fluorescence signal);

"input" refers to a Calcium fluorescence recording. The parameters such as the sampling rate of the signal, required resampling rate for spike estimation, threshold values etc. are provided inside the main function and can be edited appropriately.

Demo:
One example file for computation is provided in the file input.txt (Calcium fluorescence signal-60Hz) and spike.txt (Ground truth spike-10kHz resampled to 60 Hz) taken from publicly avilable GCaMP6f dataset (Chen. et al. 2013) at cncrs website.

input=importdata('input.txt');
spike=importdata('spike.txt');
[group_delay_output, triangulation_output, discrete_spike_estimates]= GDspike(input);

Plotting the system output:

subplot(5,1,1)
plot(input)
subplot(5,1,2)
plot(grp_delay)
subplot(5,1,3)
plot(signal)
subplot(5,1,4)
plot(GD_out)
subplot(5,1,5)
plot(spike)

Comments:

You could use the triangulation output as a analog spike inference signal or use discrete_spike_estimate as a discrete spike estimate. Tune parameters, especially the threshold for creating the discrete spike train. Evaluations can be done with respect to well-known measures such as correlation, auc and F-measure. 

If you are using GDspike for research, please cite our paper at ICASSP 2017 and the IEEE transactions on signal processing, FEB 2019.
