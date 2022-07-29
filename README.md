# SleepStates ([bioRxiv](https://www.biorxiv.org/content/10.1101/2022.04.28.489967v2))

Tracking intracortical signals from the motor cortex of non-human primates via a Utah array using the Neurochip3 (Shupe, 2021) for 20-24 hours and classifying brain states between 1) awake and moving, 2) awake and at rest, 3) REM sleep, 4) non-REM sleep. Once the states were classified, spike and LFP dynamics throughout each state were analyzed. The "Classification" folder contains code performing state classification.  

States were classified via dimensionality reduction using a stacked sparse autoencoder (inspired by Tsinalis et al. 2016) with the PyTorch package, then clustered with with k-means clustering. 

Please see the latest manuscript for further details: https://www.biorxiv.org/content/10.1101/2022.04.28.489967v2

Shupe, L. E., Miles, F. P., Jones, G., Yun, R., Mishler, J., Rembado, I., Murphy, R. L., Perlmutter, S. I., & Fetz, E. E. (2021). Neurochip3: An autonomous multichannel bidirectional brain-computer interface for closed-loop activity-dependent stimulation. Frontiers in Neuroscience, 15(August), 1–15. https://doi.org/10.3389/fnins.2021.718465

Tsinalis, O., Matthews, P. M., & Guo, Y. (2016). Automatic Sleep Stage Scoring Using Time-Frequency Analysis and Stacked Sparse Autoencoders. Annals of Biomedical Engineering, 44(5), 1587–1597. https://doi.org/10.1007/s10439-015-1444-y

<p align="center">
  <img width="651.78" height="447.63" src="https://github.com/richyyun/SleepStates/blob/main/ClassificationFigure.png">
</p>

A. Diagram of classification process. LFP and accelerometer signals are used to perform PCA and k-means clustering to determine different sleep states. B. Example of classification and spectra over time.

Code flowchart

<p align="center">
  <img src="https://github.com/richyyun/SleepStates/blob/main/CodeFlowChart.png">
</p>

## Analyses Performed
- Obtain the LFP power spectral density of every 8-second time bin using Welch's estimate for 0 to 50 Hz. Also calculate all pairwise magnitude-squared coherence.
- Sort spikes. 
  - I manually sorted each spike using two-window discrimination on ten minutes of data in the middle of the night then applied the determined parameters to the rest of the duration.
  - The fidelity of the spikes were calculated using the coefficient of determination (CoD) between spike waveforms. If the CoD between the first 1000 and last 1000 spikes is greater than the CoD of the first 1000 spikes with 1000 spikes from another recorded neuron, the spike was considered to have been stable overnight.
- Classify different brain states. 
  - Dimensionality reduction. I initially used PCA for dimensionality reduction which was able to detect sleep vs awake, but could not properly discern REM cycles and often missed naps during the day. I then applied a stacked autoencoder inspired by Tsinalis et al. 2016. Differences from the Tsinalis implementation include different normalization of the data (I normalize the spectra for each time bin rather than for each frequency to preserve the relative power between frequencies), smaller encoder/decoder (8 vs 20 layers), batch normalization applied to each intermediate layer, and ReLU activation function rather than sigmoid. 
  - Clustering. I used k-means clustering with 4 clusters. Rather than use all the data to find centroids, I calculated the average pairwise Euclidean distance between each data point and used all points under the 90th percentile of distances to avoid the influence of outliers. 
  - Majority filter. I applied a majority filter over ~30 seconds of classification to denoise the classification, remove outliers, and provide a temporal component to the classification.
- Confirm the classification using k-fold cross validation, by comparing it to body movement detected with the Kinect and eye movements detected with EOG electrodes, and observing changes in spike firing rates between states.
- Calculate cross-frequency phase-amplitude coupling between the phase of each lower frequency and amplitude of each higher frequency. I used mean vector length (MVL, Canolty et al., 2006) which transforms each phase-amplitude pair into a vector. The length of the average of the vector shows the strength of synchrony as larger vectors would mean the amplitudes are aligned to specific phases. 
  - Although the MVL has shown to be accurate and sensitive to modulations in coupling strength, it is also heavily influenced by the amplitudes. As the frequency components have different amplitudes in different states, there is a need to normalize the measure for comparisons. As a result I found the highest possible MVL by aligning the most common phases with the highest amplitudes and divided the actual MVL by it as a normalized MVL (nMVL).
  - Canolty, R. T., Edwards, E., Dalal, S. S., Soltani, M., Nagarajan, S. S., Kirsch, H. E., Berger, M. S., Barbare, N. M., & Knight, R. T. (2006). High gamma power is phase-locked to theta oscillations in human neocortex. Science, 313(5793), 1626–1628. https://doi.org/10.1126/science.1128115
- Determine changes in spiking dynamics. I calculated changes in the inter-spike intervals (ISIs) of each spike depending on the state to see how the dynamics change. I also separated the spikes in fast- or regular-spiking using the ISI distributions. 
- Calculate spike-LFP synchrony. I used phase-locking value (PLV) to calculate spike-LFP synchrony for each frequency band. PLV is calculated by creating unit vectors for each phase of the LFP at spike time and averaging the vectors together. 
