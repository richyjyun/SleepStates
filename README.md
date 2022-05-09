# SleepStates ([bioRxiv](https://www.biorxiv.org/content/10.1101/2022.04.28.489967v1))

Tracking intracortical signals from the motor cortex of non-human primates via a Utah array for 20-24 hours and classifying brain states between 1) awake and moving, 2) awake and at rest, 3) REM sleep, 4) non-REM sleep. 

States were classified via dimensionality reduction and clustering. I first implemented a PCA (SVD) oriented approach, but ran into the limitations of linear transformations. I then implemented a stacked sparse autoencoder (inspired by Tsinalis et al. 2016) with the PyTorch package which was more consistent compared to EOG signals and movements recorded with the Kinect. The autoencoder method also showed clear REM cycles. Once the states were classified, spike and LFP dynamics throughout each state were analyzed. The "Classification" folder contains code performing state classification.  

Please see the latest manuscript for further details: https://www.biorxiv.org/content/10.1101/2022.04.28.489967v1

<p align="center">
  <img width="651.78" height="447.63" src="https://github.com/richyyun/SleepStates/blob/main/ClassificationFigure.png">
</p>

A. Diagram of classification process. LFP and accelerometer signals are used to perform PCA and k-means clustering to determine different sleep states. B. Example of classification and spectra over time.

Code flowchart

<p align="center">
  <img width="900" height="450" src="https://github.com/richyyun/SleepStates/blob/main/FlowChart.png">
</p>

Tsinalis, O., Matthews, P. M., & Guo, Y. (2016). Automatic Sleep Stage Scoring Using Time-Frequency Analysis and Stacked Sparse Autoencoders. Annals of Biomedical Engineering, 44(5), 1587–1597. https://doi.org/10.1007/s10439-015-1444-y
