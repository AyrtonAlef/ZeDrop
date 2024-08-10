# ZeDrop


## ZeDropSIM
Simulation is important for understanding drop behavior in different configurations, for evaluating the effects of process parameters on their shape, and for verifying and validating new analysis routines for measuring surface properties. In the latter case, it is essential that images of the drop can be generated and sources of error present during the measurement process can be reproduced and isolated in order to evaluate the accuracy and robustness of the routines developed for evaluating real drop configurations.

Within this scope, ZeDropSIM was developed for verifying and simulating different drop configurations, such as pendant, sessile, and inclined drops. The program has several functionalities, such as two-dimensional (2D) and three-dimensional (3D) drop profiles can be generated, binary and grayscale images can be simulated, and sequences of consecutive images can be created to simulate quasi-static experiments. In addition, several sources of error can be added to the simulated image like lack of sharpness, lack of contrast, non-uniform illumination, random disturbance in the drop profile, vertical misalignment of the camera, uniform noise, Gaussian noise, impulse noise and satellite droplets. Figure 1 presents a flowchart that summarizes, according to the available drop configurations, all the functionalities offered by ZeDropSIM.

![Figure 1. ZeDropSIM flowchart](FlowchartZeDropSIM.png)
**Figure 1**. ZeDropSIM flowchart

ZeDropSIM was developed in Matlab. All the routines are present in ZeDropSIM folder. To start ZeDropSIM, MainCode_ZeDropSIM.m file must be run. 
