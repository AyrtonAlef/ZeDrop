# ZeDrop


## ZeDropSIM
<p><div style="text-align: justify">Simulation is important for understanding drop behavior in different configurations, for evaluating the effects of process parameters on their shape, and for verifying and validating new analysis routines for measuring surface properties. In the latter case, it is essential that images of the drop can be generated and sources of error present during the measurement process can be reproduced and isolated in order to evaluate the accuracy and robustness of the routines developed for evaluating real drop configurations.</div></p>

<p><div style="text-align: justify">Within this scope, ZeDropSIM was developed for verifying and simulating different drop configurations, such as pendant, sessile, and inclined drops. The program has several functionalities, such as two-dimensional (2D) and three-dimensional (3D) drop profiles can be generated, binary and grayscale images can be simulated, and sequences of consecutive images can be created to simulate quasi-static experiments. In addition, several sources of error can be added to the simulated image like lack of sharpness, lack of contrast, non-uniform illumination, random disturbance in the drop profile, vertical misalignment of the camera, uniform noise, Gaussian noise, impulse noise and satellite droplets. <strong>Figure 1</strong> presents a flowchart that summarizes, according to the available drop configurations, all the functionalities offered by ZeDropSIM.</div></p>

![Figure 1. ZeDropSIM flowchart](FlowchartZeDropSIM.png)
<p><div align="center"><strong>Figure 1. ZeDropSIM flowchart. </strong></div></p>

<p>ZeDropSIM was developed in Matlab. All the routines are present in ZeDropSIM folder. To start ZeDropSIM, <em>MainCode_ZeDropSIM.m</em> file must be run. </p>

## ZeDropACT

<p><div style="text-align: justify">The implementation of optical methods for analyzing drops involves hardware and software requirements. The hardware requirements are related to the experimental apparatus needed to perform the tests and capture images. The software must be able to analyze the captured images, identify and extract the drop profile and, using drop shape analysis methods, determine surface properties, such as liquid surface tension and contact angle. Although commercial solutions exist for determining surface properties using optical methods, these are expensive and often inaccessible.</div></p>

<p><div style="text-align: justify">Given this scenario, the ZeDropACT solution was developed. ZeDropACT aims to provide an affordable solution for performing tests involving drops to obtain reliable measurements of surface properties using optical methods. ZeDropACT consists of independent hardware and software solutions. In terms of hardware, a low-cost, modular goniometer project capable of performing different tests was designed, while in terms of software, a program capable of controlling different process variables and following different test protocols was developed.</div></p>

<p><div style="text-align: justify">With ZeDropACT software, the user has access to adjustment and test routines. The adjustment routines allow fine adjustment of the different systems that make up the goniometer, while the test routines execute previously established procedures for image acquisition and subsequent determination of surface properties by optical methods.</div></p>

<p><div style="text-align: justify">It is possible to perform intermittent and continuous tests, involving pendant, sessile and inclined drops. It is important to highlight that for certain types of tests, such as quasi-static and dynamic, it is necessary to use electronic, mechanical and optical components that allow more refined and rigorous control of the process variables.</div></p>

<p><div style="text-align: justify">After starting the ZeDropACT software, a menu is displayed showing all the program's functionalities. Figure 2 shows a flowchart containing all the options offered. Options 1 to 8 are related to the adjustment routines, while options 9 to 14 are linked to the execution of the test routines.</div></p>

![Figure 2. ZeDropACT flowchart](FlowchartZeDropACT.png)
<p><div align="center"><strong>Figure 2. ZeDropACT flowchart. </strong></div></p>

<p><div style="text-align: justify">The ZeDropACT software does not contain drop shape analysis routines. The main result produced by the software is drop image. In order to identify the drop profile, analyze its shape and determine surface properties, it is necessary to use another solution. The ZeDropEVAL sofwtare can be used for this purpose.</div></p>

