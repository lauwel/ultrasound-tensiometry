Ultrasound-tensiometry

This GitHub repository contains the MATLAB code that enables users to capture ultrasound-tensiometry data (see manuscript: https://doi.org/10.1016/j.gaitpost.2024.06.026). The goal of ultrasound-tensiometry is to assess the spatial variation in shear wave speed in response to a mechanical tap. This repository contains two sections: "Collection" and "Processing".

Clone the repository directly into the Vantage ### folder on the Verasonics computer to collect data. The processing code does NOT require any of the Verasonics functions, and can be cloned onto any other computer. 

---REQUIRED HARDWARE ----
This system is designed to function with a Verasonics Vantage 64LE Ultrasound system and a Phillips CL15-7 ultrasound probe. Other probes and Vantage implementations (128, 256, etc) may be implemented in the future. We are actively working on a 128 implementation.

To induce the tap, the Vantage system outputs an external trigger signal, which drives a function generator to produce an impulsive tap, driving a mechanical piston attached to a surface transducer (voice coil). When placed superficially to the tissue of interest, the wave is induced in the tissue and then measured using the ultrasound probe.

For processing, it can be completed on the Verasonics-supplied computer, or on a different computer that is running MATLAB. This version has only been tested for R2022a. I highly recommend that you have sufficient RAM and storage space as ultrasound-tensiometry files can be quite large and time-consuming to process. 


---Collecting ultrasound-tensiometry---
On the Verasonics computer, ensure the system is running. Open the SetUpCL15_7v_ultrasound_tensiometry.m code in MATLAB. Ensure the working directory is the Vantage folder. The first section contains the parameters to change. Based on the tissue of interest and the depth, consider changing the start depth, end depth and the desired resolution. 

Note: as you collect more "depths", the data file requirements can increase. The code will error out if you request more than 200 pixels of depth - i.e. calculated by (endDepthMM - startDepthMM) / resDesiredMM  gives the number of "depth" pixels. If exceeded, it will create unmanageable structure sizes and generally causes MATLAB to crash.

Change the number of taps (numTaps) to align with your needs (ensure it is an even number). The tap rate can change based on your collection requirements. If you exceed the possible tap rate, when you run a collection through the GUI, you will receive time out errors from Verasonics. Currently, it is a bit time-consuming to determine the maximum allowable capture rate. Generally, as you collect less data per tap, you can capture at a higher rate.