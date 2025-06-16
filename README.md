# Ultrasound-tensiometry

This GitHub repository contains the MATLAB code that enables users to capture **ultrasound-tensiometry** data
(see manuscript: [https://doi.org/10.1016/j.gaitpost.2024.06.026](https://doi.org/10.1016/j.gaitpost.2024.06.026)).

Ultrasound-tensiometry is used to assess the spatial variation in shear wave speed in response to a mechanical tap.
This repository contains two main sections:

* `CollectData/`
* `ProcessData/`

## Table of Contents

* [Installation](#installation)
* [Required Hardware](#required-hardware)
* [Required MATLAB Toolboxes](#required-matlab-toolboxes)
* [Collecting Ultrasound-Tensiometry Data](#collecting-ultrasound-tensiometry-data)
* [Processing Notes](#processing-notes)

---

## Installation

Clone the repository directly into the `Vantage #.#.#` folder on the Verasonics computer to enable data collection.
The processing code **does not** rely on Verasonics functions and can be run on any other computer.

---

## Required Hardware

* Verasonics Vantage 64LE Ultrasound System
* Philips CL15-7 Ultrasound Probe

> *Other probes and Vantage systems (e.g., 128, 256) may be supported in the future. Development of a 128 implementation is in progress.*

To generate the mechanical tap:

* The Vantage system outputs an **external trigger signal**
* This drives a **function generator**, which sends an impulsive signal to a **voice coil actuator**
* The actuator (a mechanical piston) delivers a tap to the surface of the tissue
* The ultrasound probe measures the resulting shear wave in the tissue

---

## Required MATLAB Toolboxes

| Toolbox                                 | Required?   |
| --------------------------------------- | ----------- |
| Signal Processing Toolbox               |   Yes       |
| Image Processing Toolbox                |   Yes       |
| Statistics and Machine Learning Toolbox |   Yes       |
| Parallel Computing Toolbox              |  Optional   |

---

## Collecting Ultrasound-Tensiometry Data

This code is originally based on `RunSetUpL11_5vFlashHFR_acquire.m`, the **"High Frame Rate" acquisition** code for ultrafast imaging with a Vantage 64LE system, written for **Vantage 4.4.0**. It has been modified by **Darryl Thelen** and **Lauren Welte**.

1. On the Verasonics computer, make sure the system is powered on.
2. Open `SetUpCL15_7v_ultrasound_tensiometry.m` in MATLAB.
3. Ensure the **working directory is the `Vantage #.#.#` folder**.
4. Update the [parameters](#-Parameters-to-Update) as desired
5. Click **Run** and it will create the set up file and run VSX automatically. 
---
### Parameters to Update
#### Imaging Parameters

Update parameters in the first section of the script based on:

* Target tissue
* Start and end depth (in mm)
* Desired resolution (in mm)

> *Limit the number of pixels of depth to avoid excessive memory use:*
> The number of depth pixels is calculated as:
> `(endDepthMM - startDepthMM) / resDesiredMM`
> If this exceeds **200 pixels**, MATLAB will likely crash due to oversized data structures.

#### Taps and Capture Rate

* Set the number of desired tap events in `numTaps`.
* Tap rate (frequency of tapping) can be adjusted, but exceeding the allowable rate will cause **timeout errors** from Verasonics.

> Tip: Lower data volume per tap → higher allowable tap rate.
> Unfortunately, finding the ideal rate is currently trial-and-error.

---

### Sequence Overview


* The `Event` structure controls the sequence of imaging events.
* Two primary modes:

  1. **Initial Plane Wave Imaging** for probe placement and image guidance (loops until `Start` is pressed).
  2. **Tensiometry Mode**, triggered after pressing `Start`, where:

     * A **trigger out** signal initiates an external mechanical tap.
     * The system captures `tapSamples` of plane wave data at the defined `sampleRate` (recommend 20,000 Hz).
     * This repeats for `numTaps`, using the specified `tapRate`.

>  MATLAB may **freeze** during data acquisition due to large data throughput. This is expected behaviour.
---
### After Acquisition

* The system prompts the user to save the data.
* Saving can take **several minutes** depending on depth and resolution.
* The data are stored in a `.mat` file, which includes:

  * `RF`: first frame from each tap (for reconstruction)
  * `IQ`: full imaging data
  * Parameters: `TX`, `TW`, `Trans`, `P`, `PData`, `CollectionParams`, `Resource`, `Receive`, `SeqControl`

---
### Tap Check

Use the **Tap Check GUI** to check the data quality by clicking `Tap Check`:

1. Select number of taps to review (e.g., 10 out of 20 will show taps 1, 3, 5...19).
2. Click `Select ROI`, and draw a region from superficial to deep on the ultrasound image.
3. Click `Calculate` to analyze wave propagation.
4. Navigate using:

   * **Arrow keys** (←/→): switch between taps
   * **Arrow keys** (↑/↓): switch between depths
   * Or use the on-screen buttons
   > If the arrow keys are not working, click in the background of the GUI, away from any axes or element and try again. 

This provides a quick qualitative assessment of wave propagation before full processing.


---

## Processing Notes

* Processing can be done on the Verasonics computer **or** any computer with MATLAB installed.
* This version has only been tested with **MATLAB R2022a**.
* **Large files:** Ensure sufficient **RAM** and **disk space**—ultrasound-tensiometry data is heavy.

1. Start by opening `CalculateDisplacement` by typing in the command window.
2. Select the folder with your tensiometry files.
3. Click or CTRL + Click the files you want to process.
4. Choose parameters for the Loupas algorithm under the "Calculate Displacement" checkbox. Generally you won't need to change these.
5. To calculate the displacement peaks along the depth lines, choose "Radon Transform". It is much more stable. 
6. Choose the number of taps to process. It will automatically choose "All Taps". If you want less, uncheck the box and fill in the text boxes with the desired tap range.
7. If you want a region of interest (i.e. select only an area of the image for processing), check "Select ROI". This will pop up with an image for you to select a region.
8. If you have the parallel computing toolbox installed, you will have the option to use it. It is generally much faster to use the PCT.
9. Click calculate. For the options selected, it will save out a `filename_displacement.mat` and a `filename_radon.mat` file. 
    * The displacement file contains two matrices, one with the raw displacement and one with filtered displacement with the indices: (depth, width, withintapImages,tapNumber).
    * The radon file contains  :
        * `SWS` which is a cell array with SWS{tap number}, and each one contains the shear wave speed at each depth pixel in the region of interest (or the entire image if selected).
        * `Rsq` this is a measure of the sum of the Radon transform. It gives a sense of how reliable the wave speed measure is. This can be adjusted in the visualization GUI to omit low confidence values
        * `method` will give "Radon Transform" as a string.
        * `ROI_dx`, `ROI_wx` give the values in mm of the region of interest pixels in depth (dx) and width (along transducer) (wx)

If you want to visualise the data, open `visualiseHFR_Radon` by typing it into the command window. Click `File` in the menu bar and `Load`. Choose the ORIGINAL file. It will automatically look for the _displacement and _radon files. 
