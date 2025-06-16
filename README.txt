# Ultrasound-tensiometry

This GitHub repository contains the MATLAB code that enables users to capture **ultrasound-tensiometry** data
(see manuscript: [https://doi.org/10.1016/j.gaitpost.2024.06.026](https://doi.org/10.1016/j.gaitpost.2024.06.026)).

Ultrasound-tensiometry is used to assess the spatial variation in shear wave speed in response to a mechanical tap.
This repository contains two main sections:

* `Collection/`
* `Processing/`

## Table of Contents

* [Installation](#installation)
* [Required Hardware](#required-hardware)
* [Required MATLAB Toolboxes](#required-matlab-toolboxes)
* [Collecting Ultrasound-Tensiometry Data](#collecting-ultrasound-tensiometry-data)
* [Processing Notes](#processing-notes)

---

## Installation

Clone the repository directly into the `Vantage` folder on the Verasonics computer to enable data collection.
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

1. On the Verasonics computer, make sure the system is powered on.
2. Open `SetUpCL15_7v_ultrasound_tensiometry.m` in MATLAB.
3. Ensure the **working directory is the `Vantage` folder**.

### Parameters to Update

Update parameters in the first section of the script based on:

* Target tissue
* Start and end depth (in mm)
* Desired resolution (in mm)

> *Limit depth to avoid excessive memory use:*
> The number of depth pixels is calculated as:
> `(endDepthMM - startDepthMM) / resDesiredMM`
> If this exceeds **200 pixels**, MATLAB will likely crash due to oversized data structures.

### Taps and Capture Rate

* Set `numTaps` (must be an **even** number).
* Tap rate can be adjusted, but exceeding the allowable rate will cause **timeout errors** from Verasonics.

> Tip: Lower data volume per tap → higher allowable tap rate.
> Unfortunately, finding the ideal rate is currently trial-and-error.

---

## Processing Notes

* Processing can be done on the Verasonics computer **or** any computer with MATLAB installed.
* This version has only been tested with **MATLAB R2022a**.
* **Large files:** Ensure sufficient **RAM** and **disk space**—ultrasound-tensiometry data is heavy.

