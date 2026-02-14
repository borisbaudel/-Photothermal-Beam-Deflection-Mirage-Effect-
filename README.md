# 🔬 Photothermal Beam Deflection (Mirage Effect)

![Physics](https://img.shields.io/badge/Field-Experimental%20Physics-blue)
![Optics](https://img.shields.io/badge/Topic-Laser%20Optics-red)
![Signal Processing](https://img.shields.io/badge/Method-Signal%20Processing-green)
![Detection](https://img.shields.io/badge/Technique-Lock--in%20Detection-orange)

---

## 📖 Overview

This project presents an **experimental and theoretical investigation of laser beam deflection**
induced by:

- 🌡️ **Photothermal gradients** (Mirage Effect)
- 🌀 **Acousto-optic interactions**

The study combines **precision optics**, **thermal-wave physics**, and **advanced signal processing**
to characterize beam displacement, sensitivity limits, and noise behavior.

---

## 🌡️ Physical Principle — Mirage Effect

A modulated heat source generates a **temperature gradient in air**, producing a spatial variation
of the refractive index:

$$ n(T) \quad \Rightarrow \quad \frac{\partial n}{\partial T} \neq 0 $$

This gradient deflects the laser beam:

$$ \theta = \frac{1}{n} \frac{\partial n}{\partial T}
\int \frac{\partial T}{\partial y} \, dy $$ 

where:

- \( \theta \) : beam deflection angle  
- \( n \) : refractive index  
- \( T \) : temperature  

---

## 📐 Measurement Chain

The deflection is converted through a **three-stage transduction process**:

### 1️⃣ Angular → Spatial Shift

$$ \Delta y = L \cdot \Delta \theta $$

- \( L \) : propagation distance  

---

### 2️⃣ Spatial Shift → Photodiode Signal

Beam displacement modifies the **intensity imbalance** on a split detector.

---

### 3️⃣ Photodiode Current → Voltage

Signal amplification via **transimpedance amplifier**.

---

## ⚙️ Experimental Setup

Key components:

- 🔦 Probe laser (Gaussian beam)
- 🌡️ Modulated Peltier thermal source
- 🪞 Optical alignment stage
- 📡 Split / quadrant photodiode
- 🎚️ Lock-in amplifier (phase-sensitive detection)

---

## 📡 Detection Techniques

### 🎯 Lock-in Detection

Extracts the signal at the modulation frequency:

✔ Amplitude  
✔ Phase  
✔ Noise rejection  

---

### 📊 FFT Analysis

Noise spectral decomposition:

- Low-frequency drift (1/f)
- 50 Hz electrical noise
- Mechanical vibration peaks
- Flat electronic noise floor

---

## 🌀 Thermal Wave Modeling

Temperature field under harmonic excitation:

$$ T(y,t) =
T_0 e^{-y/L_a}
\cos\left(2\pi f t - \frac{y}{\lambda_T} + \phi \right) $$

Measured quantities:

- Thermal attenuation length \( L_a \)
- Thermal wavelength \( \lambda_T \)

---

## 🧠 Signal & Data Processing

Implemented methods:

- Gaussian filtering  
- Cubic interpolation  
- Error-function fitting  
- Linear regression  
- Spectral analysis (FFT)  

---

## 🎯 Key Results

✔ Microradian-scale beam deflection measurement  
✔ Experimental estimation of thermal attenuation length  
✔ Phase evolution vs. height  
✔ Sensitivity & noise floor characterization  
✔ Validation of synchronous detection  

---

## 📚 Concepts Covered

- Mirage effect / photothermal deflection  
- Refractive index gradients  
- Gaussian beam optics  
- Knife-edge error-function model  
- Lock-in detection  
- Noise propagation  
- Sensitivity optimization  

---

## 📚 References

- D. J. Griffiths — *Introduction to Electrodynamics*  
- Saleh & Teich — *Fundamentals of Photonics*  
- Mandel & Wolf — *Optical Coherence and Quantum Optics*  
- Stanford Research Systems — Lock-in Amplifier Application Notes  
- Gonzalez & Woods — *Digital Image Processing*

---

## 🚀 Relevance & Applications

This work is relevant to:

- Precision optical metrology  
- Quantum sensing & magnetometry  
- Photothermal spectroscopy  
- Beam diagnostics  
- Low-signal detection systems  

---

## 👨‍🔬 Authors

**Boris Baudel**  
Le Mans Université  

Paul Barraud  
Le Mans Université  

Etienne Raguillat  
Le Mans Université  

---

## ✨ Repository Purpose

✔ Laboratory project documentation  
✔ Experimental physics portfolio  
✔ Signal-processing demonstration  
✔ Optical sensing study  

---
