# **MARIE — Backend Tools**

Visit the [project page](https://cmr.cloudmrhub.com/apps/marie/) for more information about **MARIE**.

## ⚙️ **Repository Structure**

### **1. 📇 `src/` — Core Source Code**

| Subfolder                 | Description                                                                                                                                         |
| ------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------- |
| `src_geometry/`           | Parsers for GMSH surface and wire coil formats|
| `src_integral_equations/` | Implements all integral formulations: **WIE**, **SIE**, **VIE**, **WVIE**, **SVIE**, and **WSVIE**|
| `src_solver/`             | Iterative solvers (**GMRES**) with support for **block-Jacobi** and **LU preconditioning**, GPU acceleration.|
| `src_mathematics/`        | Tensor decompositions (**SVD**, **Tucker**, **Tensor Train**) and **Adaptive Cross Approximation (ACA)** for matrix compression.|
| `src_physics/`            | Computation of electromagnetic fields and SNR. Electronic simulator.|
| `src_utils/`              | General utilities for file I/O, `.json` parsing, and data management.|

---

### **2. 🔌 `src_runners/` — Simulation Entry Points**

Each runner is a self-contained script that performs geometry assembly, matrix construction, system solution, and post-processing.

| Runner           | Description                                                                                 |
| ---------------- | ------------------------------------------------------------------------------------------- |
| `MARIE_runner.m` | Executes the full **WSVIE** simulation.                                                     |
| `BASIS_runner.m` | Generates an **elcectromangetic basis**, **UISNR** or the **ideal current patterns**, and the **MRGF**|
| `MRGF_runner.m`  | Runs the **MRGF-based reduced-order WSVIE** for faster simulations                 |
| `TMD_runner.m`   | Performs **automated tuning, matching, decoupling, preamplifier decoupling, and detuning** optimization of receive, transmit, and transceive coils (and their combinations) based on pre-computed Y parameters with external solvers (co-simulation is also performed for MARIE and MRGF runners by default)|
---

### **3. 🗂 `data/` — Simulation Data**

```
data/
├── bases/         # EM Bases files (include UISNR and MRGF)
├── bodies/        # Body voxel models stored in mat files
├── coils/         # Coil meshes and lumped element files
├── inputs/        # Input simulation files
└── solutions/     # Simulation results
```

---

### **4. 📒 `doc/` — Documentation**
(Coming Soon...)
```
doc/
```

---

## 🚀 **Getting Started**

1. Install MATLAB R2023a or later.
2. **GMSH** (2.85) is provided with the repository (2.8.5).
3. Clone the repository and initialize:

   ```bash
   git clone https://github.com/cloudmrhub/marie-tools
   cd MARIE
   ```
4. From MATLAB to install all c++ routines
   ```
   installer.m
   startup.m
   ```
5. Run simulations from the home directory:

   ```matlab
   MARIE_runner
   ```

   Input parameters are defined in `/data/inputs/`.

---

## 💬 **Acknowledgments**

This work is supported in part by **NIH K99 EB035163**, **NIH R01 EB036483**, and **NIH R01 EB024536**, and performed under the **Center for Advanced Imaging Innovation and Research (CAI²R)**, an NIH P41 National Center for Biomedical Imaging and Bioengineering (P41 EB017183).

---

## ⚠️ **Disclaimer**
This is an in-house electromagnetic simulator suite and it is intended for research purposes. While every effort has been made to ensure its quality, the code may contain bugs. We do not assume responsibility for any errors or issues resulting from the software's use. We strongly encourage the users to post any errors, warnings, questions, etc, in the Issues page, or send directly to the corresponding author (ilias.giannakopoulos@nyulangone.org).

---

## 📃 **License**

MARIE is released under the **MIT License**. 
Please cite the following publication when using this software in academic or commercial work:
Giannakopoulos, Ilias I., et al. "An Open‐Source Software Toolbox for Rapid Radiofrequency Coil Design and Evaluation in MRI." Magnetic Resonance in Medicine (2026)
