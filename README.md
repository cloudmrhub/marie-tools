# **MARIE 3.0 — Magnetic Resonance Integral Equation Suite**

**MARIE 3.0** (Magnetic Resonance Integral Equation Suite) is an open-source MATLAB toolbox for the **electromagnetic simulation** of radiofrequency (RF) coils in Magnetic Resonance Imaging (MRI).

---

## 🧠 **Overview**
**MARIE 3.0** extends the functionalities of **MARIE** and **MARIE 2.0** by introducing 
* A **wire-surface-volume integral equation (WSVIE)** coupled solver that allows modeling of any type of coil array in MRI
* **Tensor decomposition-based compression** of the WSVIE Green's function operators leading to significantly reduced memory footprint and acceleration through GPU programming
* A **fully automated circuit co-simulation pipeline** that mimics the RF lab coil prototyping process
* The **Magnetic Resonance Green's Function (MRGF)**, a reduced-order model of the WSVIE for ultra-fast, patient-specific simulations
* Routines to compute the **Ultimate Intrinsic SNR (UISNR)**, **Ideal Current Pattern**, and **Ultimate Intrinsic Efficiency** computations.

---

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
| `TMD_runner.m`   | Performs **automated tuning, matching, and decoupling** optimization of, based on Y parameters generated with external solvers (co-simulation is also performed for MARIE and MRGF runners by default)|
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
4. From MATLAB to install all c++ accelerated assembly routines
   ```
   startup.m
   installer.m
   ```
5. Run your first simulation:

   ```matlab
   MARIE_runner
   ```

   Input parameters are defined in `/data/inputs/`.

---

## ⚠️ **Disclaimer**
This is an in-house electromagnetic simulator suire and it is intended for research purposes. While every effort has been made to ensure its quality, the code may contain bugs. We do not assume responsibility for any errors or issues resulting from the software's use. We strongly encourage the users to post any errors, warnings, questions, etc, in the Issues page, or send directly to the corresponding authors.

---

## 📚 **References**

If you use MARIE 3.0, please cite:

> Giannakopoulos, I. I., Zhang, B., Serrallés Cruz, J. E., Brown, R., Lattanzi, R.
> *An Open-Source Software Toolbox for Rapid Radiofrequency Coil Design and Evaluation in MRI*,
> Magnetic Resonance in Medicine, 2025.

### MARIE 3.0 inherited a significant amount of routines from the following repositories

1. Geuzaine, Christophe, and Jean‐François Remacle. "Gmsh: A 3‐D finite element mesh generator with built‐in pre‐and post‐processing facilities." International journal for numerical methods in engineering 79.11 (2009): 1309-1331.
2. Oseledets, I. *TT-Toolbox: Tensor Train Decomposition Library.* [https://github.com/oseledets/TT-Toolbox](https://github.com/oseledets/TT-Toolbox)
3. Polimeridis, A. **MARIE** [https://github.com/thanospol/MARIE](https://github.com/thanospol/MARIE)
4. Guryev, G. **MARIE 2.0** [https://github.com/georgyguryev/MARIE_2.0](https://github.com/georgyguryev/MARIE_2.0)
5. Polimeridis, A. **DIRECTFN** [https://github.com/thanospol/DIRECTFN](https://github.com/thanospol/DIRECTFN)
6. Polimeridis, A. **DEMCEM** [https://github.com/thanospol/DEMCEM](https://github.com/thanospol/DEMCEM)

---

### Papers on MARIE (please cite explicitly if needed)
1. **Piecewise-constant approximation of VIE**: Polimeridis, Athanasios G., et al. "Stable FFT-JVIE solvers for fast analysis of highly inhomogeneous dielectric objects." Journal of Computational Physics 269 (2014): 280-296.
2. **Power Calculations**: Polimeridis, Athanasios G., et al. "On the computation of power in volume integral equation formulations." IEEE Transactions on Antennas and Propagation 63.2 (2014): 611-620.
2. **MARIE**: Villena, Jorge Fernandez, et al. "MARIE–a MATLAB-based open source software for the fast electromagnetic analysis of MRI systems." Proceedings of the 23rd Annual Meeting of ISMRM, Toronto, Canada. 2015.
3. **Magnetic Resonance Green's function**: Villena, Jorge Fernández, et al. "Fast electromagnetic analysis of MRI transmit RF coils based on accelerated integral equation methods." IEEE Transactions on Biomedical Engineering 63.11 (2016): 2250-2261.
4. **Evaluation of Signular Integrals**: Tambova, Alexandra A., et al. "On the generalization of DIRECTFN for singular integrals over quadrilateral patches." IEEE Transactions on Antennas and Propagation 66.1 (2017): 304-314.
5. **Tucker decomposition-based compresion of VIE**: Giannakopoulos, Ilias I., Mikhail S. Litsarev, and Athanasios G. Polimeridis. "Memory footprint reduction for the FFT-based volume integral equation method via tensor decompositions." IEEE Transactions on Antennas and Propagation 67.12 (2019): 7476-7486.
6. **Precorrected Fast Fourier Transform-based SVIE**: Guryev, Georgy D., et al. "Fast field analysis for complex coils and metal implants in MARIE 2.0." Proc. ISMRM. 2019.
7. **Piecewise-linear approximation of VIE**: Georgakis, Ioannis P., et al. "A fast volume integral equation solver with linear basis functions for the accurate computation of EM fields in MRI." IEEE Transactions on Antennas and Propagation 69.7 (2020): 4020-4032.
8. **Ultimate Transmit Efficiency**: Georgakis, Ioannis P., Athanasios G. Polimeridis, and Riccardo Lattanzi. "A formalism to investigate the optimal transmit efficiency in radiofrequency shimming." NMR in Biomedicine 33.11 (2020): e4383.
9. **Tucker decomposition-based compression of SVIE**: Giannakopoulos, Ilias I., et al. "Compression of volume-surface integral equation matrices via Tucker decomposition for magnetic resonance applications." IEEE transactions on antennas and propagation 70.1 (2021): 459-471.
10. **Tensor Train-based SVIE**: Giannakopoulos, Ilias I., et al. "A tensor train compression scheme for remote volume-surface integral equation interactions." 2021 International Applied Computational Electromagnetics Society Symposium (ACES). IEEE, 2021.
11. **Ultimate Electromagnetic Basis and Ultimate SNR**: Georgakis, Ioannis P., et al. "Novel numerical basis sets for electromagnetic field expansion in arbitrary inhomogeneous objects." IEEE transactions on antennas and propagation 70.9 (2022): 8227-8241.
12. **Hybrid cross-Tensor Train and Precorrected Fast Fourier Transform-based SVIE**: Giannakopoulos, Ilias I., et al. "A hybrid volume-surface integral equation method for rapid electromagnetic simulations in MRI." IEEE Transactions on Biomedical Engineering 70.1 (2022): 105-114.
13. **Perturbation Matrix-based SVIE**: Guryev, Georgy D., et al. "MARIE 2.0: a perturbation matrix based patient-specific MRI field simulator." IEEE Transactions on Biomedical Engineering 70.5 (2022): 1575-1586.
14. **Ideal Current Patterns**: Giannakopoulos, Ilias I., et al. "Computational methods for the estimation of ideal current patterns in realistic human models." Magnetic resonance in medicine 91.2 (2024): 760-772.

## 💬 **Acknowledgments**

This work is currently supported in part by **NIH K99 EB035163**, **R01 EB036483**, and **R01 EB024536**, and performed under the **Center for Advanced Imaging Innovation and Research (CAI²R)** — an NIH P41 National Center for Biomedical Imaging and Bioengineering (P41 EB017183).

---

## 🧑‍🔬 **Contact**

**Ilias I. Giannakopoulos, PhD**
Center for Biomedical Imaging, NYU Grossman School of Medicine
📧 [ilias.giannakopoulos@nyulangone.org](mailto:ilias.giannakopoulos@nyulangone.org)
🌐 [www.cai2r.net](https://www.cai2r.net)

---

## 🧑‍💻 **Contibutors**

**Ilias I. Giannakopoulos, PhD**, Center for Biomedical Imaging, NYU Grossman School of Medicine

**José E. Cruz Serrallés, PhD**, Center for Biomedical Imaging, NYU Grossman School of Medicine

**Georgy D. Guryev, PhD**, Industry

**Ioannis P. Georgakis, PhD**, Industry

**Alexandra Tambova, PhD**, Industry

**Mikhail S. Litsarev, PhD**, Industry

**Jorge F. Villena, PhD**, Industry

**Luca Daniel, PhD**, Electrical Engineering and Computer Science, Massachusetts Institute of Technology

**Jacob K. White, PhD**, Electrical Engineering and Computer Science, Massachusetts Institute of Technology

**Athanasios G. Polimeridis, PhD**, Industry

**Riccardo Lattanzi, PhD**, Center for Biomedical Imaging, NYU Grossman School of Medicine

---

## 📃 **License**

MARIE 3.0 is released under the **MIT License**.
Please cite the above publication when using this software in academic or commercial work.
