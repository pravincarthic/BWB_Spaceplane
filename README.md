# Numerical Investigation of Discrete Symmetry based Passive Micro-cavity Patterns for Boundary-Layer Control in a Hypersonic Blended Wing Body Spaceplane

## Abstract

This repository contains the computational fluid dynamics (CFD) setup, geometric definition files, solver configurations, and mathematical frameworks for a novel hypersonic spaceplane architecture. The vehicle features a unique design that integrates the **Blended Wing Body (BWB)** configuration with the **flying wing concept**, combining the volumetric efficiency and smooth aerodynamic blending of a BWB with the pure, tailless lifting planform efficiency of a flying wing.

In this architecture, there is no conventional distinction between the aft fuselage and outer wings; the entire aft section functions as a single, continuous lifting wing. This hybrid outer mold line (OML) is augmented with an upstream array of **passive micro-cavities (acoustic metasurfaces)** to damp acoustic instabilities, delay laminar-to-turbulent boundary layer transition, and enhance aerodynamic efficiency during post-reentry atmospheric glide (Mach 8 to Mach 10). Flow modeling is performed as **implicit LES (ILES)** in stock **OpenFOAM v2412** using **`rhoCentralFoam`**, with a calorically imperfect frozen-composition air model (`janaf` thermo and `polynomial` transport, both valid to 6000 K) over a 30.375 ms physical time horizon.

---

## Key Physical Concepts & Innovation

### 1. Integrated BWB and Flying Wing Reentry Geometry
Rather than relying on a conventional fuselage or a standard tube-and-wing layout, this architecture uses a unique hybrid design integrating the Blended Wing Body and flying wing concepts:
* **All-Lifting Aft Structure:** Eliminates traditional tail surfaces, discrete wing-body joints, and non-lifting aft body sections. The entire aft portion of the vehicle is contoured as a unified, continuous lifting surface that generates uniform spanwise lift.
* **Seamless Volume and Lift Integration:** Merges the thick, volume-efficient centerbody of a BWB with the continuous lifting planform of a flying wing, maximizing usable internal volume for cryogenic propellants while eliminating parasitic non-lifting structures.
* **Wetted Area Reduction:** Reduces global wetted surface area by up to 33% relative to equivalent payload volume tube-and-wing designs, directly decreasing hypersonic skin-friction drag.
* **Compression Lift:** Eliminates distinct wing-body junction corners, creating a continuous compression ramp on the windward surface that captures the high-pressure bow shock wave evenly.
* **Spanwise Load Distribution:** Distributes aerodynamic lift across the wide centerbody (carrying 30 to 45% of total lift), alleviating root bending moments and reducing structural weight.
* **Target Aspect Ratio:** Optimized at approximately 1.275:1 (Length to Span), placing the vehicle in a high-efficiency hypersonic glide envelope while maintaining transonic and subsonic static pitch stability.

### 2. Acoustic Metasurfaces for Boundary Layer Stabilization
In hypersonic boundary layers ($M \ge 4$), the primary mechanism triggering turbulent transition is the acoustic **Mack second-mode instability**.
* **Resonant Energy Absorption:** Subwavelength micro-cavity arrays act as passive acoustic absorbers, trapping and dissipating high-frequency acoustic wave energy within the boundary layer.
* **Acoustic Impedance & Phase Matching:** The sub-surface micro-pockets modify the wall acoustic admittance phase angle, forcing destructive interference on reflecting pressure waves trapped between the wall and the relative sonic line.
* **Vortex Trapping & Fluidic Damping:** Fluid shear layers spanning each micro-cavity drive stable, localized recirculation zones that absorb turbulent fluctuations and prevent cross-flow vortex amplification along the forebody transition zone.
* **Strategic Placement Strategy:** Metasurface patches are positioned strictly upstream in the linearly unstable second-mode growth region. To prevent skin friction penalties induced by shear-dissipation overshoots, the metasurface terminates prior to the non-linear breakdown zone.

---

## Aerodynamic & Aerothermal Design Specifications

| Parameter | Baseline Value / Target | Physical Significance |
| :--- | :--- | :--- |
| **Freestream Regime** | Mach 8.0 to 10.0 (Post-Reentry Glide) | Severe hypersonic continuum flow with localized air ionization |
| **Physical Simulation Horizon** | 30 ms | Sufficient time for flow-field development and acoustic wave passage |
| **Aft Surface Configuration** | Fully Integrated All-Lifting Wing | Converts entire rear airframe into a continuous lifting body |
| **Aspect Ratio (L:S)** | 1.275 : 1 | Optimal trade-off between volume, cross-range, and stability |
| **Micro-Cavity Scale** | Sub-millimeter to millimeter ($5\text{ mm} \times 10\text{ mm}$) | Scaled to match acoustic second-mode wavelengths ($\lambda_2 \approx 2\delta$) |
| **Gas Model** | 11-species reacting air ($N_2, O_2, NO, N, O, N_2^+, O_2^+, NO^+, N^+, O^+, e^-$) | Thermo-chemical non-equilibrium (multi-temperature $T, T_{vib}, T_e$) |

---

## Computational Fluid Dynamics (CFD) Architecture: OpenFOAM v2412

All numerical simulations are performed using stock **OpenCFD OpenFOAM v2412**
with **`rhoCentralFoam`**, the density-based central-upwind solver. No custom
solver or shared library is required.

The configuration previously targeted the hyStrath suite (`hy2Foam`, 11-species
reacting air, two-temperature non-equilibrium). It was replaced because the
peak temperature this vehicle actually reaches is 4331 K, below the point where
dissociation materially changes the aerodynamics of interest, while the caloric
imperfection and the high-temperature viscosity - which do matter - are
captured by `janaf` and a fitted `polynomial` transport model. The reasoning is
set out in `scripts/OpenFOAM configs/hypersonic_case/docs/SETTINGS.md`.

### Solver Configuration Summary
* **Simulation Framework:** OpenFOAM v2412 (`rhoCentralFoam`), 384 MPI ranks, collated I/O.
* **Fidelity Level:** **Implicit LES (ILES)** on a wall-resolved (y+ below 1) hybrid mesh. `laminar` / `Stokes` stress model: the vanLeer limiters supply the subgrid dissipation, with no explicit SGS model.
* **Governing Equations:** 3D compressible Navier-Stokes, single static temperature, frozen-composition air. `sensibleInternalEnergy`, `perfectGas`, `janaf` cp(T) and `polynomial` mu(T)/kappa(T) valid 200-6000 K.
* **Convective Flux Discretization:** `Kurganov` (Kurganov-Noelle-Petrova central-upwind) with `vanLeer` / `vanLeerV` TVD limiters on every reconstruction.
* **Temporal Integration:** Explicit multi-stage Runge-Kutta (RK3) with adaptive time-stepping bound by $Co \le 0.4$ to resolve micro-cavity acoustic frequencies ($\Delta t \approx 10^{-9}\text{ s}$).
* **Boundary Layer Mesh Quality:** Unstructured/hybrid prism-hex grids targeting $y^+ < 1$ with a growth rate $\le 1.12$ near wall surfaces and cavity apertures.
