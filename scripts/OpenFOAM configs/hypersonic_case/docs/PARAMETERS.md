# Parameters

## Units

**The OpenFOAM case is entirely in SI: metres, seconds, kilograms, kelvin,
pascals.** Every coordinate below - probe stations, rake lengths, cutting
plane bounds, `CofR`, `lRef`, `Aref` - is in metres.

The mesh is not. `scripts/final_meshing.py` scales the STEP geometry by 1000
so that the sizing constants in `scripts/gmsh_config.json` (`size_near` 0.5,
`dist_max` 25000, `far_field_length_mm` 200000) can be written in
millimetres, so the `.msh` it writes is in **millimetres**.

`gmshToFoam` does not rescale on its own. The conversion must therefore be
run as `gmshToFoam <mesh> -scale 0.001`, which `scripts/setup.sh` does by
default (`MESH_SCALE`). Two guards exist because getting this wrong is silent
- the solver runs happily on a 48 km vehicle:

- `scripts/setup.sh` step 3a parses the `checkMesh` bounding box and aborts if
  the largest coordinate exceeds 1000 (metres would give about 100)
- `scripts/make_pse_probes.py` measures the wall surface and refuses to write
  probes if the vehicle is not about 48 m long

| Reads as | In metres | In the mm mesh |
|---|---|---|
| Vehicle length | 48.061 | 48061 |
| Far-field box | 200 | 200000 |
| Minimum cell | 5.0e-4 | 0.5 |
| First rake point | 2.0e-5 | 0.02 |

Every number this case runs on, with the working behind it. Cross-check
against the dictionaries; if the two disagree, the dictionary is what runs.

## Solver

| Item | Value | Where |
|---|---|---|
| Solver | `rhoCentralFoam` | `system/controlDict` |
| OpenFOAM | v2412 (stock OpenCFD) | - |
| Extra libraries | none | - |

Density-based, explicit, central-upwind. No hyStrath/hy2Foam component is
used and none of the `libstrath*.so` libraries are loaded.

## Freestream (130,000 ft, 39.6 km)

| Quantity | Value |
|---|---|
| Static pressure `p_inf` | 292.12 Pa |
| Static temperature `T_inf` | 250 K |
| Velocity `u_inf` | 3169.7 m/s, axial |
| Gas constant `R` | 287.102 J/kg/K (molWeight 28.95999) |
| Density `rho_inf` | 0.00406992 kg/m3 |
| `gamma` at 250 K | 1.4035 (janaf, not 1.4) |
| Speed of sound `a_inf` | 317.39 m/s |
| Mach `M_inf` | 9.987 |
| Viscosity `mu_inf` | 1.6014e-05 kg/m/s |
| Unit Reynolds number | 8.056e5 per metre |
| `Re_L` (L = 48.061 m) | 3.872e7 |

`M_inf` is 9.987 rather than 10.000 because the calorically imperfect
`gamma` at 250 K is 1.4035, not 1.4. `u_inf` is held at the specified
3169.7 m/s rather than being adjusted to force `M = 10` exactly.

Angle of attack is -5 degrees and is baked into the meshed geometry
(the vehicle is rotated in Gmsh). The freestream vector is therefore purely
axial and the lift/drag axes stay aligned with the domain axes.

## Thermophysical model

| Item | Value |
|---|---|
| `type` | `hePsiThermo` (required by `rhoCentralFoam`) |
| `mixture` | `pureMixture` - one frozen-composition air pseudo-species |
| `thermo` | `janaf`, Tlow 200 K, **Thigh 6000 K**, Tcommon 1000 K |
| `transport` | `polynomial`, degree 7, valid 200-6000 K |
| `equationOfState` | `perfectGas` |
| `energy` | **`sensibleInternalEnergy`** |
| `molWeight` | 28.95999 |

Composition used to build the coefficients: X(N2) 0.78157, X(O2) 0.20964,
X(Ar) 0.00935 (dry air, CO2 dropped and renormalised).

The janaf coefficients are the mole-fraction-weighted sum of the Burcat N2,
O2 and Ar polynomials, which is exact for a molar NASA-7 fit of an ideal
mixture at fixed composition.

| T [K] | cp [J/kg/K] | gamma | mu [kg/m/s] | kappa [W/m/K] | Pr |
|---|---|---|---|---|---|
| 250 | 998.7 | 1.4035 | 1.601e-05 | 2.315e-02 | 0.691 |
| 300 | 1003.6 | 1.4007 | 1.842e-05 | 2.671e-02 | 0.692 |
| 1000 | 1140.6 | 1.3364 | 4.149e-05 | 6.759e-02 | 0.700 |
| 2000 | 1251.9 | 1.2976 | 6.611e-05 | 1.179e-01 | 0.702 |
| 4000 | 1319.5 | 1.2781 | 1.090e-04 | 2.041e-01 | 0.705 |
| 6000 | 1357.2 | 1.2683 | 1.492e-04 | 2.861e-01 | 0.708 |

cp is continuous across Tcommon to 0.19 percent.

### Why Thigh 6000 K is the right ceiling

Normal shock at these conditions, solved with the thermally perfect janaf
gas at frozen composition:

| Quantity | Value | Ratio |
|---|---|---|
| `T2` | 4256 K | 17.0 x `T_inf` |
| `p2` | 35,445 Pa | 121.3 x `p_inf` |
| `rho2` | 0.02901 kg/m3 | 7.13 x `rho_inf` |
| `u2` | 444.7 m/s | - |
| Stagnation `T0` | 4331 K | - |

The whole thermal range of the run is 250 K to 4331 K. Thigh 6000 K is
therefore a genuine 1670 K safety margin, not a fitted limit. This matters
because janaf **aborts** the run if any cell leaves [Tlow, Thigh].

### Transport fit provenance

`polynomial`, not `sutherland`, because Sutherland underpredicts air
viscosity by 25 percent at 3000 K and 34 percent at 6000 K, which would
directly corrupt both the wall heat flux and the boundary-layer thickness
that sets the second Mack mode frequency.

Reference data the degree-7 fits were built against:

- T below 500 K: Sutherland, `1.458e-6 T^1.5/(T+110.4)`
- T above 1000 K: Blottner curve fits for N2/O2/Ar with the Wilke mixing rule
- 500-1000 K: smoothstep blend (the two models agree to 0.7 percent at 500 K
  and 0.6 percent at 1000 K)
- kappa: modified Eucken, `kappa = mu*Cv*(1.32 + 1.77*R/Cv)`, with Cv from
  the janaf polynomials

Fit accuracy over 200-6000 K: mu max relative error 1.06 percent, kappa max
1.41 percent, both worst at the 200 K endpoint. Both polynomials stay
strictly positive and monotonic over 150-6500 K.

## Boundary conditions

| Patch | p | U | T |
|---|---|---|---|
| `Inlet` | `freestreamPressure` 292.12 | `freestream` (3169.7 0 0) | `inletOutlet` **250 K** |
| `Atmosphere` | `freestreamPressure` 292.12 | `freestream` (3169.7 0 0) | `inletOutlet` **250 K** |
| `Outlet` | `waveTransmissive` 292.12 | `inletOutlet` | `inletOutlet` **250 K** |
| `Solid_Walls` | `zeroGradient` | `noSlip` | `fixedValue` **1500 K** isothermal |

The 1500 K cold wall against a 4331 K stagnation temperature is what makes
the second Mack mode the dominant instability here: wall cooling damps the
first mode and amplifies the second.

## Numerics

| Item | Value |
|---|---|
| `fluxScheme` | **`Kurganov`** (Kurganov-Noelle-Petrova central-upwind, KNP) |
| `reconstruct(rho)` | **`vanLeer`** |
| `reconstruct(U)` | **`vanLeerV`** |
| `reconstruct(T)` | **`vanLeer`** |
| `ddtSchemes` | `Euler` |
| `laplacianSchemes` | `Gauss linear corrected` |
| Turbulence | **laminar / Stokes, ILES** |

`Kurganov` rather than `Tadmor`: KNP uses one-sided local wave speeds, so its
dissipation collapses towards zero in smooth regions. That is what lets the
second-mode wavepacket survive advection down 48 m of boundary layer instead
of being numerically damped out. The vanLeer limiters are the only subgrid
dissipation in the run - that is the ILES closure, and it is why no
Smagorinsky or WALE model must be added.

## Time stepping

| Item | Value |
|---|---|
| `adjustTimeStep` | **no** - fixed step |
| `deltaT` | **4.5e-08 s (45 ns)** |
| Total steps | **675,000** |
| `endTime` | **0.030375 s** |
| Target CFL | **0.37** |
| Actual CFL at 0.5 mm | 0.314 |

Courant working:

```
fastest wave speed   |u| + a = 3169.7 + 317.4 = 3487 m/s
mesh minimum size    0.5 mm   (Mesh.MeshSizeMin, scripts/gmsh_config.json)
Co = 3487 * 45e-9 / 0.5e-3 = 0.314          about 15 percent under target

inverting for the target:
dx_min = 3487 * 45e-9 / 0.37 = 4.241e-4 m = 0.424 mm
```

**The 0.5 mm mesh floor is what makes the fixed step safe.** Refine the mesh
below 0.424 mm anywhere and the timestep must come down with it - there is no
adaptive step to catch it. `maxCo 0.37` in `controlDict` is documentation
only while `adjustTimeStep no`; the running value is logged every 100 steps
by the `CourantNo` / `courantMonitor` pair.

The fixed step is required by the PSE analysis, not merely convenient: an
adaptive step produces unevenly spaced probe samples, which cannot be FFT'd
without resampling and which smears the narrow second-mode peak.

Physical coverage: 0.030375 s over the 48.061 m body at 3169.7 m/s is
**2.0 flow-through times**. Enough to establish the shock layer and give a
30.4 ms probe record, but short for converged turbulence statistics.

## Output

| Output | Cadence | Count |
|---|---|---|
| Volume fields | every 33,750 steps (1.51875 ms) | **20** |
| Shock surfaces and planes | every **1000** steps (45 us) | 675 |
| Wall surface (p, heat flux, shear) | every 1000 steps | 675 |
| PSE wall probes | every **10** steps (450 ns) | 67,500 rows |
| BL rakes | every 100 steps (4.5 us) | 6,750 rows |
| Forces, coefficients, extrema | every 100 steps | 6,750 rows |

`writeFormat binary`, `writeCompression off`, `purgeWrite 0`.

## Parallel and file handling

| Item | Value |
|---|---|
| `numberOfSubdomains` | **384** |
| `method` | `ptscotch` |
| Cells per rank | about 81,000 (31M / 384) |
| `fileHandler` | **`collated`** |
| `maxThreadFileBufferSize` | 0 (synchronous writes) |

`ptscotch` over a geometric split because cell density varies by more than
three orders of magnitude across this mesh (0.5 mm prisms against
metre-scale far-field tets), so a bounding-box decomposition cannot balance
it. The `preservePatches` constraint the previous configuration applied to
`Atmosphere` has been removed: it exists for coupled patches (cyclic,
cyclicAMI) and had nothing to act on for an ordinary external boundary.

Collated turns roughly 7,700 files per volume write into 20, removes the
metadata storm on a parallel filesystem at 384 ranks, and means
`reconstructPar` is never run - ParaView reads the collated case directly.

## PSE sampling

```
sample interval  dt_s = 10 * 45 ns = 450 ns
sampling rate    fs   = 2.2222 MHz
NYQUIST          fN   = 1.1111 MHz
record length    T    = 0.030375 s
resolution       df   = 32.92 Hz
samples          N    = 67,500 per probe
```

Second Mack mode frequency `f2 = u_e / (2*delta)`, with `delta` from the
Eckert reference-temperature method for a cold-wall laminar layer
(`Tw` 1500 K, `Te` 250 K, `Me` 9.99, giving `T*` = 1794 K and
`Re*/x` = 2.90e4 per metre):

| x [m] | delta [mm] | f2 [kHz] |
|---|---|---|
| 0.10 | 9.3 | 170 |
| 0.25 | 14.7 | 108 |
| 1.0 | 29.4 | 54 |
| 3.0 | 50.9 | 31 |
| 10.0 | 92.9 | 17 |
| 24.0 | 143.9 | 11 |
| 47.0 | 201.3 | 7.9 |

Band to capture: about 8 kHz to 250 kHz. **Nyquist at 1.111 MHz is at least
4.4x the top of that band**, so no anti-alias filtering or resampling is
needed before the FFT, and `df` = 32.9 Hz resolves the narrow second-mode
peak with hundreds of bins.

Spatial check: the shortest resolved wave is `lambda = 0.9*u_inf/f`. At
250 kHz that is 11.4 mm, or 23 cells at the 0.5 mm mesh minimum. Above about
500 kHz the mesh, not the sampling rate, becomes the limit.

Probe counts: 26 windward wall stations, 26 leeward, 14 spanwise (two rakes
of 7), and 8 wall-normal boundary-layer rakes of 241 points each.

## Geometry reference

| Quantity | Value |
|---|---|
| Vehicle X extent | 0.057 to 49.643 m (`lRef` 48.061 m as-meshed) |
| Vehicle Y extent | -22.459 to 14.213 m |
| Vehicle Z extent | -0.022 to 2.600 m |
| `CofR` | (25.717, 0.000, 3.022) |
| `Aref` | 912.1 m2 - **confirm against CAD** |
| Far-field box | 200 m, from x = -50 m |
| Mesh | about 31M cells, hybrid prism/tet, y+ below 1 |
