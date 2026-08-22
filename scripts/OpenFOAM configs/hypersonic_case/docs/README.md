# hy2Foam Case: Mach 10 Spaceplane, 130,000 ft, Implicit LES, 11-Species Real Gas

## What this is

A hy2Foam (hyStrath suite, via https://github.com/ivanZanardi/hypersonicfoam)
OpenFOAM case for a hypersonic vehicle at Mach 10, 130,000 ft altitude,
-5 deg angle of attack, using:

- Implicit LES (ILES) only - no explicit subgrid closure
- KNP (Kurganov-Noelle-Petrova) flux scheme, Minmod limiter on U,
  vanLeer limiter on everything else
- Isothermal wall at 1500 K
- 11-species real gas thermochemistry (Park 1993 mechanism), with
  translational/vibrational temperature non-equilibrium (Tt, Tv)
- Mesh: hybrid (prism boundary layer + tet core), y+ < 1

The actual hypersonicfoam repo was cloned and inspected to build this case,
so the file structures, dictionary keys, and chemistry/thermo data below are
checked against the real solver source rather than assumed. See
`SETTINGS.md` for a full list of corrections made after that check.

The mesh itself was not uploaded (upload size limits), so this package
ships without `constant/polyMesh/`. You must convert your own mesh before
running - see "The one placeholder left" below.

## Freestream conditions used

| Quantity | Value |
|---|---|
| Altitude | 130,000 ft (39.6 km) |
| Pressure P_inf | 292.12 Pa |
| Density rho_inf | 0.004071 kg/m^3 |
| Temperature T_inf (= Tt = Tv) | 249 K |
| Speed of sound a_inf | 316.97 m/s |
| Mach number | 10 |
| Velocity u_inf | 3169.7 m/s |
| Wall temperature | 1500 K (isothermal) |
| AoA | -5 deg (built into mesh geometry, not the velocity vector) |

Full derivation and rationale in `PARAMETERS.md` and `SETTINGS.md`.

## The one placeholder left: the mesh

Everything else in this package - chemistry mechanism, species thermo data,
V-T relaxation model, transport/collision data - is now the REAL data
copied verbatim from the hypersonicfoam repo, not a placeholder. The only
thing you must edit is:

**`scripts/setup.sh`** - `MESH_FILE` variable near the top. Set this to
the actual path of your `.msh` file, then run `bash scripts/setup.sh`.

After conversion, also check `constant/polyMesh/boundary` and confirm the
patch names match what `case/0/*` boundary conditions assume (`farfield`,
`wall`, `outlet`) - the repo's own example case uses different names
(`inlet`, `cylinder`) for its own geometry, so don't assume ours match
your mesh either. Rename patches in Gmsh, or edit the `0/` files, as
needed.

## What ships pre-filled (not placeholders)

- `constant/chemDicts/hTCReactionsEarth93` - real Park (1993) 11-species
  air reaction mechanism, copied from the repo's `genericCase` example
- `constant/thermoDEM` - real per-species thermodynamic data
- `constant/thermo2TModel` - real Millikan-White/Park V-T relaxation model
- `constant/hTCProperties`, `constant/chemistryProperties` - real config,
  copied from the repo
- `constant/transportProperties` - real Gupta/Yos/Thompson 11-species
  collision-integral data (rarefaction diagnostics turned off since this
  is dense continuum flow, not the near-continuum-breakdown regime those
  diagnostics target)
- `constant/thermophysicalProperties` references all of the above via
  `$FOAM_CASE`-relative paths, so it works immediately once you've placed
  this case anywhere and run `setup.sh`

## Running

```bash
tar xzf hypersonic_case_*.tar.gz
cd hypersonic_case
# Edit scripts/setup.sh: MESH_FILE only
bash scripts/setup.sh
bash scripts/run.sh 8       # 8 MPI ranks, adjust as needed
bash scripts/postProcess.sh
```

## Package contents

```
case/0/            Initial conditions (U, p, Tt, Tv, 11 species Y fields)
case/constant/      transportProperties, thermophysicalProperties,
                    chemistryProperties, hTCProperties, thermoDEM,
                    thermo2TModel, chemDicts/hTCReactionsEarth93,
                    turbulenceProperties, polyMesh/ (empty - see above)
case/system/        controlDict, fvSchemes, fvSolution
scripts/            setup.sh, run.sh, postProcess.sh, cleanup.sh
docs/               this file, PARAMETERS.md, SETTINGS.md, VALIDATION.md
logs/               run logs written here by run.sh
```

## Important caveats

- Boundary patch names (`farfield`, `wall`, `outlet`) are placeholders -
  confirm against your actual converted mesh.
- `endTime` in `controlDict` is a placeholder judgment call (0.02 s) - set
  it based on your vehicle's reference length and desired flow-through
  time once known.
- Mesh classification (hybrid, y+<1) and complexity level (complex,
  leaning down from "complex/extreme") were set based on your description
  of a "smooth spaceplane" body - see `SETTINGS.md` for the reasoning.
- This case targets OpenFOAM v1706, which is what hyStrath (and therefore
  hy2Foam) is built against per the repo's own README - `fvConstraints`
  (introduced later, ~v1812) is intentionally not used here since this
  solver predates it.
