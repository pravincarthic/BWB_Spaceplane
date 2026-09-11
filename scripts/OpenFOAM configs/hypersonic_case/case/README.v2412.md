This case targets OpenCFD OpenFOAM v2412 with the fully ported
ivanZanardi/hypersonicfoam hyStrath suite. It requires `hy2Foam`,
`libstrathFiniteVolume.so`, and `libstrathForces.so` built for that installation.
The v2412 port was not supplied for inspection; the custom dictionary keys and
library names follow the original fork. Any renaming in the local port must be
reflected here.

The migration repairs 30 unterminated header comments, selects the laminar
model with `model Stokes;`, enables compression with `writeCompression on;`,
and loads the fork's force objects, which support `multi2Thermo`. Explicit
`isoMethod point;` retains the legacy shock surface algorithm because v2412
defaults to `topo`.

`FoamFile.version` remains `2.0`: it identifies the file format. The case keeps
the custom `chemistrySolver`, `chemistryThermo`, `hTC2Model`, thermo types,
`thermo:psi`, all 11 species and their `Tv_*` fields. Initial values, boundary
conditions, reaction and transport data, numerical schemes, time controls and
decomposition settings are preserved. The existing `surfaces` dictionary form
and `libs` entries are supported by v2412.

After sourcing the v2412 environment and the ported solver environment, run
these read checks from this case directory in Bash:

```bash
foamVersion
command -v hy2Foam
export FOAM_CASE="$PWD"
(
    set -e
    for dictionary in 0/* constant/chemDicts/* constant/chemistryProperties \
        constant/hTCProperties constant/thermo2TModel constant/thermoDEM \
        constant/thermophysicalProperties constant/transportProperties \
        constant/turbulenceProperties system/*
    do
        foamDictionary "$dictionary" -keywords >/dev/null
    done
)
```

Native parsing and solver execution were not performed during migration:
OpenFOAM and the ported solver were unavailable, and `constant/polyMesh`
contains only a placeholder. Once the mesh is available, run
`checkMesh -allTopology -allGeometry` and verify that its patches include
`Inlet`, `Atmosphere`, `Solid_Walls` and `Outlet`. Check a short run in a copy
of the case, including force output and shock surface output. A full
decomposition uses the existing 384 subdomains; the sibling run wrapper's
default of 8 ranks does not match this dictionary.

Source checks used the original fork at
[`98e8db2`](https://github.com/ivanZanardi/hypersonicfoam/tree/98e8db22ceca2b816714abb13fcc7e51ed6c09e0)
and OpenCFD v2412 source mirrored at
[`88fd2b2`](https://github.com/JiangHuShrimp/OpenFoam-v2412/tree/88fd2b21dad72fef3fecb6bd79b7c1947c186e6c):

* `src/TurbulenceModels/turbulenceModels/laminar/laminarModel/laminarModel.C`
  reads `model` and accepts `laminarModel` as a compatibility alias.
* `src/OpenFOAM/db/options/IOstreamOption.C` reads compression as a switch.
* `src/sampling/sampledSurface/isoSurface/sampledIsoSurface.H` documents the
  default algorithm and available `isoMethod` values.
* The fork's `hyStrath/src/functionObjects/forces/Make/files` names
  `libstrathForces`; its `forces/forces.C` handles `multi2Thermo`.

Older notes in the sibling `docs` directory describe the v1706 case.
