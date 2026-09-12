# Settings - why each choice was made

Numbers live in `docs/PARAMETERS.md`. This file is the reasoning.

## 1. Why rhoCentralFoam and not hy2Foam

The previous configuration ran `hy2Foam` from the hyStrath suite: 11-species
reacting air, two-temperature (translational plus vibrational) thermal
non-equilibrium, Park 1993 chemistry, Blottner-Eucken transport, and three
custom shared libraries that had to be built against v2412 by hand.

That model is the right one when the shock layer is hot enough to dissociate
air and slow enough to relax. Neither holds here. Solving the normal shock
with a thermally perfect janaf gas gives a post-shock temperature of 4256 K
and a stagnation temperature of 4331 K. Oxygen dissociation becomes
significant above roughly 2500 K and is substantial by 4000 K, so there is
some O2 dissociation in the stagnation region - but it is confined to a small
pocket, and nitrogen (79 percent of the gas) is essentially undissociated
below 4000 K. Across the boundary layer that the PSE analysis actually cares
about, the wall is at 1500 K and the edge is far cooler than the stagnation
point; that flow is chemically frozen.

What is *not* negligible is the caloric imperfection: cp rises from
999 J/kg/K at 250 K to 1319 J/kg/K at 4000 K, and gamma falls from 1.404 to
1.278. A constant-gamma perfect gas would get the shock stand-off and the
shock-layer temperature visibly wrong. `janaf` captures exactly that, at a
fraction of the cost and none of the fragility of an 11-species reacting
solve, and it does it with stock OpenFOAM that needs no custom build.

The trade is explicit: this case will underpredict the stagnation-region
temperature somewhat, because in reality dissociation absorbs energy there
that a frozen gas has nowhere to put. If the stagnation-point heating number
is the deliverable, that limitation matters and a reacting model is needed.
If the deliverable is the boundary-layer stability and the shock structure -
which is what this case is instrumented for - the frozen calorically
imperfect gas is the appropriate model.

`rhoCentralFoam` also has to be the solver because a density-based
central-upwind scheme is what captures a Mach 10 bow shock without a
pressure-based solver's acoustic stiffness, and because it is the only stock
OpenFOAM solver whose flux scheme is directly selectable as `Kurganov`.

## 2. Thermo: janaf, Thigh 6000

`Thigh 6000` is not arbitrary and it is not a fitted limit. It is 1670 K
above the highest temperature the case physically reaches (4331 K stagnation).
That headroom is deliberate: **janaf aborts the run** if any cell leaves
`[Tlow, Thigh]`, and a transient overshoot during shock formation in the
first few thousand steps is entirely plausible on an explicit scheme. 6000 K
is also where the standard NASA-7 high-temperature range ends, so extending
the ceiling further would mean extrapolating the polynomial rather than using
it.

`Tlow 200` covers the 250 K freestream with margin for an expansion
overshoot around the wing leading edges.

The coefficients are the mole-fraction-weighted sum of the Burcat N2, O2 and
Ar polynomials. That weighting is exact rather than approximate, because a
NASA-7 fit is molar and the mixture is ideal at fixed composition. cp is
continuous across Tcommon to 0.19 percent.

## 3. Transport: polynomial, valid to 6000 K

Sutherland is the default choice and it is wrong for this case. It
underpredicts air viscosity by 25 percent at 3000 K and 34 percent at 6000 K.
Two things depend directly on getting mu(T) right in the shock layer:

- the wall heat flux, which is the TPS sizing output of the run
- the boundary-layer thickness, which sets the second Mack mode frequency
  through `f2 = u_e/(2*delta)` - a 25 percent viscosity error moves the
  frequency the PSE solver is looking for by more than 10 percent

The degree-7 polynomials were fitted against Sutherland below 500 K (where
it is accurate), Blottner curve fits with the Wilke mixing rule above 1000 K
(the standard high-temperature air model), and a smoothstep blend between.
The two reference models agree to better than 1 percent across the blend
region, so the blend introduces no artefact.

kappa comes from the modified Eucken relation with Cv from the janaf
polynomials, which keeps the Prandtl number physically consistent with the
thermo rather than pinning it to a constant. Pr comes out at 0.691 at 250 K
rising to 0.708 at 6000 K.

Both fits hold to better than 1.5 percent over the full 200-6000 K range and
stay strictly positive over 150-6500 K, so a transient excursion just outside
the window cannot produce a negative transport coefficient before janaf traps
it.

## 4. Energy: sensibleInternalEnergy

`rhoCentralFoam` transports `rhoE`, a total *internal* energy density. With
`sensibleInternalEnergy` the thermo works in the same variable, so the
conversion between the conserved quantity and the thermodynamic state is
direct. `sensibleEnthalpy` would work but adds a `p/rho` round-trip at every
cell every step, for no benefit in a non-reacting case. `absolute` variants
are only needed when species with different heats of formation are being
created or destroyed, which is not happening here.

The practical consequence: the solved energy field is named `e`, which is why
`system/fvSolution` has an `"(e|h)"` block.

## 5. Flux: Kurganov with vanLeer

`fluxScheme Kurganov` is the Kurganov-Noelle-Petrova central-upwind scheme.
The alternative stock option is `Tadmor` (Kurganov-Tadmor), which uses a
symmetric wave speed and is uniformly more dissipative.

KNP is required here for a specific reason beyond general shock quality. It
estimates one-sided local wave speeds, so in smooth flow the upwind bias -
and with it the numerical dissipation - collapses towards zero. The second
Mack mode is a small-amplitude wave that has to survive advection down 48 m
of boundary layer to reach the downstream probes. Under KT's extra
dissipation it would be damped out before it got there, and the run would
report a stable boundary layer that is not physically stable.

`vanLeer` on every reconstruction: smooth, symmetric, second-order TVD.
Monotone through the shock, and far less dissipative in smooth flow than
`minmod`. `vanLeerV` on velocity is the vector form - it applies one limiter
value to all three components based on the most-limited direction, so the
limiter cannot rotate the reconstructed velocity vector.

## 6. Turbulence: laminar / Stokes, ILES

`simulationType laminar` with `model Stokes` selects the plain Newtonian
stress model. The deviatoric stress is built purely from the molecular mu(T)
above. No eddy viscosity, no subgrid model, no wall function, and no `nut` or
`alphat` field exists.

The only other dissipation in the run is the vanLeer limiter dissipation from
section 5. **That limiter dissipation is the implicit subgrid model.** Adding
a Smagorinsky or WALE block would double-count it.

There is a second reason this has to be laminar, specific to this case: PSE
linearises about a *laminar* mean flow. A modelled eddy viscosity in the
baseflow would change the mean profile the stability analysis is performed
on, and the resulting N-factors would not mean anything.

## 7. Temperatures

Freestream, inlet and outlet at 250 K is the ISA value at 130,000 ft.
`inletOutlet` rather than `fixedValue` on all three: the inflow is supersonic
so the full state is imposed, but `inletOutlet` keeps the boundary well posed
if a start-up transient briefly reverses a face, and reverts to `fixedValue`
the moment inflow is restored.

The wall at 1500 K isothermal is a radiative-equilibrium TPS surface
estimate. Plain `fixedValue` is correct - the flow is dense continuum at
130,000 ft on a y+ below 1 mesh, so no temperature-jump or slip condition is
needed. (The previous configuration noted this too, in rejecting the
rarefied `nonEqSmoluchowskiJumpT` boundary condition from the hyStrath
examples.)

The 1500 K wall is also not incidental to the PSE work. Against a 4331 K
stagnation temperature this is a strongly cooled wall, and wall cooling damps
the first Mack mode while amplifying the second. That is precisely why the
second mode is the mode worth instrumenting on this vehicle.

## 8. Time stepping: fixed 45 ns

A fixed step is a deliberate requirement, not a fallback. The probe time
series that feed the PSE analysis must be uniformly sampled: an adaptive step
gives unevenly spaced samples that cannot be FFT'd without resampling, and
resampling smears the narrow second-mode peak that the whole exercise is
trying to measure.

The step is safe because of the mesh, and only because of the mesh:

```
fastest wave speed   |u| + a = 3169.7 + 317.4 = 3487 m/s
mesh minimum size    0.5 mm
Co = 3487 * 45e-9 / 0.5e-3 = 0.314
```

against the 0.37 target, with about 15 percent of headroom. Inverting the
target gives a hard floor: **no cell may be smaller than 0.424 mm**. If the
mesh is ever refined past that, `deltaT` must come down with it. Nothing in
the case will catch this automatically - `maxCo` is documentation while
`adjustTimeStep` is `no`. The `CourantNo` / `courantMonitor` pair logs the
running maximum every 100 steps; that log is the safety net.

675,000 steps is 0.030375 s, or 2.0 flow-through times over the 48.061 m
body. That is enough to establish the shock layer and to give a 30.4 ms probe
record with 32.9 Hz frequency resolution. It is short for converged
turbulence statistics, which is why `fieldAverage` here is documented as a
laminar mean baseflow for PSE and not as an LES statistical average.

## 9. Output cadences

Three different rates, each set by what it feeds:

**Volume fields, 20 writes.** These are 31M-cell snapshots and they exist for
three-dimensional inspection and for restart, not for time-resolved analysis.
20 is what the request specified and it is the right order: at 1.52 ms apart
they sample the run coarsely but cover it completely.

**Surfaces, every 1000 steps.** 45 us per frame, 675 frames. Surfaces are two
orders of magnitude cheaper than a volume write, which is the whole reason to
sample the shock this way rather than reconstructing it from 20 snapshots. At
45 us the shock system convects 143 mm between frames, so unsteady shock
motion and shock-shock interaction at the leading edges are time-resolved.

**Probes, every 10 steps.** Set by Nyquist, covered in section 10.

## 10. PSE probes and the second Mack mode

The second Mack mode is a trapped acoustic wave in the boundary layer, with
`f2 ~ u_e/(2*delta)`. Because delta grows down the body, the frequency sweeps
downward: using the Eckert reference-temperature method for this cold-wall
layer (`T*` = 1794 K, `Re*/x` = 2.90e4 per metre), f2 runs from about 170 kHz
at x = 0.1 m to 7.9 kHz at x = 47 m. The band to capture is roughly
8 kHz to 250 kHz.

Sampling every 10 steps at 45 ns gives fs = 2.2222 MHz and a **Nyquist
frequency of 1.1111 MHz** - at least 4.4x the top of that band. No
anti-aliasing filter and no resampling is needed. The 30.375 ms record gives
df = 32.9 Hz, which resolves the narrow second-mode peak with hundreds of
bins, and 67,500 samples per probe, which is a comfortable number for Welch
averaging.

The spatial limit is the mesh, not the sampling rate: at 250 kHz the
convected wavelength is 11.4 mm, or 23 cells at 0.5 mm. Above about 500 kHz
the mesh runs out first.

Wall pressure is the sampled quantity because it is the standard experimental
observable for the second mode (what a PCB or PVDF sensor measures in a quiet
tunnel) and the quantity PSE N-factors are conventionally reported against.
`patchProbes` rather than `probes`: it snaps to the nearest face centre on
`Solid_Walls`, so the sample is exactly on the wall regardless of how well
the nominal coordinate matches the real surface. A plain `probes` object
would silently land inside the solid or out in the shock layer.

The boundary-layer rakes are sampled at 100 steps instead of 10. They exist
for wall-normal structure - the mean profiles PSE linearises about, and the
eigenfunction shape that confirms the disturbance really is the second mode -
not for spectra. 8 rakes of 241 points at the 10-step rate would be 1.6e8
samples for no analytical gain.

Spanwise rakes at x = 20 m and x = 35 m measure the spanwise wavenumber beta.
`pse_solver/config/mach10_130kft.json` currently runs `beta = 0`, a purely
two-dimensional second mode. If those rakes show coherent phase variation
across the span, that assumption is wrong for this geometry and `n_beta`
needs to be opened up.

## 11. Collated output

At 384 ranks, uncollated writing produces roughly 7,700 files per volume
write. Collated produces 20. On a parallel filesystem the metadata operations
from that many small files dominate the write time - the bytes are not the
problem, the file count is.

Collated also removes `reconstructPar` from the workflow entirely. ParaView
reads the collated case directly through `case.foam`, and `postProcessing/`
output is already gathered on rank 0. Reconstructing a 31M cell case 20 times
would cost hours and produce nothing that was not already readable.

`maxThreadFileBufferSize 0` keeps the writes synchronous. Non-zero enables a
background write thread, which is faster but requires `MPI_THREAD_MULTIPLE`
from the MPI build - do not enable it without knowing the site MPI supports
that.

Every parallel utility must agree on the handler, so it is set in
`system/controlDict` rather than per-command, and `-fileHandler collated` is
passed explicitly by the scripts as well.

## 12. Shock sampling for ParaView

Four different shock representations are written, because no single one is
sufficient:

**`|grad rho|` isosurfaces** are the primary capture. Unlike a pressure
isovalue, the density-gradient magnitude does not need retuning when the
freestream changes, and it picks up the weak oblique shocks off the wing
leading edges as well as the strong bow shock. Two thresholds are written
(5.0 and 1.0 kg/m4) because the shock weakens and the cells coarsen going
aft, so one threshold cannot follow it the whole way. Both values are
estimates from the 7.13x normal-shock density jump and should be retuned
against the first written frame.

**Pressure isosurfaces** at 2,000 Pa and 17,800 Pa. The first is the leading
edge of the compression, which is what shock stand-off is measured from. The
second is the mid-jump `(p1+p2)/2`, which is the conventional definition of
"the shock surface" when comparing against inviscid shock fitting or a Billig
correlation.

**The sonic surface** `Ma = 1` bounds the subsonic pocket behind the
near-normal part of the bow shock. It is cheap, small, and it is the surface
that reveals whether the nose region has gone unsteady.

**Cutting planes**: the y = 0 symmetry plane (the single most useful view -
bow shock, shock layer, boundary layer and wake in one image), a waterline
plane through the vehicle mid-thickness for the planform shock pattern, and
six cross-flow stations at 8 m spacing for shock shape and spanwise
shock-shock interaction. Every plane carries `bounds` to clip it to the near
field; without that, each would drag along a large sheet of undisturbed
freestream from the 200 m far-field box.

Format is `vtp` (VTK XML PolyData, binary) rather than legacy ASCII `vtk`:
roughly 4x smaller, and ParaView does not have to re-parse it every frame.

Budget disk for this. 675 frames of five isosurfaces through a 31M cell mesh
can reach tens of gigabytes. If that is too much, the cheapest reductions in
order are: drop `shockIsoP_outer`, tighten the plane `bounds`, then halve the
frame rate.
