# Parameters and Derivations

## Freestream state (as supplied, not from ISA table lookup)

130,000 ft = 39,624 m (39.6 km). This falls between the 30 km and 40 km rows of the standard 1976 US Standard Atmosphere table:

```
P_inf   = 292.12 Pa
rho_inf = 0.004071 kg/m^3
T_inf   = 249 K
a_inf   = 316.97 m/s
```

### Consistency check (ideal gas law)

```
R = P_inf / (rho_inf * T_inf)
  = 292.12 / (0.004071 * 249)
  = 292.12 / 1.0137
  = 288.16 J/(kg.K)
```

Standard air R = 287.05 J/(kg.K). The 0.4% difference is well within normal real-atmosphere variation and is not a concern - the inputs are self-consistent.

### Speed of sound cross-check

```
a = sqrt(gamma * R * T) = sqrt(1.4 * 287.05 * 249) = 316.33 m/s
```

vs. assumed 316.97 m/s - a 0.2% difference, consistent with the same minor R variation above. The supplied value (316.97 m/s) was used directly in all downstream calculations rather than the recomputed one.

### Velocity

```
u_inf = M * a_inf = 10 * 316.97 = 3169.7 m/s
```
That's mach 10

### Dynamic viscosity (Sutherland's law)

```
mu(T) = mu_ref * (T/T_ref)^1.5 * (T_ref + S) / (T + S)

mu_ref = 1.716e-5 Pa.s   (at T_ref = 288.15 K)
S      = 110.4 K          (Sutherland constant, air)
T      = 249 K

(T/T_ref)       = 249 / 288.15 = 0.8641
(T/T_ref)^1.5   = 0.8032
(T_ref+S)/(T+S) = 398.55 / 359.4 = 1.1089

mu_inf = 1.716e-5 * 0.8032 * 1.1089 = 1.529e-5 Pa.s
```

### Kinematic viscosity

```
nu_inf = mu_inf / rho_inf = 1.529e-5 / 0.004071 = 3.756e-3 m^2/s
```

### Reynolds number

Not computed now, however can be computed using following formula:

```
Re_L = (rho_inf * u_inf * L) / mu_inf
```

## Post-shock estimate (documentation only, not used in initial conditions)

Normal-shock relation (calorically perfect gas, for order-of-magnitude reference only - the actual 11-species real-gas will differ):

```
T2/T1 = [(2*gamma*M1^2 - (gamma-1)) * ((gamma-1)*M1^2 + 2)] / [(gamma+1)^2 * M1^2]

M1 = 10, gamma = 1.4

Numerator   = (2*1.4*100 - 0.4) * (0.4*100 + 2) = 279.6 * 42 = 11,743.2
Denominator = (2.4)^2 * 100 = 576

T2/T1 = 20.39
T2 = 249 * 20.39 = 5077 K
```

This ~5000 K estimate is used in `SETTINGS.md` to reason about whether ionization is likely to matter at this flight condition (probably not significantly - ionization onset is typically above 8000-10000 K).

## Species composition (freestream)

Standard undissociated air, mass basis:

```
Y_N2 = 0.767
Y_O2 = 0.233
all other species (NO, N, O, N2+, O2+, NO+, N+, O+, e-) = 0.0
```

Sum = 1.000, satisfying mass fraction conservation. Field file names in
`case/0/` match these species names exactly (including the literal `+`
and `-` characters), because that's what `constant/chemDicts/
hTCReactionsEarth93` and `constant/thermoDEM` (both copied verbatim from
the real hy2Foam repo) use.

## Wall temperature

1500 K, isothermal, applied via `fixedValue` to both `Tt` and `Tv` at the `wall` patch (full vibrational accommodation assumed at the wall).

## AoA handling

-5 deg AoA is represented by the mesh geometry (the vehicle is rotated relative to the domain's freestream axis during meshing). The freestream velocity vector in `0/U` is therefore purely axial: `(3169.7 0 0)`.
