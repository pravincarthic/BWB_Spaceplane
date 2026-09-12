#!/usr/bin/env python3
"""Regenerate case/system/probesPSE from the real Solid_Walls surface.

The probe coordinates shipped in case/system/probesPSE are nominal: they are
derived from the vehicle bounding box and rely on patchProbes snapping each
point to the nearest wall face. That is good enough to start a run, but the
boundary-layer rakes in particular need to leave the wall along the true
surface normal, which a bounding box cannot give.

This script extracts the Solid_Walls patch as an STL, reads the real facet
centroids and normals, and writes an exact probesPSE:

  - patchProbes points sitting on real windward and leeward wall faces at each
    requested streamwise station
  - wall-normal rakes that start a fraction of a millimetre off the wall and
    run outward along the true surface normal, with points geometrically
    clustered near the wall so the viscous sublayer is resolved

Run it from the case package root, after scripts/setup.sh has converted the
mesh:

    python3 scripts/make_pse_probes.py

Options worth knowing:

    --dry-run       print what would be written, change nothing
    --stations ...  override the streamwise station list, in metres
    --rake-length   rake length in metres (default 0.60)
    --stl PATH      skip the extraction and use an STL you already have

Only the Python standard library is used, so this runs anywhere the OpenFOAM
environment does.
"""

import argparse
import math
import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.dirname(HERE)
CASE = os.path.join(PKG, "case")

# Streamwise stations, metres. Clustered near the nose: that is where the
# second Mack mode first destabilises and where its frequency is highest, so
# it is where the station spacing has to be finest.
DEFAULT_STATIONS = [
    0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0,
    10.0, 12.0, 14.0, 17.0, 20.0, 23.0, 26.0, 29.0, 32.0, 35.0,
    38.0, 41.0, 44.0, 47.0,
]

# Stations carrying a spanwise rake, and the spanwise offsets on each.
SPAN_STATIONS = [20.0, 35.0]
SPAN_OFFSETS = [-12.0, -8.0, -4.0, 0.0, 4.0, 8.0, 12.0]

# Stations carrying a wall-normal boundary-layer rake.
RAKE_STATIONS = [1.0, 3.0, 6.0, 10.0, 16.0, 24.0, 34.0, 44.0]

RAKE_POINTS = 241
RAKE_FIRST = 2.0e-5      # first point 20 microns off the wall
RAKE_LENGTH = 0.60       # covers the 201 mm boundary layer at x=47 m plus the
                         # entropy layer above it

WALL_PATCH = "Solid_Walls"


# ---------------------------------------------------------------------------
# STL extraction and parsing
# ---------------------------------------------------------------------------

def extract_stl(case_dir, out_path):
    """Write the wall patch to an ASCII STL using whichever utility exists."""
    candidates = [
        ["surfaceMeshExtract", "-patches", "(%s)" % WALL_PATCH,
         "-latestTime", out_path],
        ["surfaceMeshTriangulate", "-patches", "(%s)" % WALL_PATCH,
         "-latestTime", out_path],
    ]
    for cmd in candidates:
        if shutil.which(cmd[0]) is None:
            continue
        print("running: %s" % " ".join(cmd))
        res = subprocess.run(cmd, cwd=case_dir)
        if res.returncode == 0 and os.path.exists(out_path):
            return out_path
        print("  %s failed with code %d, trying the next one"
              % (cmd[0], res.returncode))
    raise SystemExit(
        "ERROR: could not extract the wall surface.\n"
        "Neither surfaceMeshExtract nor surfaceMeshTriangulate ran "
        "successfully.\n"
        "Source the OpenFOAM v2412 environment first, or pass an STL you "
        "already have with --stl."
    )


def read_stl(path):
    """Return [(centroid, unit normal, area)] from an ASCII or binary STL."""
    with open(path, "rb") as fh:
        head = fh.read(5)
        fh.seek(0)
        raw = fh.read()
    if head == b"solid" and b"facet normal" in raw[:4096]:
        return _read_stl_ascii(raw.decode("utf-8", "replace"))
    return _read_stl_binary(raw)


def _facet(v0, v1, v2):
    cx = ((v0[0] + v1[0] + v2[0]) / 3.0,
          (v0[1] + v1[1] + v2[1]) / 3.0,
          (v0[2] + v1[2] + v2[2]) / 3.0)
    ax = (v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2])
    bx = (v2[0] - v0[0], v2[1] - v0[1], v2[2] - v0[2])
    n = (ax[1] * bx[2] - ax[2] * bx[1],
         ax[2] * bx[0] - ax[0] * bx[2],
         ax[0] * bx[1] - ax[1] * bx[0])
    mag = math.sqrt(n[0] ** 2 + n[1] ** 2 + n[2] ** 2)
    if mag <= 0.0:
        return None
    area = 0.5 * mag
    return cx, (n[0] / mag, n[1] / mag, n[2] / mag), area


def _read_stl_ascii(text):
    facets, verts = [], []
    for line in text.splitlines():
        parts = line.split()
        if not parts:
            continue
        if parts[0] == "vertex" and len(parts) == 4:
            verts.append(tuple(float(v) for v in parts[1:4]))
        elif parts[0] == "endfacet":
            if len(verts) == 3:
                f = _facet(*verts)
                if f:
                    facets.append(f)
            verts = []
    return facets


def _read_stl_binary(raw):
    import struct
    if len(raw) < 84:
        raise SystemExit("ERROR: STL file is truncated.")
    count = struct.unpack("<I", raw[80:84])[0]
    facets = []
    off = 84
    for _ in range(count):
        if off + 50 > len(raw):
            break
        vals = struct.unpack("<12f", raw[off:off + 48])
        f = _facet(vals[3:6], vals[6:9], vals[9:12])
        if f:
            facets.append(f)
        off += 50
    return facets


# ---------------------------------------------------------------------------
# Station selection
# ---------------------------------------------------------------------------

def pick(facets, x, y, side, x_tol):
    """Nearest facet to (x, y) on the given side.

    side is +1 for leeward (outward normal has a positive z component) and -1
    for windward. Facets whose normal is nearly parallel to the flow are
    rejected: those are leading/trailing edge slivers, and a probe on one of
    them would be measuring an edge rather than a boundary layer.
    """
    best, best_d = None, None
    for c, n, area in facets:
        if side * n[2] <= 0.25:      # want a surface facing up or down
            continue
        if abs(c[0] - x) > x_tol:
            continue
        d = (c[0] - x) ** 2 + (c[1] - y) ** 2
        if best_d is None or d < best_d:
            best, best_d = (c, n, area), d
    return best


def cluster(first, total, npoints):
    """Geometric point spacing: first cell `first`, total length `total`."""
    lo, hi = 1.0 + 1e-9, 1.5
    for _ in range(200):                      # bisect on the growth ratio
        mid = 0.5 * (lo + hi)
        s = first * (mid ** (npoints - 1) - 1.0) / (mid - 1.0)
        if s < total:
            lo = mid
        else:
            hi = mid
    r = 0.5 * (lo + hi)
    out, d = [], 0.0
    for i in range(npoints):
        out.append(d)
        d += first * (r ** i)
    return out, r


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

def fmt_pts(points, indent=12):
    pad = " " * indent
    return "\n".join("%s(%.6f %.6f %.6f)" % (pad, p[0], p[1], p[2])
                     for p in points)


def build(facets, stations, rake_length, x_tol):
    xmin = min(c[0] for c, _, _ in facets)
    xmax = max(c[0] for c, _, _ in facets)
    print("wall surface: %d facets, x from %.3f to %.3f m"
          % (len(facets), xmin, xmax))

    wind, lee, missing = [], [], []
    for x in stations:
        if not (xmin <= x <= xmax):
            missing.append((x, "outside the wall x range"))
            continue
        w = pick(facets, x, 0.0, -1, x_tol)
        l = pick(facets, x, 0.0, +1, x_tol)
        if w:
            wind.append(w[0])
        else:
            missing.append((x, "no windward facet"))
        if l:
            lee.append(l[0])
        else:
            missing.append((x, "no leeward facet"))

    span = []
    for x in SPAN_STATIONS:
        for y in SPAN_OFFSETS:
            hit = pick(facets, x, y, -1, x_tol)
            if hit:
                span.append(hit[0])

    offs, ratio = cluster(RAKE_FIRST, rake_length, RAKE_POINTS)
    print("rake spacing: first %.1f um, last %.2f mm, growth ratio %.4f"
          % (RAKE_FIRST * 1e6, (offs[-1] - offs[-2]) * 1e3, ratio))

    rakes = []
    for x in RAKE_STATIONS:
        hit = pick(facets, x, 0.0, -1, x_tol)
        if hit is None:
            missing.append((x, "no windward facet for a rake"))
            continue
        c, n, _ = hit
        # Step outward along the true wall normal. The STL normal points
        # out of the solid and into the fluid, so -n would walk into the
        # body; the sign below is chosen to move away from the wall.
        pts = [(c[0] + n[0] * d, c[1] + n[1] * d, c[2] + n[2] * d)
               for d in offs]
        rakes.append((x, pts))

    return wind, lee, span, rakes, missing


HEADER = """// -*- C++ -*- ---------------------------------------------------------------
// system/probesPSE
//
// GENERATED by scripts/make_pse_probes.py from the real Solid_Walls surface.
// Edit the station lists at the top of that script and re-run it rather than
// editing this file by hand.
//
// Included from the functions{{ }} block of system/controlDict.
// NO FoamFile header: every top-level entry becomes a function object.
//
// Sampling every 10 solver steps at the fixed 45 ns step:
//     dt_s = 450 ns    fs = 2.2222 MHz    NYQUIST = 1.1111 MHz
//     T    = 0.030375 s   df = 32.92 Hz   N = 67,500 samples per probe
// The second Mack mode band on this vehicle runs from about 7.9 kHz at the
// trailing edge to about 170 kHz at x = 0.10 m (f2 = u_e/2delta, Eckert
// reference-temperature laminar boundary layer, Tw 1500 K, Te 250 K,
// Me 9.99). The Nyquist frequency is therefore at least 6x the highest
// frequency of interest, so no anti-alias filtering or resampling is needed
// before the FFT.
//
// Wall stations used: {nwind} windward, {nlee} leeward, {nspan} spanwise.
// Boundary-layer rakes: {nrake}, each {npts} points over {rlen:.3f} m along the
// true wall normal, first point {first:.1f} um off the wall.
// ---------------------------------------------------------------------------
"""

PROBE_BLOCK = """
{name}
{{
    type            patchProbes;
    libs            (sampling);
    patches         (Solid_Walls);

    executeControl  timeStep;
    executeInterval 10;
    writeControl    timeStep;
    writeInterval   10;

    fields          (p T U rho);

    // {comment}
    probeLocations
    (
{points}
    );
}}
"""

RAKE_HEADER = """
// ===========================================================================
// Boundary-layer rakes - PSE baseflow profiles and mode shapes
// ===========================================================================
// Sampled every 100 steps (4.5 us, fs 222 kHz, Nyquist 111 kHz). That is
// still above the second-mode frequency at every station aft of about x = 3 m
// and keeps the output to a manageable size; use the wall probes above for
// spectra and these for wall-normal structure.
//
// Points are geometrically clustered so the first sample sits inside the
// viscous sublayer while the rake still reaches past the entropy layer.

pseBLRakes
{
    type            sets;
    libs            (sampling);

    executeControl  timeStep;
    executeInterval 100;
    writeControl    timeStep;
    writeInterval   100;

    setFormat       raw;
    interpolationScheme cellPoint;

    fields          (p T U rho);

    sets
    (
"""

RAKE_SET = """        rake_x{tag}
        {{
            type        points;
            axis        distance;
            ordered     true;
            points
            (
{points}
            );
        }}
"""


def write_probes(path, wind, lee, span, rakes, rake_length):
    out = [HEADER.format(nwind=len(wind), nlee=len(lee), nspan=len(span),
                         nrake=len(rakes), npts=RAKE_POINTS,
                         rlen=rake_length, first=RAKE_FIRST * 1e6)]

    out.append(PROBE_BLOCK.format(
        name="pseWallWindward", points=fmt_pts(wind),
        comment=("Windward (lower) surface. Wall pressure is the standard "
                 "second-mode observable.")))
    out.append(PROBE_BLOCK.format(
        name="pseWallLeeward", points=fmt_pts(lee),
        comment=("Leeward (upper) surface. At -5 deg this layer is thicker, "
                 "so its second mode sits lower in frequency.")))
    out.append(PROBE_BLOCK.format(
        name="pseWallSpanwise", points=fmt_pts(span),
        comment=("Spanwise rakes, windward. These measure the spanwise "
                 "wavenumber beta; the PSE run currently assumes beta = 0.")))

    out.append(RAKE_HEADER)
    for x, pts in rakes:
        tag = ("%.3f" % x).rstrip("0").rstrip(".").replace(".", "p")
        out.append(RAKE_SET.format(tag=tag, points=fmt_pts(pts, indent=16)))
    out.append("    );\n}\n")
    out.append("\n// "
               "***************************************************"
               "************************\n")

    with open(path, "w") as fh:
        fh.write("".join(out))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--stl", help="use this STL instead of extracting one")
    ap.add_argument("--stations", type=float, nargs="+",
                    default=DEFAULT_STATIONS,
                    help="streamwise probe stations in metres")
    ap.add_argument("--rake-length", type=float, default=RAKE_LENGTH,
                    help="wall-normal rake length in metres")
    ap.add_argument("--x-tol", type=float, default=0.75,
                    help="how far in x a facet may be from a station, metres")
    ap.add_argument("--out", default=os.path.join(CASE, "system", "probesPSE"))
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    stl = args.stl
    if stl is None:
        stl = os.path.join(CASE, "pse_wall_surface.stl")
        extract_stl(CASE, stl)

    facets = read_stl(stl)
    if not facets:
        raise SystemExit("ERROR: no facets read from %s" % stl)

    wind, lee, span, rakes, missing = build(
        facets, args.stations, args.rake_length, args.x_tol)

    print("windward probes  %d" % len(wind))
    print("leeward probes   %d" % len(lee))
    print("spanwise probes  %d" % len(span))
    print("BL rakes         %d" % len(rakes))
    if missing:
        print("\nstations that could not be placed:")
        for x, why in missing:
            print("  x = %.3f m: %s" % (x, why))
        print("Widen --x-tol or drop those stations if this is unexpected.")

    if args.dry_run:
        print("\ndry run, nothing written")
        return

    write_probes(args.out, wind, lee, span, rakes, args.rake_length)
    print("\nwritten: %s" % args.out)
    print("Check it into git so the run is reproducible.")


if __name__ == "__main__":
    main()
