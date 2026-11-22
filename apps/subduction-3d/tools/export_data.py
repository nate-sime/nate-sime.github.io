# This file is part of the Subduction-3D app, a viewer for a 3D subduction-zone run.
# Copyright (c) 2026 Nathan Sime
# SPDX-License-Identifier: MIT
"""Export the 3D subduction-zone run of Sime et al. (2023), Fig. 16, for the viewer.

Run it with ParaView's Python, which bundles VTK, numpy, h5py and matplotlib:

    pvpython tools/export_data.py <run-dir> public/data [--cache <dir>]

<run-dir> is `2a_pushforward_n100_fine_xc75.0`: 100 output steps (0..99) over
11 Myr, each with the slab-top surface (surface_data_*), and the temperature and
velocity of the whole domain (temperature_*, velocity_*). The domain moves with
the slab and is remeshed every step, so nothing here assumes a fixed mesh.

What is written, all little-endian, gzip-compressed, described by manifest.json:

  surface   the slab top on one fixed grid, so the viewer can blend between
            keyframes: radius, slab-top temperature, and a signed distance to the
            edge of the slab (the slab grows, and its edge is the zero contour).
            The grid is gnomonic -- (R0 x/z, R0 y/z) -- in which the slab top,
            dipping at most ~52 degrees, is single valued along every radius.
  slices    temperature on the planes y = -200, 0, 200 km, as in the paper figure.
  outline   the domain's feature edges, every step, as simplified segments.
  tracers   ParaView's ParticleTracer through the velocity at every step, with
            the seeds and settings of the published animation (statev2temp.pvsm).

Each extraction is cached as .npz under --cache, so the packing can be redone
without reading the 23 GB run again.
"""

import argparse
import gzip
import json
import os
import time

import numpy as np

# --- what the published animation used (statev2temp.pvsm) -------------------
N_STEPS = 100                    # output steps 0..99
MYR_TOTAL = 11.0                 # simulated time spanned by steps 0..99
KM_PER_STEP_SCALE = 5e-5 * 11e6 / 99   # u (units of 5 cm/yr) -> km per step
SEED_ORIGIN = (0.0, -300.0, 6175.0)
SEED_POINT1 = (0.0, 300.0, 6175.0)
SEED_POINT2 = (0.0, -300.0, 6325.0)
SEED_RES = (10, 6)               # 11 x 7 = 77 seeds
REINJECT_EVERY = 25              # fresh seeds at steps 0, 25, 50, 75
SLICE_Y = (-200.0, 0.0, 200.0)

# --- what the viewer gets ----------------------------------------------------
KEY_EVERY = 3                    # keyframes 0, 3, ..., 99 (34 of them)
R0 = 6371.0                      # gnomonic scale: grid units are ~km at the surface
SURFACE_H = 4.0                  # slab grid spacing, km
SDF_SCALE = 0.25                 # km per int8 unit of the slab-edge distance
SLICE_H = 2.5                    # slice texel size, km
SLICE_X = (-400.0, 500.0)
SLICE_Z = (5897.5, 6372.5)
T_RANGE = (272.87, 1575.93)      # the colour map's range in the published figure
OUTLINE_TOL = 0.25               # Douglas-Peucker tolerance for the outline, km


def keyframes():
    return list(range(0, N_STEPS, KEY_EVERY)) + ([N_STEPS - 1] if (N_STEPS - 1) % KEY_EVERY else [])


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


# =============================================================================
# Slab surface
# =============================================================================

def read_surface(run, step):
    import h5py
    with h5py.File(os.path.join(run, f"surface_data_{step:05d}.h5"), "r") as f:
        x = f["Mesh/mesh/geometry"][:]
        tri = f["Mesh/mesh/topology"][:]
        key = next(iter(f["Function/T"]))
        T = f[f"Function/T/{key}"][:, 0]
    return x, tri, T


def gnomonic(x):
    return R0 * x[:, 0] / x[:, 2], R0 * x[:, 1] / x[:, 2]


def signed_distance(inside, h):
    """Distance (km) from each cell to the edge of `inside`, positive inside."""
    ny, nx = inside.shape
    pad = np.pad(inside, 1, constant_values=False)
    nb_out = ~pad[:-2, 1:-1] | ~pad[2:, 1:-1] | ~pad[1:-1, :-2] | ~pad[1:-1, 2:]
    nb_in = pad[:-2, 1:-1] | pad[2:, 1:-1] | pad[1:-1, :-2] | pad[1:-1, 2:]
    jj, ii = np.mgrid[0:ny, 0:nx]
    P = np.stack([ii.ravel(), jj.ravel()], 1).astype(np.float64) * h
    edge_in = P[(inside & nb_out).ravel()]      # inside cells touching outside
    edge_out = P[(~inside & nb_in).ravel()]     # outside cells touching inside
    flat = inside.ravel()
    d = np.empty(len(P))
    for sel, edge, sign in ((flat, edge_out, 1.0), (~flat, edge_in, -1.0)):
        idx = np.flatnonzero(sel)
        if len(edge) == 0:
            d[idx] = sign * 1e3
            continue
        for c in range(0, len(idx), 2048):
            q = P[idx[c:c + 2048]]
            dd = np.sqrt(((q[:, None, :] - edge[None, :, :]) ** 2).sum(-1)).min(1)
            d[idx[c:c + 2048]] = sign * (dd - 0.5 * h)
    return d.reshape(ny, nx)


def dilate_fill(v, valid, passes):
    """Extend v outward from `valid` by averaging neighbours, `passes` cells."""
    v = np.where(valid, v, 0.0)
    w = valid.astype(np.float64)
    for _ in range(passes):
        vs = np.zeros_like(v)
        ws = np.zeros_like(w)
        for sl_dst, sl_src in (
            ((slice(1, None), slice(None)), (slice(None, -1), slice(None))),
            ((slice(None, -1), slice(None)), (slice(1, None), slice(None))),
            ((slice(None), slice(1, None)), (slice(None), slice(None, -1))),
            ((slice(None), slice(None, -1)), (slice(None), slice(1, None))),
        ):
            vs[sl_dst] += (v * w)[sl_src]
            ws[sl_dst] += w[sl_src]
        grow = (w == 0) & (ws > 0)
        v = np.where(grow, vs / np.maximum(ws, 1e-12), v)
        w = np.where(grow, 1.0, w)
    rest = w == 0
    if rest.any():
        v[rest] = v[~rest].mean()
    return v


def export_surface(run, cache):
    import matplotlib.tri as mtri
    keys = keyframes()
    xs = []
    lo, hi = np.full(2, np.inf), np.full(2, -np.inf)
    for s in keys:
        x, tri, T = read_surface(run, s)
        gx, gy = gnomonic(x)
        lo = np.minimum(lo, [gx.min(), gy.min()])
        hi = np.maximum(hi, [gx.max(), gy.max()])
        xs.append((x, tri, T))
    h = SURFACE_H
    x0 = np.floor(lo / h) * h - 2 * h
    n = (np.ceil((hi - x0) / h) + 3).astype(int)
    nx, ny = int(n[0]), int(n[1])
    GX, GY = np.meshgrid(x0[0] + h * np.arange(nx), x0[1] + h * np.arange(ny))
    R = np.zeros((len(keys), ny, nx), np.float32)
    TT = np.zeros_like(R)
    SD = np.zeros_like(R)
    for k, (s, (x, tri, T)) in enumerate(zip(keys, xs)):
        t0 = time.time()
        gx, gy = gnomonic(x)
        r = np.linalg.norm(x, axis=1)
        tr = mtri.Triangulation(gx, gy, tri)
        ri = mtri.LinearTriInterpolator(tr, r)(GX, GY)
        ti = mtri.LinearTriInterpolator(tr, T)(GX, GY)
        inside = ~np.ma.getmaskarray(ri)
        R[k] = dilate_fill(np.ma.getdata(ri), inside, 12)
        TT[k] = dilate_fill(np.ma.getdata(ti), inside, 12)
        SD[k] = signed_distance(inside, h)
        log(f"surface step {s:3d}: {inside.sum()} of {nx * ny} cells on the slab ({time.time() - t0:.1f}s)")
    np.savez(os.path.join(cache, "surface.npz"), keys=keys, x0=x0, h=h, R=R, T=TT, SD=SD)


# =============================================================================
# Temperature slices and the domain outline (both from the temperature mesh)
# =============================================================================

def chain_segments(P, lines, tol):
    """Join 2-point cells into polylines, simplify each, return segments (m, 2, 3)."""
    key = np.round(P / 1e-3).astype(np.int64)
    _, uid = np.unique(key, axis=0, return_inverse=True)
    uid = uid.ravel()
    pts = np.zeros((uid.max() + 1, 3))
    pts[uid] = P
    edges = {tuple(sorted((uid[a], uid[b]))) for a, b in lines if uid[a] != uid[b]}
    adj = {}
    for a, b in edges:
        adj.setdefault(a, []).append(b)
        adj.setdefault(b, []).append(a)
    used = set()
    polylines = []

    def walk(a, b):
        line = [a, b]
        used.add((min(a, b), max(a, b)))
        while len(adj[b]) == 2:
            c = adj[b][0] if adj[b][0] != line[-2] else adj[b][1]
            e = (min(b, c), max(b, c))
            if e in used:
                break
            used.add(e)
            line.append(c)
            b = c
        return line

    for a in adj:
        if len(adj[a]) != 2:
            for b in adj[a]:
                if (min(a, b), max(a, b)) not in used:
                    polylines.append(walk(a, b))
    for a, b in edges:                       # closed loops with no junctions
        if (a, b) not in used:
            polylines.append(walk(a, b))

    segs = []
    for line in polylines:
        q = pts[line]
        keep = douglas_peucker(q, tol)
        q = q[keep]
        segs.extend(zip(q[:-1], q[1:]))
    return np.array(segs, np.float32).reshape(-1, 2, 3)


def douglas_peucker(q, tol):
    keep = np.zeros(len(q), bool)
    keep[0] = keep[-1] = True
    stack = [(0, len(q) - 1)]
    while stack:
        i, j = stack.pop()
        if j <= i + 1:
            continue
        a, b = q[i], q[j]
        ab = b - a
        L = np.linalg.norm(ab)
        w = q[i + 1:j] - a
        if L < 1e-9:
            d = np.linalg.norm(w, axis=1)
        else:
            d = np.linalg.norm(np.cross(w, ab / L), axis=1)
        m = int(np.argmax(d))
        if d[m] > tol:
            keep[i + 1 + m] = True
            stack += [(i, i + 1 + m), (i + 1 + m, j)]
    return np.flatnonzero(keep)


def export_volume(run, cache):
    from paraview import simple as pv
    from vtkmodules.numpy_interface import dataset_adapter as dsa

    files = [os.path.join(run, f"temperature_{s:05d}.xdmf") for s in range(N_STEPS)]
    reader = pv.Xdmf3ReaderS(FileName=files)
    merged = pv.MergeBlocks(Input=reader)
    edges = pv.FeatureEdges(Input=pv.ExtractSurface(Input=merged), BoundaryEdges=1, FeatureEdges=1,
                            NonManifoldEdges=1, ManifoldEdges=0, FeatureAngle=30.0)
    nx = int(round((SLICE_X[1] - SLICE_X[0]) / SLICE_H)) + 1
    nz = int(round((SLICE_Z[1] - SLICE_Z[0]) / SLICE_H)) + 1
    planes = [pv.Plane(Origin=[SLICE_X[0], y, SLICE_Z[0]], Point1=[SLICE_X[1], y, SLICE_Z[0]],
                       Point2=[SLICE_X[0], y, SLICE_Z[1]], XResolution=nx - 1, YResolution=nz - 1)
              for y in SLICE_Y]
    probe = pv.ResampleWithDataset(SourceDataArrays=merged, DestinationMesh=pv.AppendDatasets(Input=planes))
    probe.PassPointArrays = 0

    keys = set(keyframes())
    slices = np.zeros((len(keys), len(SLICE_Y), nz, nx), np.float32)
    seg_list = []
    k = 0
    for s in range(N_STEPS):
        t0 = time.time()
        edges.UpdatePipeline(float(s))
        e = dsa.WrapDataObject(edges.GetClientSideObject().GetOutputDataObject(0))
        P = np.array(e.Points, np.float64)
        lines = vtk_lines(e.VTKObject)
        seg_list.append(chain_segments(P, lines, OUTLINE_TOL))
        if s in keys:
            probe.UpdatePipeline(float(s))
            o = dsa.WrapDataObject(probe.GetClientSideObject().GetOutputDataObject(0))
            T = np.array(o.PointData["T"], np.float32)
            ok = np.array(o.PointData["vtkValidPointMask"]).astype(bool)
            T[~ok] = np.nan
            slices[k] = T.reshape(len(SLICE_Y), nz, nx)
            k += 1
        log(f"volume step {s:3d}: {len(seg_list[-1])} outline segments ({time.time() - t0:.1f}s)")
    counts = np.array([len(a) for a in seg_list])
    np.savez(os.path.join(cache, "volume.npz"), slices=slices, nx=nx, nz=nz,
             outline=np.concatenate(seg_list), outline_counts=counts)


def vtk_lines(poly):
    from vtkmodules.util.numpy_support import vtk_to_numpy
    ca = poly.GetLines()
    off = vtk_to_numpy(ca.GetOffsetsArray())
    con = vtk_to_numpy(ca.GetConnectivityArray())
    segs = []
    for i in range(len(off) - 1):
        ids = con[off[i]:off[i + 1]]
        segs.extend(zip(ids[:-1], ids[1:]))
    return np.array(segs, np.int64).reshape(-1, 2)


# =============================================================================
# Tracers
# =============================================================================

def export_tracers(run, cache):
    from paraview import simple as pv
    from vtkmodules.numpy_interface import dataset_adapter as dsa

    files = [os.path.join(run, f"velocity_{s:05d}.xdmf") for s in range(N_STEPS)]
    reader = pv.Xdmf3ReaderS(FileName=files)
    calc = pv.Calculator(Input=pv.MergeBlocks(Input=reader), Function=f"{KM_PER_STEP_SCALE!r} * u",
                         ResultArrayName="Result")
    seeds = pv.Plane(Origin=list(SEED_ORIGIN), Point1=list(SEED_POINT1), Point2=list(SEED_POINT2),
                     XResolution=SEED_RES[0], YResolution=SEED_RES[1])
    tracer = pv.ParticleTracer(Input=calc, SeedSource=seeds)
    tracer.SelectInputVectors = ["POINTS", "Result"]
    tracer.ForceReinjectionEveryNSteps = REINJECT_EVERY
    tracer.StaticSeeds = 0
    tracer.MeshOverTime = "Different"
    tracer.InterpolatorType = "Interpolator with Cell Locator"
    tracer.ComputeVorticity = 0

    # ParticleId restarts at 0 with every reinjection: a tracer is (injection step, id).
    track = {}
    for s in range(N_STEPS):
        t0 = time.time()
        tracer.UpdatePipeline(float(s))
        o = dsa.WrapDataObject(tracer.GetClientSideObject().GetOutputDataObject(0))
        if o.GetNumberOfPoints():
            ids = np.array(o.PointData["ParticleId"]).astype(np.int64)
            inj = np.array(o.PointData["InjectionStepId"]).astype(np.int64)
            P = np.array(o.Points, np.float64)
            for key, p in zip(zip(inj.tolist(), ids.tolist()), P):
                track.setdefault(key, np.full((N_STEPS, 3), np.nan))[s] = p
        log(f"tracer step {s:3d}: {o.GetNumberOfPoints()} alive ({time.time() - t0:.1f}s)")
    keys = sorted(track)
    np.savez(os.path.join(cache, "tracers.npz"), ids=np.array(keys), pos=np.stack([track[k] for k in keys]))


# =============================================================================
# Packing for the viewer
# =============================================================================

def write_gz(out, name, arr):
    data = gzip.compress(np.ascontiguousarray(arr).tobytes(), compresslevel=9, mtime=0)
    with open(os.path.join(out, name), "wb") as f:
        f.write(data)
    return {"file": name, "dtype": str(arr.dtype), "shape": list(arr.shape), "bytes": len(data)}


def delta_over_frames(q):
    """Store frame k as frame k minus frame k-1, wrapping: the frames change slowly."""
    u = q.view({1: np.uint8, 2: np.uint16}[q.dtype.itemsize])
    d = u.copy()
    d[1:] = u[1:] - u[:-1]
    return d.view(q.dtype)  # the manifest records the dtype the frames decode to


def pack(run, cache, out):
    os.makedirs(out, exist_ok=True)
    keys = keyframes()
    man = {
        "source": "Sime et al. (2023), arXiv:2309.12895, Fig. 16: run " + os.path.basename(os.path.normpath(run)),
        "steps": N_STEPS,
        "myrPerStep": MYR_TOTAL / (N_STEPS - 1),
        "keyframes": keys,
        "temperatureRange": list(T_RANGE),
        "encoding": "Each file is gzip of a little-endian array. 'delta' arrays store frame k "
                    "minus frame k-1 (wrapping): undo with a running sum over frames.",
    }

    s = np.load(os.path.join(cache, "surface.npz"))
    R, T, SD = s["R"], s["T"], s["SD"]
    r_lo, r_hi = float(np.floor(R.min())), float(np.ceil(R.max()))
    qR = np.round((R - r_lo) / (r_hi - r_lo) * 65535).astype(np.uint16)
    qT = np.round(np.clip((T - T_RANGE[0]) / (T_RANGE[1] - T_RANGE[0]), 0, 1) * 255).astype(np.uint8)
    qS = np.clip(np.round(SD / SDF_SCALE), -127, 127).astype(np.int8)
    man["surface"] = {
        "nx": int(R.shape[2]), "ny": int(R.shape[1]),
        "origin": [float(v) for v in s["x0"]], "spacing": float(s["h"]), "R0": R0,
        "radius": dict(write_gz(out, "surface-radius.bin.gz", delta_over_frames(qR)), range=[r_lo, r_hi], delta=True),
        "temperature": dict(write_gz(out, "surface-temperature.bin.gz", delta_over_frames(qT)), delta=True),
        "edgeDistance": dict(write_gz(out, "surface-edge.bin.gz", delta_over_frames(qS)), kmPerUnit=SDF_SCALE, delta=True),
        "about": "Grid node (i, j) is the radial line through (R0 x/z, R0 y/z) = origin + spacing (i, j). "
                 "radius is quantised over range to 0..65535, temperature over temperatureRange to 0..255, "
                 "edgeDistance is the signed distance to the slab's edge (positive on the slab).",
    }

    v = np.load(os.path.join(cache, "volume.npz"))
    sl = v["slices"]
    q = np.where(np.isnan(sl), 0, 1 + np.round(np.clip((sl - T_RANGE[0]) / (T_RANGE[1] - T_RANGE[0]), 0, 1) * 254))
    man["slices"] = dict(
        write_gz(out, "slices.bin.gz", delta_over_frames(q.astype(np.uint8))), delta=True,
        y=list(SLICE_Y), nx=int(v["nx"]), nz=int(v["nz"]), x=list(SLICE_X), z=list(SLICE_Z),
        about="Texel (i, j) of slice m is the point (x0 + (x1 - x0) i / (nx - 1), y[m], z0 + (z1 - z0) j / (nz - 1)); "
              "0 is outside the domain, 1..255 spans temperatureRange.")

    man["outline"] = dict(write_gz(out, "outline.bin.gz", v["outline"].astype(np.float32)),
                          counts=[int(c) for c in v["outline_counts"]],
                          about="Segments (2 points of x, y, z in km) of the domain's feature edges; "
                                "counts[step] of them per step, in step order.")

    t = np.load(os.path.join(cache, "tracers.npz"))
    pos = t["pos"].astype(np.float32)
    man["tracers"] = dict(write_gz(out, "tracers.bin.gz", pos),
                          about="Position (km) of each tracer at every step; NaN before it is seeded "
                                "and after it leaves the domain.",
                          seeds={"origin": SEED_ORIGIN, "point1": SEED_POINT1, "point2": SEED_POINT2,
                                 "resolution": SEED_RES, "reinjectEvery": REINJECT_EVERY})

    man["camera"] = {  # the published animation's view
        "position": [1419.0518151872818, -680.7404666277585, 6788.765304396278],
        "target": [-587.6020777552321, 426.0370709601656, 5627.056511118317],
        "up": [-0.3583189295605358, 0.28388097729773726, 0.8893925654327641],
        "fov": 30.0,
    }
    with open(os.path.join(out, "manifest.json"), "w") as f:
        json.dump(man, f, indent=1)
    total = sum(os.path.getsize(os.path.join(out, n)) for n in os.listdir(out))
    log(f"wrote {out}: {total / 1e6:.2f} MB")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run")
    ap.add_argument("out")
    ap.add_argument("--cache", default=None, help="where extractions are kept (default: tools/.export-cache)")
    ap.add_argument("--only", default="surface,volume,tracers,pack",
                    help="comma-separated stages: surface, volume, tracers, pack")
    a = ap.parse_args()
    cache = a.cache or os.path.join(os.path.dirname(os.path.abspath(__file__)), ".export-cache")
    os.makedirs(cache, exist_ok=True)
    stages = a.only.split(",")
    if "surface" in stages:
        export_surface(a.run, cache)
    if "volume" in stages:
        export_volume(a.run, cache)
    if "tracers" in stages:
        export_tracers(a.run, cache)
    if "pack" in stages:
        pack(a.run, cache, a.out)


if __name__ == "__main__":
    main()
