"""
author: Martin Heistermann
Chazelle Polyhedron Generator — Wavefront OBJ format
=====================================================
Constructs the Chazelle polyhedron: a cube with wedge-shaped notches
removed from alternating sides, whose inner saddle surface is a
hyperbolic paraboloid (z = x*y).

Construction (Chazelle 1984):
  - Take a unit cube subdivided into n vertical "fin" columns along X.
  - Column i (x in [i/n, (i+1)/n]) has:
      * A top notch (wedge cut from z=1) whose depth along Y grows
        linearly with i  →  top_depth(i) = MAX_DEPTH * i/(n-1)
      * A bottom notch (wedge cut from z=0) whose depth is complementary
        →  bot_depth(i) = MAX_DEPTH - top_depth(i)
  - The inner boundary between the two notched faces traces out a
    twisted saddle — the hyperbolic paraboloid z = x*y — giving the
    characteristic shape visible on Wikipedia.

Usage:
    python chazelle_polyhedron.py [n] [output.obj]

    n       — number of fin columns (default: 6, range 2–30)
    output  — filename (default: chazelle.obj)
"""

import sys
from typing import Dict, List, Tuple

# ---------------------------------------------------------------------------
# Mesh helper
# ---------------------------------------------------------------------------

class Mesh:
    def __init__(self):
        self.vertices: List[Tuple[float,float,float]] = []
        self.faces:    List[Tuple[int,...]]            = []
        self._cache:   Dict[tuple,int]                 = {}

    def v(self, x: float, y: float, z: float) -> int:
        """Add (or reuse) a vertex; return 1-based index."""
        key = (round(x, 9), round(y, 9), round(z, 9))
        if key not in self._cache:
            self.vertices.append((round(x,6), round(y,6), round(z,6)))
            self._cache[key] = len(self.vertices)
        return self._cache[key]

    def quad(self, a, b, c, d): self.faces.append((a,b,c,d))
    def tri (self, a, b, c):    self.faces.append((a,b,c))


# ---------------------------------------------------------------------------
# Chazelle polyhedron builder
# ---------------------------------------------------------------------------

def build(n: int = 6) -> Mesh:
    """
    Build the Chazelle polyhedron.

    The unit cube [0,1]^3 is sliced into n columns along X.
    Column i occupies x in [i/n, (i+1)/n].

    Inside column i, the solid cross-section at height y (0=front,1=back):

      z_ceiling(y) = min(1, y / td)        if td > 0, else 1
      z_floor(y)   = max(0, (y-ybs)/bd)   if bd > 0, else 0

    where:
      td  = top_depth  = MAX_DEPTH * i/(n-1)   (notch from top,  y-side near 0)
      bd  = bot_depth  = MAX_DEPTH*(n-1-i)/(n-1) (notch from bot, y-side near 1)
      ybs = 1 - bd  (where bottom notch starts)

    We sweep slabs along Y and emit a closed quad mesh per column.
    Adjacent columns share their X-boundary faces (emitted once each).
    """
    if n < 2:
        raise ValueError("n must be >= 2")

    m   = Mesh()
    v   = m.v
    MAX = 0.85          # max notch depth (fraction of Y span); keep < 1
    NY  = max(n * 4, 16)  # Y slices per column for smooth saddle

    for i in range(n):
        x0 = i / n
        x1 = (i + 1) / n

        td = MAX * i / (n - 1) if n > 1 else 0.0       # top notch depth
        bd = MAX * (n-1-i) / (n-1) if n > 1 else MAX    # bottom notch depth
        ybs = 1.0 - bd   # y where bottom notch starts

        def ceil_z(y: float) -> float:
            """Top surface of solid at this y."""
            if td < 1e-12: return 1.0
            return min(1.0, y / td)

        def floor_z(y: float) -> float:
            """Bottom surface of solid at this y."""
            if bd < 1e-12: return 0.0
            if y <= ybs:   return 0.0
            return min(1.0, (y - ybs) / bd)

        ys = [j / NY for j in range(NY + 1)]

        # Build ring of 4 verts at each y-slice, then stitch into quads
        rings: List[Tuple[int,int,int,int]] = []  # (bl, br, tr, tl)
        for y in ys:
            zt = ceil_z(y)
            zb = floor_z(y)
            if zt < zb: zt = zb  # collapsed (shouldn't happen with valid params)
            bl = v(x0, y, zb)
            br = v(x1, y, zb)
            tr = v(x1, y, zt)
            tl = v(x0, y, zt)
            rings.append((bl, br, tr, tl))

        # Lateral faces between consecutive rings
        for j in range(NY):
            bl0,br0,tr0,tl0 = rings[j]
            bl1,br1,tr1,tl1 = rings[j+1]

            # Bottom face (floor, normal roughly -Z)
            m.quad(bl0, br0, br1, bl1)
            # Top face    (ceiling, normal roughly +Z)
            m.quad(tl0, tl1, tr1, tr0)
            # Left wall   (x=x0, normal -X)
            m.quad(tl0, bl0, bl1, tl1)
            # Right wall  (x=x1, normal +X)
            m.quad(br0, tr0, tr1, br1)

        # Front cap (y=0)
        bl,br,tr,tl = rings[0]
        m.quad(bl, tl, tr, br)

        # Back cap (y=1)
        bl,br,tr,tl = rings[-1]
        m.quad(bl, br, tr, tl)

    return m


# ---------------------------------------------------------------------------
# OBJ writer & stats
# ---------------------------------------------------------------------------

def write_obj(mesh: Mesh, path: str, name: str = "chazelle") -> None:
    with open(path, "w") as f:
        f.write("# Chazelle polyhedron — cube with alternating wedge notches\n")
        f.write("# Inner saddle surface = hyperbolic paraboloid (z = x*y)\n")
        f.write(f"# Vertices: {len(mesh.vertices)}  Faces: {len(mesh.faces)}\n\n")
        f.write(f"o {name}\n\n")
        for x,y,z in mesh.vertices:
            f.write(f"v {x} {y} {z}\n")
        f.write("\n")
        for face in mesh.faces:
            f.write("f " + " ".join(str(i) for i in face) + "\n")
    print(f"Wrote {len(mesh.vertices)} vertices, {len(mesh.faces)} faces → {path}")


def stats(n: int, mesh: Mesh) -> None:
    xs=[w[0] for w in mesh.vertices]; ys=[w[1] for w in mesh.vertices]; zs=[w[2] for w in mesh.vertices]
    print(f"\n  Chazelle polyhedron  (n={n} fin columns)")
    print(f"  Vertices : {len(mesh.vertices)}")
    print(f"  Faces    : {len(mesh.faces)}")
    print(f"  Bbox: x=[{min(xs):.3f},{max(xs):.3f}]  y=[{min(ys):.3f},{max(ys):.3f}]  z=[{min(zs):.3f},{max(zs):.3f}]")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    n, out = 6, "chazelle.obj"
    args = sys.argv[1:]
    if args:
        try:
            n = int(args[0])
            assert 2 <= n <= 30
        except:
            print("Error: n must be 2..30"); sys.exit(1)
    if len(args) >= 2:
        out = args[1]

    print(f"\nBuilding Chazelle polyhedron (n={n} columns) …")
    mesh = build(n)
    stats(n, mesh)
    write_obj(mesh, out, name=f"chazelle_n{n}")
    print("Done. Import into Blender, MeshLab, or any 3D viewer.\n")

if __name__ == "__main__":
    main()
