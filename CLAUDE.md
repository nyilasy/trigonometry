# trigonometry

Single-script project: `triangle_trig.py` takes three triangle side
lengths and computes an extensive set of trigonometric/geometric
properties (angles, area, Euler line, exradii, Brocard angle, identity
checks, etc.), then renders a 2x2 matplotlib visualization. See
`README.md` for the full feature list and usage.

## Structure

- `triangle_trig.py` — the entire project; one file, no package layout.
  - `compute_all(a, b, c)` builds a results dict with every computed
    quantity, cross-checked via an independent coordinate-geometry model.
  - `print_report(...)` prints the numeric report to the terminal.
  - `visualize(...)` draws the 2x2 matplotlib figure (requires
    `matplotlib`; the script degrades gracefully to text-only output if
    it isn't installed).

## Conventions

- Keep this as a single script unless it grows enough to justify
  splitting into modules — no premature package structure.
- Any new geometric quantity should be added to `compute_all` and
  wired into both `print_report` and, where it's a distance/center
  computable from coordinates, cross-verified in the coordinate model
  the same way existing ones (OG, OH, OI, IN, etc.) already are.
- Test changes by running the script with a scalene/acute triangle
  (e.g. 6,7,9), a right triangle (3,4,5), and an obtuse triangle
  (e.g. 2,3,4) — the obtuse case is what catches sign/placement bugs in
  altitudes and the orthocenter.
