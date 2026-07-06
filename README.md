# Triangle Trigonometry Explorer

A Python script that computes a comprehensive set of trigonometric and geometric properties for any triangle given its three side lengths.

## Features

- Perimeter, semi-perimeter, and area (Heron's formula)
- Interior angles, triangle classification (by sides and angles)
- Altitudes, medians, and angle bisector lengths
- Inradius, circumradius, and exradii
- Euler line: circumcenter, centroid, orthocenter, nine-point center, and key distances (including Feuerbach's theorem)
- Brocard angle, Law of Tangents, and Mollweide's formula identity checks
- Coordinate-geometry model for independent verification of all computed values
- 2×2 matplotlib visualization showing the triangle with labeled sides/angles, incircle, circumcircle, nine-point circle, Euler line, cevians, and exradii circles

## Usage

```bash
python triangle_trig.py
```

Enter three side lengths when prompted. The full numeric report prints to the terminal; the visualization requires `matplotlib`:

```bash
pip install matplotlib
```
