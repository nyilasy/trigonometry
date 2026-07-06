"""
Triangle Trigonometry Explorer
==============================

Given the three side lengths of a triangle, this script computes an
extensive set of trigonometric and geometric characteristics:

  - Perimeter / semi-perimeter
  - Area (Heron's formula, and cross-checks via inradius/circumradius)
  - Interior angles (degrees and radians)
  - Triangle classification (by sides and by angles)
  - Altitudes, medians, angle bisector lengths
  - Inradius, circumradius, exradii (+ their circles' areas)
  - Distances from vertices to the incenter, incircle tangent lengths
  - Euler line: circumcenter O, centroid G, orthocenter H, nine-point
    center N, and the classic distances OG, GH, OH, OI, IN (Feuerbach)
  - Brocard angle
  - Law of Tangents and Mollweide's formula identity checks
  - A coordinate-geometry model of the triangle used to independently
    verify every center/distance computed from closed-form formulas
  - A 2x2 matplotlib visualization of the triangle: labeled sides/angles,
    incircle/circumcircle/nine-point circle with the Euler line, the
    cevians (medians/altitudes/bisectors), and the exradii circles

Run it and enter three side lengths when prompted. Visualization requires
matplotlib (`pip install matplotlib`); the script still prints the full
numeric report without it.
"""

import math

try:
    import matplotlib.pyplot as plt
    from matplotlib.patches import Circle
    HAVE_MPL = True
except ImportError:
    HAVE_MPL = False


# --------------------------------------------------------------------------
# Input handling
# --------------------------------------------------------------------------

def get_positive_float(prompt):
    while True:
        raw = input(prompt).strip()
        try:
            value = float(raw)
        except ValueError:
            print("  -> Please enter a valid number.")
            continue
        if value <= 0:
            print("  -> Length must be a positive number.")
            continue
        return value


def input_triangle_sides():
    print("Enter the lengths of the three sides of the triangle.")
    while True:
        a = get_positive_float("Side a: ")
        b = get_positive_float("Side b: ")
        c = get_positive_float("Side c: ")
        if a + b > c and b + c > a and a + c > b:
            return a, b, c
        print("  -> These lengths do not satisfy the triangle inequality "
              "(sum of any two sides must exceed the third). Try again.\n")


# --------------------------------------------------------------------------
# Core trigonometric computations
# --------------------------------------------------------------------------

def law_of_cosines_angles(a, b, c):
    """Return interior angles A, B, C (radians), opposite sides a, b, c."""
    cos_a = (b ** 2 + c ** 2 - a ** 2) / (2 * b * c)
    cos_b = (a ** 2 + c ** 2 - b ** 2) / (2 * a * c)
    cos_c = (a ** 2 + b ** 2 - c ** 2) / (2 * a * b)
    # Clamp for floating-point safety before acos
    cos_a = max(-1.0, min(1.0, cos_a))
    cos_b = max(-1.0, min(1.0, cos_b))
    cos_c = max(-1.0, min(1.0, cos_c))
    return math.acos(cos_a), math.acos(cos_b), math.acos(cos_c)


def classify_by_sides(a, b, c):
    if math.isclose(a, b) and math.isclose(b, c):
        return "Equilateral"
    if math.isclose(a, b) or math.isclose(b, c) or math.isclose(a, c):
        return "Isosceles"
    return "Scalene"


def classify_by_angles(A, B, C):
    """A, B, C in radians."""
    right = math.pi / 2
    angles = (A, B, C)
    if any(math.isclose(ang, right, abs_tol=1e-9) for ang in angles):
        return "Right"
    if any(ang > right for ang in angles):
        return "Obtuse"
    return "Acute"


# --------------------------------------------------------------------------
# Coordinate-geometry model (used to verify formula-based centers)
# --------------------------------------------------------------------------

def build_coordinates(a, b, c, A):
    """Place the triangle in the plane: A=(0,0), B=(c,0), C from angle A."""
    Ax, Ay = 0.0, 0.0
    Bx, By = c, 0.0
    Cx, Cy = b * math.cos(A), b * math.sin(A)
    return (Ax, Ay), (Bx, By), (Cx, Cy)


def circumcenter(P1, P2, P3):
    (x1, y1), (x2, y2), (x3, y3) = P1, P2, P3
    d = 2 * (x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2))
    ux = ((x1 ** 2 + y1 ** 2) * (y2 - y3) +
          (x2 ** 2 + y2 ** 2) * (y3 - y1) +
          (x3 ** 2 + y3 ** 2) * (y1 - y2)) / d
    uy = ((x1 ** 2 + y1 ** 2) * (x3 - x2) +
          (x2 ** 2 + y2 ** 2) * (x1 - x3) +
          (x3 ** 2 + y3 ** 2) * (x2 - x1)) / d
    return ux, uy


def dist(P, Q):
    return math.hypot(P[0] - Q[0], P[1] - Q[1])


def midpoint(P, Q):
    return ((P[0] + Q[0]) / 2, (P[1] + Q[1]) / 2)


# --------------------------------------------------------------------------
# Master computation
# --------------------------------------------------------------------------

def compute_all(a, b, c):
    A, B, C = law_of_cosines_angles(a, b, c)

    perimeter = a + b + c
    s = perimeter / 2
    area = math.sqrt(s * (s - a) * (s - b) * (s - c))  # Heron's formula

    # Altitudes (from each vertex to the opposite side)
    h_a = 2 * area / a
    h_b = 2 * area / b
    h_c = 2 * area / c

    # Medians (from each vertex to the midpoint of the opposite side)
    m_a = 0.5 * math.sqrt(2 * b ** 2 + 2 * c ** 2 - a ** 2)
    m_b = 0.5 * math.sqrt(2 * a ** 2 + 2 * c ** 2 - b ** 2)
    m_c = 0.5 * math.sqrt(2 * a ** 2 + 2 * b ** 2 - c ** 2)

    # Internal angle bisector lengths (vertex to opposite side)
    t_a = math.sqrt(b * c * ((b + c) ** 2 - a ** 2)) / (b + c)
    t_b = math.sqrt(a * c * ((a + c) ** 2 - b ** 2)) / (a + c)
    t_c = math.sqrt(a * b * ((a + b) ** 2 - c ** 2)) / (a + b)

    # Inradius / circumradius
    r = area / s
    R = (a * b * c) / (4 * area)

    # Exradii
    r_a = area / (s - a)
    r_b = area / (s - b)
    r_c = area / (s - c)

    # Distances from vertices to the incenter
    AI = r / math.sin(A / 2)
    BI = r / math.sin(B / 2)
    CI = r / math.sin(C / 2)

    # Incircle tangent lengths from each vertex (touch-point distances)
    tangent_a = s - a  # from vertex A along sides AB/AC
    tangent_b = s - b  # from vertex B
    tangent_c = s - c  # from vertex C

    # Euler line quantities (formula-based)
    OI = math.sqrt(max(0.0, R ** 2 - 2 * R * r))                # Euler's formula
    OH_sq = 9 * R ** 2 - (a ** 2 + b ** 2 + c ** 2)
    OH = math.sqrt(max(0.0, OH_sq))
    OG = OH / 3
    GH = 2 * OH / 3
    nine_point_radius = R / 2
    IN_feuerbach = nine_point_radius - r                          # Feuerbach's theorem

    # Brocard angle: cot(w) = cot A + cot B + cot C
    cot_sum = (1 / math.tan(A)) + (1 / math.tan(B)) + (1 / math.tan(C))
    brocard = math.atan(1 / cot_sum) if cot_sum != 0 else math.pi / 2

    # Area cross-checks
    area_via_r = r * s
    area_via_R = (a * b * c) / (4 * R)
    area_via_exradii = math.sqrt(r_a * r_b * r_c * r)

    # Law of Tangents identity check: (a-b)/(a+b) == tan((A-B)/2) / tan((A+B)/2)
    lot_lhs = (a - b) / (a + b)
    lot_rhs = math.tan((A - B) / 2) / math.tan((A + B) / 2)

    # Mollweide's formula identity check: (a+b)/c == cos((A-B)/2) / sin(C/2)
    moll_lhs = (a + b) / c
    moll_rhs = math.cos((A - B) / 2) / math.sin(C / 2)

    # Law of Sines constant 2R, verified three ways
    two_r_a = a / math.sin(A)
    two_r_b = b / math.sin(B)
    two_r_c = c / math.sin(C)

    # Orthic triangle sides (only geometrically meaningful for acute triangles)
    orthic_a = a * abs(math.cos(A))
    orthic_b = b * abs(math.cos(B))
    orthic_c = c * abs(math.cos(C))

    # --- Coordinate-geometry cross-verification ---
    P_A, P_B, P_C = build_coordinates(a, b, c, A)
    O_pt = circumcenter(P_A, P_B, P_C)
    G_pt = ((P_A[0] + P_B[0] + P_C[0]) / 3, (P_A[1] + P_B[1] + P_C[1]) / 3)
    # H = A + B + C - 2*O   (since centroid divides OH in ratio 1:2)
    H_pt = (P_A[0] + P_B[0] + P_C[0] - 2 * O_pt[0],
            P_A[1] + P_B[1] + P_C[1] - 2 * O_pt[1])
    perim = a + b + c
    I_pt = ((a * P_A[0] + b * P_B[0] + c * P_C[0]) / perim,
            (a * P_A[1] + b * P_B[1] + c * P_C[1]) / perim)
    N_pt = midpoint(O_pt, H_pt)

    coord_checks = {
        "OG (coords)": dist(O_pt, G_pt),
        "OH (coords)": dist(O_pt, H_pt),
        "GH (coords)": dist(G_pt, H_pt),
        "OI (coords)": dist(O_pt, I_pt),
        "IN (coords)": dist(I_pt, N_pt),
        "Area (shoelace, coords)": 0.5 * abs(
            P_A[0] * (P_B[1] - P_C[1]) +
            P_B[0] * (P_C[1] - P_A[1]) +
            P_C[0] * (P_A[1] - P_B[1])
        ),
    }

    return {
        "angles_rad": (A, B, C),
        "angles_deg": (math.degrees(A), math.degrees(B), math.degrees(C)),
        "perimeter": perimeter,
        "semiperimeter": s,
        "area": area,
        "class_sides": classify_by_sides(a, b, c),
        "class_angles": classify_by_angles(A, B, C),
        "altitudes": (h_a, h_b, h_c),
        "medians": (m_a, m_b, m_c),
        "bisectors": (t_a, t_b, t_c),
        "inradius": r,
        "circumradius": R,
        "exradii": (r_a, r_b, r_c),
        "incircle_area": math.pi * r ** 2,
        "circumcircle_area": math.pi * R ** 2,
        "vertex_to_incenter": (AI, BI, CI),
        "incircle_tangent_lengths": (tangent_a, tangent_b, tangent_c),
        "euler": {
            "OI": OI, "OH": OH, "OG": OG, "GH": GH,
            "nine_point_radius": nine_point_radius,
            "IN_feuerbach": IN_feuerbach,
        },
        "brocard_angle_rad": brocard,
        "brocard_angle_deg": math.degrees(brocard),
        "area_checks": (area_via_r, area_via_R, area_via_exradii),
        "law_of_tangents": (lot_lhs, lot_rhs),
        "mollweide": (moll_lhs, moll_rhs),
        "law_of_sines_2R": (two_r_a, two_r_b, two_r_c),
        "orthic_sides": (orthic_a, orthic_b, orthic_c),
        "coords": {"A": P_A, "B": P_B, "C": P_C, "O": O_pt, "G": G_pt,
                   "H": H_pt, "I": I_pt, "N": N_pt},
        "coord_checks": coord_checks,
    }


# --------------------------------------------------------------------------
# Reporting
# --------------------------------------------------------------------------

def fmt(x, nd=6):
    return f"{x:.{nd}f}"


def print_report(a, b, c, res):
    A, B, C = res["angles_rad"]
    Ad, Bd, Cd = res["angles_deg"]

    print("\n" + "=" * 70)
    print(f" TRIANGLE  a={a}, b={b}, c={c}")
    print("=" * 70)

    print(f"\nClassification : {res['class_sides']} / {res['class_angles']}")

    print("\n--- Basic measures ---")
    print(f"Perimeter          : {fmt(res['perimeter'])}")
    print(f"Semi-perimeter (s) : {fmt(res['semiperimeter'])}")
    print(f"Area (Heron)       : {fmt(res['area'])}")

    print("\n--- Interior angles ---")
    print(f"A = {fmt(Ad)} deg  = {fmt(A)} rad  (opposite side a)")
    print(f"B = {fmt(Bd)} deg  = {fmt(B)} rad  (opposite side b)")
    print(f"C = {fmt(Cd)} deg  = {fmt(C)} rad  (opposite side c)")
    print(f"Sum check          : {fmt(Ad + Bd + Cd)} deg (should be 180)")

    h_a, h_b, h_c = res["altitudes"]
    print("\n--- Altitudes ---")
    print(f"h_a = {fmt(h_a)}   h_b = {fmt(h_b)}   h_c = {fmt(h_c)}")

    m_a, m_b, m_c = res["medians"]
    print("\n--- Medians ---")
    print(f"m_a = {fmt(m_a)}   m_b = {fmt(m_b)}   m_c = {fmt(m_c)}")

    t_a, t_b, t_c = res["bisectors"]
    print("\n--- Internal angle bisector lengths ---")
    print(f"t_a = {fmt(t_a)}   t_b = {fmt(t_b)}   t_c = {fmt(t_c)}")

    print("\n--- Inscribed / circumscribed circles ---")
    print(f"Inradius  r        : {fmt(res['inradius'])}")
    print(f"Circumradius R     : {fmt(res['circumradius'])}")
    print(f"Incircle area      : {fmt(res['incircle_area'])}")
    print(f"Circumcircle area  : {fmt(res['circumcircle_area'])}")

    r_a, r_b, r_c = res["exradii"]
    print("\n--- Exradii (escribed circles) ---")
    print(f"r_a = {fmt(r_a)}   r_b = {fmt(r_b)}   r_c = {fmt(r_c)}")

    AI, BI, CI = res["vertex_to_incenter"]
    print("\n--- Vertex-to-incenter distances ---")
    print(f"AI = {fmt(AI)}   BI = {fmt(BI)}   CI = {fmt(CI)}")

    ta, tb, tc = res["incircle_tangent_lengths"]
    print("\n--- Incircle tangent lengths from each vertex ---")
    print(f"from A = {fmt(ta)} (=s-a)   from B = {fmt(tb)} (=s-b)   from C = {fmt(tc)} (=s-c)")

    eu = res["euler"]
    print("\n--- Euler line & nine-point circle ---")
    print(f"OG (circumcenter-centroid)   : {fmt(eu['OG'])}")
    print(f"GH (centroid-orthocenter)    : {fmt(eu['GH'])}")
    print(f"OH (circumcenter-orthocenter): {fmt(eu['OH'])}")
    print(f"OI (circumcenter-incenter, Euler's formula) : {fmt(eu['OI'])}")
    print(f"Nine-point circle radius     : {fmt(eu['nine_point_radius'])}")
    print(f"IN (incenter to nine-point center, Feuerbach): {fmt(eu['IN_feuerbach'])}")

    print("\n--- Brocard angle ---")
    print(f"omega = {fmt(res['brocard_angle_deg'])} deg = {fmt(res['brocard_angle_rad'])} rad")

    av_r, av_R, av_ex = res["area_checks"]
    print("\n--- Area cross-checks (should all match Heron's area) ---")
    print(f"Area = r*s              : {fmt(av_r)}")
    print(f"Area = abc/(4R)         : {fmt(av_R)}")
    print(f"Area = sqrt(r*ra*rb*rc) : {fmt(av_ex)}")

    lot_l, lot_r = res["law_of_tangents"]
    print("\n--- Law of Tangents identity check ---")
    print(f"(a-b)/(a+b)              = {fmt(lot_l)}")
    print(f"tan((A-B)/2)/tan((A+B)/2) = {fmt(lot_r)}")

    mo_l, mo_r = res["mollweide"]
    print("\n--- Mollweide's formula identity check ---")
    print(f"(a+b)/c                   = {fmt(mo_l)}")
    print(f"cos((A-B)/2)/sin(C/2)     = {fmt(mo_r)}")

    twoRa, twoRb, twoRc = res["law_of_sines_2R"]
    print("\n--- Law of Sines constant (2R), verified per side ---")
    print(f"a/sinA = {fmt(twoRa)}   b/sinB = {fmt(twoRb)}   c/sinC = {fmt(twoRc)}")

    oa, ob, oc = res["orthic_sides"]
    print("\n--- Orthic triangle side lengths (meaningful for acute triangles) ---")
    print(f"a' = a*|cosA| = {fmt(oa)}   b' = b*|cosB| = {fmt(ob)}   c' = c*|cosC| = {fmt(oc)}")

    print("\n--- Coordinate-geometry model (independent verification) ---")
    coords = res["coords"]
    for label in ("A", "B", "C", "O", "G", "H", "I", "N"):
        x, y = coords[label]
        print(f"{label} = ({fmt(x)}, {fmt(y)})")

    cc = res["coord_checks"]
    print("\nDistances recomputed directly from coordinates (cross-check vs. formulas above):")
    for k, v in cc.items():
        print(f"{k:28s}: {fmt(v)}")

    print("\n" + "=" * 70 + "\n")


# --------------------------------------------------------------------------
# Visualization helpers
# --------------------------------------------------------------------------

def foot_perpendicular(P, Q1, Q2):
    """Foot of the perpendicular dropped from P onto line Q1-Q2."""
    dx, dy = Q2[0] - Q1[0], Q2[1] - Q1[1]
    t = ((P[0] - Q1[0]) * dx + (P[1] - Q1[1]) * dy) / (dx ** 2 + dy ** 2)
    return (Q1[0] + t * dx, Q1[1] + t * dy)


def point_on_segment(P, Q, t):
    return (P[0] + t * (Q[0] - P[0]), P[1] + t * (Q[1] - P[1]))


def unit(vx, vy):
    n = math.hypot(vx, vy)
    return (vx / n, vy / n)


def set_bounds(ax, points, radii=None, margin_ratio=0.25):
    xs = [p[0] for p in points]
    ys = [p[1] for p in points]
    if radii:
        for (cx, cy), r in radii:
            xs += [cx - r, cx + r]
            ys += [cy - r, cy + r]
    xmin, xmax, ymin, ymax = min(xs), max(xs), min(ys), max(ys)
    span = max(xmax - xmin, ymax - ymin, 1e-9)
    margin = span * margin_ratio
    cx, cy = (xmin + xmax) / 2, (ymin + ymax) / 2
    half = span / 2 + margin
    ax.set_xlim(cx - half, cx + half)
    ax.set_ylim(cy - half, cy + half)
    ax.set_aspect("equal", adjustable="box")


def draw_triangle_outline(ax, P_A, P_B, P_C, **kwargs):
    xs = [P_A[0], P_B[0], P_C[0], P_A[0]]
    ys = [P_A[1], P_B[1], P_C[1], P_A[1]]
    ax.plot(xs, ys, color=kwargs.get("color", "black"),
            linewidth=kwargs.get("linewidth", 1.5), zorder=2)


def label_vertices(ax, P_A, P_B, P_C, centroid):
    for label, P in (("A", P_A), ("B", P_B), ("C", P_C)):
        dx, dy = unit(P[0] - centroid[0], P[1] - centroid[1])
        ax.annotate(label, (P[0] + dx * 0.28, P[1] + dy * 0.28),
                    ha="center", va="center", fontsize=12, fontweight="bold")
        ax.plot(*P, "o", color="black", markersize=4, zorder=3)


# --------------------------------------------------------------------------
# Main visualization: 2x2 figure
# --------------------------------------------------------------------------

def visualize(a, b, c, res):
    coords = res["coords"]
    P_A, P_B, P_C = coords["A"], coords["B"], coords["C"]
    O_pt, G_pt, H_pt = coords["O"], coords["G"], coords["H"]
    I_pt, N_pt = coords["I"], coords["N"]
    centroid = G_pt
    A, B, C = res["angles_rad"]
    Ad, Bd, Cd = res["angles_deg"]

    fig, axes = plt.subplots(2, 2, figsize=(13, 13))
    fig.suptitle(f"Triangle a={a}, b={b}, c={c}  "
                 f"({res['class_sides']} / {res['class_angles']})",
                 fontsize=14, fontweight="bold")

    # ---- Panel 1: sides + angles ----
    ax = axes[0, 0]
    ax.set_title("Sides and interior angles")
    draw_triangle_outline(ax, P_A, P_B, P_C)
    label_vertices(ax, P_A, P_B, P_C, centroid)

    for (P1, P2, length, name) in (
        (P_B, P_C, a, "a"), (P_A, P_C, b, "b"), (P_A, P_B, c, "c")
    ):
        mid = ((P1[0] + P2[0]) / 2, (P1[1] + P2[1]) / 2)
        dx, dy = unit(mid[0] - centroid[0], mid[1] - centroid[1])
        ax.annotate(f"{name}={length:.2f}",
                    (mid[0] + dx * 0.22, mid[1] + dy * 0.22),
                    ha="center", va="center", fontsize=9, color="darkblue")

    arc_r = 0.15 * min(a, b, c)
    for vertex, P1, P2, angle_deg in (
        (P_A, P_B, P_C, Ad), (P_B, P_A, P_C, Bd), (P_C, P_A, P_B, Cd)
    ):
        d1 = unit(P1[0] - vertex[0], P1[1] - vertex[1])
        d2 = unit(P2[0] - vertex[0], P2[1] - vertex[1])
        bis = unit(d1[0] + d2[0], d1[1] + d2[1])
        t = [i / 20 for i in range(21)]
        ang1 = math.atan2(d1[1], d1[0])
        ang2 = math.atan2(d2[1], d2[0])
        # shortest angular path from ang1 to ang2
        diff = (ang2 - ang1 + math.pi) % (2 * math.pi) - math.pi
        arc_pts = [(vertex[0] + arc_r * math.cos(ang1 + diff * s),
                    vertex[1] + arc_r * math.sin(ang1 + diff * s)) for s in t]
        ax.plot([p[0] for p in arc_pts], [p[1] for p in arc_pts],
                color="crimson", linewidth=1.3, zorder=2)
        ax.annotate(f"{angle_deg:.1f} deg",
                    (vertex[0] + bis[0] * arc_r * 2.3,
                     vertex[1] + bis[1] * arc_r * 2.3),
                    ha="center", va="center", fontsize=8.5, color="crimson")

    set_bounds(ax, [P_A, P_B, P_C])
    ax.axis("off")

    # ---- Panel 2: circles + centers + Euler line ----
    ax = axes[0, 1]
    ax.set_title("Incircle, circumcircle, nine-point circle & Euler line")
    draw_triangle_outline(ax, P_A, P_B, P_C, color="gray", linewidth=1)

    R = res["circumradius"]
    r = res["inradius"]
    npr = res["euler"]["nine_point_radius"]
    ax.add_patch(Circle(O_pt, R, fill=False, color="tab:blue", linewidth=1.6,
                         label=f"Circumcircle (R={R:.2f})"))
    ax.add_patch(Circle(I_pt, r, fill=False, color="tab:green", linewidth=1.6,
                         label=f"Incircle (r={r:.2f})"))
    ax.add_patch(Circle(N_pt, npr, fill=False, color="tab:orange", linewidth=1.6,
                         linestyle="--", label=f"Nine-point circle (R/2={npr:.2f})"))

    euler_pts = sorted([O_pt, G_pt, H_pt], key=lambda p: p[0])
    ax.plot([p[0] for p in euler_pts], [p[1] for p in euler_pts],
            color="black", linewidth=1, linestyle=":", zorder=1, label="Euler line")

    for label, P, color in (("O", O_pt, "tab:blue"), ("G", G_pt, "black"),
                             ("H", H_pt, "tab:red"), ("I", I_pt, "tab:green"),
                             ("N", N_pt, "tab:orange")):
        ax.plot(*P, "o", color=color, markersize=6, zorder=4)
        ax.annotate(label, (P[0], P[1]), textcoords="offset points",
                    xytext=(6, 6), fontsize=10, fontweight="bold", color=color)

    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.02), fontsize=8, ncol=2)
    set_bounds(ax, [P_A, P_B, P_C, O_pt, G_pt, H_pt, I_pt, N_pt],
               radii=[(O_pt, R), (I_pt, r), (N_pt, npr)])
    ax.axis("off")

    # ---- Panel 3: medians, altitudes, angle bisectors ----
    ax = axes[1, 0]
    ax.set_title("Cevians: medians, altitudes, angle bisectors")
    draw_triangle_outline(ax, P_A, P_B, P_C, color="gray", linewidth=1)
    label_vertices(ax, P_A, P_B, P_C, centroid)

    # Medians: vertex to midpoint of opposite side
    medians = (
        (P_A, ((P_B[0] + P_C[0]) / 2, (P_B[1] + P_C[1]) / 2)),
        (P_B, ((P_A[0] + P_C[0]) / 2, (P_A[1] + P_C[1]) / 2)),
        (P_C, ((P_A[0] + P_B[0]) / 2, (P_A[1] + P_B[1]) / 2)),
    )
    for i, (V, foot) in enumerate(medians):
        ax.plot([V[0], foot[0]], [V[1], foot[1]], color="tab:purple",
                linewidth=1.3, label="Median" if i == 0 else None)

    # Altitudes: vertex to foot of perpendicular on opposite side
    altitudes = (
        (P_A, foot_perpendicular(P_A, P_B, P_C)),
        (P_B, foot_perpendicular(P_B, P_A, P_C)),
        (P_C, foot_perpendicular(P_C, P_A, P_B)),
    )
    for i, (V, foot) in enumerate(altitudes):
        ax.plot([V[0], foot[0]], [V[1], foot[1]], color="tab:red",
                linewidth=1.3, linestyle="--", label="Altitude" if i == 0 else None)

    # Angle bisectors: vertex to foot on opposite side (angle bisector theorem)
    bis_a = point_on_segment(P_B, P_C, c / (b + c))
    bis_b = point_on_segment(P_A, P_C, c / (a + c))
    bis_c = point_on_segment(P_A, P_B, b / (a + b))
    for i, (V, foot) in enumerate(((P_A, bis_a), (P_B, bis_b), (P_C, bis_c))):
        ax.plot([V[0], foot[0]], [V[1], foot[1]], color="tab:green",
                linewidth=1.3, linestyle=":", label="Angle bisector" if i == 0 else None)

    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.02), fontsize=9, ncol=3)
    set_bounds(ax, [P_A, P_B, P_C])
    ax.axis("off")

    # ---- Panel 4: exradii (escribed circles) ----
    ax = axes[1, 1]
    ax.set_title("Exradii (escribed circles) + incircle")
    draw_triangle_outline(ax, P_A, P_B, P_C, color="gray", linewidth=1)
    label_vertices(ax, P_A, P_B, P_C, centroid)

    r_a, r_b, r_c = res["exradii"]
    perim = a + b + c
    # Excenters via barycentric weights (-a:b:c), (a:-b:c), (a:b:-c)
    I_A = ((-a * P_A[0] + b * P_B[0] + c * P_C[0]) / (-a + b + c),
           (-a * P_A[1] + b * P_B[1] + c * P_C[1]) / (-a + b + c))
    I_B = ((a * P_A[0] - b * P_B[0] + c * P_C[0]) / (a - b + c),
           (a * P_A[1] - b * P_B[1] + c * P_C[1]) / (a - b + c))
    I_C = ((a * P_A[0] + b * P_B[0] - c * P_C[0]) / (a + b - c),
           (a * P_A[1] + b * P_B[1] - c * P_C[1]) / (a + b - c))

    ax.add_patch(Circle(I_pt, r, fill=False, color="tab:green", linewidth=1.6,
                         label=f"Incircle (r={r:.2f})"))
    for label, center, radius, color in (
        ("I_A", I_A, r_a, "tab:red"),
        ("I_B", I_B, r_b, "tab:orange"),
        ("I_C", I_C, r_c, "tab:brown"),
    ):
        ax.add_patch(Circle(center, radius, fill=False, color=color, linewidth=1.4,
                             label=f"{label} (r={radius:.2f})"))
        ax.plot(*center, "x", color=color, markersize=6)

    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.02), fontsize=8, ncol=2)
    set_bounds(ax, [P_A, P_B, P_C, I_A, I_B, I_C],
               radii=[(I_A, r_a), (I_B, r_b), (I_C, r_c)])
    ax.axis("off")

    fig.tight_layout(rect=[0, 0.03, 1, 0.96])
    plt.show()


# --------------------------------------------------------------------------
# Entry point
# --------------------------------------------------------------------------

def main():
    a, b, c = input_triangle_sides()
    results = compute_all(a, b, c)
    print_report(a, b, c, results)

    if HAVE_MPL:
        visualize(a, b, c, results)
    else:
        print("(Install matplotlib with `pip install matplotlib` to see the "
              "visual diagrams: triangle + circles/centers + cevians + exradii.)")


if __name__ == "__main__":
    main()
