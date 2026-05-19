"""
NUMA01 — Final Project: 2D Integral över triangulerat område Ω.

Definierar klassen Mesh som följer uppgifts-PDF:en
(final_project/project_2Dintegral (1).pdf):

  - Task 1: __init__(mesh) där mesh = (coords, nodes).
  - Task 2: determinant(i) — Jacobi-determinanten för element i.
  - Task 3: min_angle(i)  — minsta vinkeln (radianer) i element i.
  - Task 4: _compute_determinants() — validerar och förberäknar; körs i __init__.
  - Task 5: integrate(f)  — ∫_Ω f dx dy med 3-punkts hörnkvadratur.
  - Task 6: area()        — total area; ska matcha integrate(f≡1).
  - Task 7: plot(ax)      — varje triangel ritas separat.

Matematisk bakgrund (Appendix A i PDF:en):
  - Variabelbyte från godtycklig triangel T till enhetstriangeln T_unit
    ger ∫_T f(x,y) dx dy = ∫_T_unit f̄(ξ,η) |D| dξ dη, där D är
    determinanten av Jacobianen J = [[x2-x1, x3-x1], [y2-y1, y3-y1]].
  - 3-punkts hörnkvadratur på T_unit ger ∫_T f ≈ |D|/6 · (f(P1)+f(P2)+f(P3)).
"""

import numpy as np
from numpy import sqrt
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection


class Mesh:
    """Triangelmesh över ett område Ω ⊂ R²."""

    # Tröskel för "tunn triangel"-validering. 1° är en pragmatisk gräns:
    # tillräckligt lågt för normala FEM-mesher (5-10° vanligt), men fångar
    # nära-kollineära element som degraderar kvadraturen.
    MIN_ANGLE_RAD = np.deg2rad(1.0)

    def __init__(self, mesh):
        """Initialisera från en mesh-tuple.

        mesh : tuple (coords, nodes)
            coords : (2, N)-array. Rad 0 = x-koordinater, rad 1 = y-koordinater.
            nodes  : (3, M)-array av 1-indexerade nodnummer (PDF-konvention).

        Vid initialisering valideras alla element via _compute_determinants;
        ValueError kastas om någon triangel är för smal.
        """
        coords, nodes = mesh
        self.coords = np.asarray(coords, dtype=float)
        # 1-indexerat → 0-indexerat sker här, en enda gång.
        self.elements = np.asarray(nodes, dtype=int) - 1

        # Förberäkna per-element hörnkoordinater för vektoriserade operationer.
        # Varje attribut är en (M,)-array.
        i1, i2, i3 = self.elements
        self.x1, self.y1 = self.coords[0, i1], self.coords[1, i1]
        self.x2, self.y2 = self.coords[0, i2], self.coords[1, i2]
        self.x3, self.y3 = self.coords[0, i3], self.coords[1, i3]

        # Task 4-krav: validering + förberäkning vid initialisering.
        self.determinants = self._compute_determinants()

    @classmethod
    def from_files(cls, coord_path, node_path):
        """Ladda mesh från textfiler enligt PDF:ens format.

        Filformat:
          coord_path : 2 rader (x, y).
          node_path  : 3 rader, 1-indexerade nodnummer.
        """
        coords = np.loadtxt(coord_path)
        # Nodfilerna är skrivna som flyttal ("1.0000000e+00"); läs som float
        # och konvertera till int eftersom np.loadtxt(dtype=int) inte parsar
        # exponentnotation.
        nodes = np.loadtxt(node_path).astype(int)
        return cls((coords, nodes))

    # ---- Task 2 ----
    def determinant(self, i):
        """Determinanten D av Jacobianen för element i.

        D = (x2-x1)(y3-y1) - (x3-x1)(y2-y1).
        |D| är areaförstoringen från enhetstriangeln (area 1/2) till
        triangel i (area |D|/2).
        """
        return (
            (self.x2[i] - self.x1[i]) * (self.y3[i] - self.y1[i])
            - (self.x3[i] - self.x1[i]) * (self.y2[i] - self.y1[i])
        )

    # ---- Task 3 ----
    def min_angle(self, i):
        """Minsta vinkeln (radianer) i triangel i.

        Beräknas via skalärprodukt mellan kantvektorer som UTGÅR FRÅN
        varje hörn — annars får man supplementvinkeln. np.clip skyddar
        mot flyttalsavrundning ut ur [-1, 1] som annars ger NaN i arccos.
        """
        v12 = np.array([self.x2[i] - self.x1[i], self.y2[i] - self.y1[i]])
        v13 = np.array([self.x3[i] - self.x1[i], self.y3[i] - self.y1[i]])
        v23 = np.array([self.x3[i] - self.x2[i], self.y3[i] - self.y2[i]])

        def angle(u, w):
            cos_theta = np.dot(u, w) / (np.linalg.norm(u) * np.linalg.norm(w))
            return np.arccos(np.clip(cos_theta, -1.0, 1.0))

        return min(
            angle(v12, v13),    # vinkeln vid P1: vektorer P1→P2 och P1→P3
            angle(-v12, v23),   # vinkeln vid P2: vektorer P2→P1 och P2→P3
            angle(-v13, -v23),  # vinkeln vid P3: vektorer P3→P1 och P3→P2
        )

    # ---- Task 4 ----
    def _compute_determinants(self):
        """Vektoriserad lista av determinanter över alla element.

        Validerar samtidigt att ingen triangel har en vinkel under
        MIN_ANGLE_RAD; kastar ValueError annars.
        """
        D = (
            (self.x2 - self.x1) * (self.y3 - self.y1)
            - (self.x3 - self.x1) * (self.y2 - self.y1)
        )

        for idx in range(len(D)):
            theta = self.min_angle(idx)
            if theta < self.MIN_ANGLE_RAD:
                raise ValueError(
                    f"Element {idx}: minsta vinkeln "
                    f"{np.rad2deg(theta):.3f}° < tröskel "
                    f"{np.rad2deg(self.MIN_ANGLE_RAD):.3f}°"
                )
        return D

    # ---- Task 5 ----
    def integrate(self, f):
        """∫_Ω f(x,y) dx dy ≈ Σ_i |D_i|/6 · (f(P1_i) + f(P2_i) + f(P3_i)).

        f ska vara vektoriserbar: f(x, y) ska kunna ta NumPy-arrayer
        av samma form och returnera en array.
        """
        f1 = f(self.x1, self.y1)
        f2 = f(self.x2, self.y2)
        f3 = f(self.x3, self.y3)
        return np.sum(np.abs(self.determinants) * (f1 + f2 + f3) / 6.0)

    # ---- Task 6 ----
    def area(self):
        """Total area = Σ |D_i|/2. Ska matcha integrate(f≡1)."""
        return 0.5 * np.sum(np.abs(self.determinants))

    # ---- Task 7 ----
    def plot(self, ax=None, edgecolor="black", facecolor="none", linewidth=0.5):
        """Rita mesh:en. Varje triangel syns separat via PolyCollection.

        PolyCollection är snabbare än per-triangel plt.plot och låter
        oss färgsätta per element om vi vill.
        """
        if ax is None:
            _, ax = plt.subplots()
        triangles = np.stack(
            [
                np.column_stack([self.x1, self.y1]),
                np.column_stack([self.x2, self.y2]),
                np.column_stack([self.x3, self.y3]),
            ],
            axis=1,
        )
        collection = PolyCollection(
            triangles,
            edgecolors=edgecolor,
            facecolors=facecolor,
            linewidths=linewidth,
        )
        ax.add_collection(collection)
        ax.set_aspect("equal")
        ax.autoscale_view()
        return ax
