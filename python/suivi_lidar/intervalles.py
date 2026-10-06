"""Arithmétique d'intervalles, en remplacement d'INTLAB.

Un objet Intervalle contient deux tableaux numpy de même forme, `inf` et `sup`.
Un vecteur d'intervalles de taille 9 est donc une boîte de l'espace d'état,
et un tableau 9 x N contient N boîtes (une par colonne), comme dans MATLAB.

Comme INTLAB, chaque opération arrondit vers l'extérieur : la borne inférieure
est descendue et la borne supérieure montée d'au moins un ulp, pour que le
résultat contienne toujours la vraie valeur malgré les arrondis de la machine.

Correspondance avec INTLAB :
    midrad(m, r)   -> Intervalle.midrad(m, r)
    infsup(a, b)   -> Intervalle.infsup(a, b)
    intval(A) * x  -> produit(A, x)       (A matrice réelle, x intervalle)
    mid(x), rad(x) -> x.mid(), x.rad()
    inf(x), sup(x) -> x.inf, x.sup
    intersect(a,b) -> intersect(a, b)     (NaN si l'intersection est vide, comme INTLAB)
"""

import numpy as np

_EPS = np.finfo(float).eps


def _bas(x):
    return np.nextafter(x, -np.inf)


def _haut(x):
    return np.nextafter(x, np.inf)


class Intervalle:
    __slots__ = ("inf", "sup")
    __array_priority__ = 1000  # pour que numpy laisse la main à nos opérateurs

    def __init__(self, inf, sup=None):
        inf = np.array(inf, dtype=float)
        sup = inf.copy() if sup is None else np.array(sup, dtype=float)
        if inf.shape != sup.shape:
            raise ValueError("inf et sup doivent avoir la même forme")
        self.inf = inf
        self.sup = sup

    # ------------------------------------------------------------------ constructeurs
    @staticmethod
    def midrad(m, r):
        m = np.asarray(m, dtype=float)
        r = np.broadcast_to(np.asarray(r, dtype=float), m.shape)
        return Intervalle(_bas(m - r), _haut(m + r))

    @staticmethod
    def infsup(a, b):
        return Intervalle(a, b)

    @staticmethod
    def zeros(forme):
        return Intervalle(np.zeros(forme), np.zeros(forme))

    def copy(self):
        return Intervalle(self.inf.copy(), self.sup.copy())

    # ------------------------------------------------------------------ propriétés
    @property
    def shape(self):
        return self.inf.shape

    def __len__(self):
        return len(self.inf)

    def mid(self):
        return 0.5 * (self.inf + self.sup)

    def rad(self):
        return _haut(0.5 * (self.sup - self.inf))

    def largeur(self):
        return self.sup - self.inf

    def contient_nan(self):
        return bool(np.any(np.isnan(self.inf)) or np.any(np.isnan(self.sup)))

    # ------------------------------------------------------------------ indexation (comme X(1:3) ou X(:, b))
    def __getitem__(self, cle):
        return Intervalle(self.inf[cle], self.sup[cle])

    def __setitem__(self, cle, valeur):
        if isinstance(valeur, Intervalle):
            self.inf[cle] = valeur.inf
            self.sup[cle] = valeur.sup
        else:
            v = np.asarray(valeur, dtype=float)
            self.inf[cle] = v
            self.sup[cle] = v

    # ------------------------------------------------------------------ opérations
    @staticmethod
    def _bornes(x):
        if isinstance(x, Intervalle):
            return x.inf, x.sup
        v = np.asarray(x, dtype=float)
        return v, v

    def __add__(self, autre):
        a, b = self._bornes(autre)
        return Intervalle(_bas(self.inf + a), _haut(self.sup + b))

    __radd__ = __add__

    def __neg__(self):
        return Intervalle(-self.sup, -self.inf)

    def __sub__(self, autre):
        a, b = self._bornes(autre)
        return Intervalle(_bas(self.inf - b), _haut(self.sup - a))

    def __rsub__(self, autre):
        return (-self) + autre

    def __mul__(self, autre):
        if isinstance(autre, Intervalle):
            produits = np.stack([self.inf * autre.inf, self.inf * autre.sup,
                                 self.sup * autre.inf, self.sup * autre.sup])
            return Intervalle(_bas(produits.min(axis=0)), _haut(produits.max(axis=0)))
        k = np.asarray(autre, dtype=float)
        p1, p2 = self.inf * k, self.sup * k
        return Intervalle(_bas(np.minimum(p1, p2)), _haut(np.maximum(p1, p2)))

    __rmul__ = __mul__

    def __truediv__(self, autre):
        if isinstance(autre, Intervalle):
            if np.any((autre.inf <= 0) & (autre.sup >= 0)):
                raise ZeroDivisionError("division par un intervalle qui contient zéro")
            return self * Intervalle(_bas(1.0 / autre.sup), _haut(1.0 / autre.inf))
        k = np.asarray(autre, dtype=float)
        q1, q2 = self.inf / k, self.sup / k
        return Intervalle(_bas(np.minimum(q1, q2)), _haut(np.maximum(q1, q2)))

    def carre(self):
        """x^2 exact : [0, max] quand l'intervalle contient zéro."""
        a2, b2 = self.inf ** 2, self.sup ** 2
        haut = np.maximum(a2, b2)
        bas = np.where((self.inf <= 0) & (self.sup >= 0), 0.0, np.minimum(a2, b2))
        return Intervalle(np.maximum(_bas(bas), 0.0), _haut(haut))

    def racine(self):
        return Intervalle(np.maximum(_bas(np.sqrt(np.maximum(self.inf, 0.0))), 0.0),
                          _haut(np.sqrt(np.maximum(self.sup, 0.0))))

    def somme(self, axis=0):
        n = self.inf.shape[axis]
        bas = self.inf.sum(axis=axis)
        haut = self.sup.sum(axis=axis)
        erreur = n * _EPS * np.maximum(np.abs(self.inf), np.abs(self.sup)).sum(axis=axis)
        return Intervalle(_bas(bas - erreur), _haut(haut + erreur))

    def __repr__(self):
        return f"Intervalle(inf={self.inf!r}, sup={self.sup!r})"


def produit(A, x):
    """Produit d'une matrice réelle par un vecteur (ou une matrice) d'intervalles.

    Équivalent de intval(A) * x dans INTLAB. Le résultat est l'enveloppe exacte,
    élargie de l'erreur d'arrondi possible du produit scalaire.
    """
    A = np.asarray(A, dtype=float)
    Ap, Am = np.maximum(A, 0.0), np.minimum(A, 0.0)
    bas = Ap @ x.inf + Am @ x.sup
    haut = Ap @ x.sup + Am @ x.inf
    n = A.shape[1]
    erreur = n * _EPS * (np.abs(A) @ np.maximum(np.abs(x.inf), np.abs(x.sup)))
    return Intervalle(_bas(bas - erreur), _haut(haut + erreur))


def intersect(a, b):
    """Intersection composante par composante. Composante vide -> NaN (comme INTLAB)."""
    a_inf, a_sup = Intervalle._bornes(a)
    b_inf, b_sup = Intervalle._bornes(b)
    bas = np.maximum(a_inf, b_inf)
    haut = np.minimum(a_sup, b_sup)
    vide = bas > haut
    bas = np.where(vide, np.nan, bas)
    haut = np.where(vide, np.nan, haut)
    return Intervalle(bas, haut)


def concatener_colonnes(liste):
    """Assemble plusieurs boîtes (vecteurs d'intervalles) en un tableau n x k."""
    return Intervalle(np.column_stack([b.inf for b in liste]), np.column_stack([b.sup for b in liste]))
