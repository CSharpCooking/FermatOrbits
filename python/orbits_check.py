#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Проверка результатов статьи
Р. Р. Шарипов, Р. Ф. Гибадуллин, «Орбиты и стабилизаторы группы симметрий решений
сравнения Ферма по простому модулю».

Два режима.

1) orbits — для всех нечётных n <= NMAX и простых p <= PMAX с p > n и непустым Omega_{n,p}
   строится множество проективных точек Omega-bar (в координатах (x:y:w), w = -z, x = 1),
   вычисляются орбиты и стабилизаторы группы K_{n,p} = (mu_n x mu_n) ⋊ S_3, и результат
   сравнивается с теоремой 1 (длины орбит, порядки стабилизаторов, число орбит
   k = N/(6d^2) + eps/2 + 2 delta/3) и следствием 2 (критерий транзитивности):

       python3 orbits_check.py orbits 15 1500

2) gl3 — полный поиск всех матриц M из GL_3(F_p), сохраняющих Omega_{3,p} (метод
   образов базиса: выбираются три линейно независимых вектора из Omega, перебираются
   все |Omega|^3 троек их образов):

       python3 orbits_check.py gl3 19

Зависимости: только стандартная библиотека Python 3 (режим gl3 быстрее с NumPy, если он
установлен; без него работает тот же алгоритм на чистом Python).
"""
import sys, math, time, itertools
from collections import defaultdict, Counter

# ---------------------------------------------------------------- общие функции

def primes(lo, hi):
    s = bytearray([1]) * (hi + 1)
    s[0:2] = b'\0\0'
    for i in range(2, int(hi ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = bytearray(len(s[i * i::i]))
    return [q for q in range(lo, hi + 1) if s[q]]

# P_sigma (x1,x2,x3) = (x_{sigma(1)}, x_{sigma(2)}, x_{sigma(3)})
PERMS = {'id': (0, 1, 2), '(12)': (1, 0, 2), '(13)': (2, 1, 0), '(23)': (0, 2, 1),
         '(123)': (1, 2, 0), '(132)': (2, 0, 1)}

# ---------------------------------------------------------------- режим orbits

def analyse(n, p):
    """Орбиты K_{n,p} на Omega-bar_{n,p} и сравнение с теоремой 1."""
    d = math.gcd(n, p - 1)
    pw = [pow(t, n, p) for t in range(p)]
    roots = defaultdict(list)               # v -> все t из F_p^* с t^n = v
    for t in range(1, p):
        roots[pw[t]].append(t)
    mu = roots[1]                           # mu_n, |mu_n| = d
    pts = {(1, y, w) for y in range(1, p) for w in roots.get((-1 - pw[y]) % p, ())}
    N = len(pts)
    eps = 1 if (p - 2) % p in roots else 0                      # -2 — вычет n-й степени?
    om = [t for t in range(2, p) if pow(t, 3, p) == 1]          # примитивные кубические корни
    delta = 1 if (p % 3 == 1 and om[0] in roots) else 0
    inv = [0] + [pow(t, p - 2, p) for t in range(1, p)]
    K = [(s, a, b) for s in PERMS for a in mu for b in mu]      # элемент (sigma, (a, b, 1))

    def act(g, v):                                              # diag(a,b,1) P_sigma, нормировка x = 1
        s = PERMS[g[0]]
        u = (v[s[0]], v[s[1]], v[s[2]])
        i = inv[u[0] * g[1] % p]
        return (1, u[1] * g[2] * i % p, u[2] * i % p)

    seen, orbits = set(), []
    for v in pts:
        if v in seen:
            continue
        orb = {act(g, v) for g in K}
        seen |= orb
        orbits.append((len(orb), 6 * d * d // len(orb)))        # (длина орбиты, порядок стабилизатора)
    orbits.sort()
    free = N - 3 * d * d * eps - 2 * d * d * delta
    pred = None
    if free >= 0 and free % (6 * d * d) == 0:
        pred = sorted([(3 * d * d, 2)] * eps + [(2 * d * d, 3)] * delta + [(6 * d * d, 1)] * (free // (6 * d * d)))
    k_formula = N / (6 * d * d) + eps / 2 + 2 * delta / 3
    transitive = len(orbits) == 1
    criterion = (N == 6 * d * d) or (eps == 1 and N == 3 * d * d) or (delta == 1 and N == 2 * d * d)
    return dict(n=n, p=p, d=d, N=N, eps=eps, delta=delta, k=len(orbits), orbits=orbits,
                theorem_ok=(pred == orbits and abs(k_formula - len(orbits)) < 1e-9),
                criterion_ok=(transitive == criterion), transitive=transitive,
                hypothesis=(p > n and N > n * n))

def run_orbits(nmax, pmax):
    t0 = time.time()
    R = [analyse(n, p) for n in range(3, nmax + 1, 2) for p in primes(5, pmax) if p > n]
    R = [r for r in R if r['N'] > 0]
    hyp = [r for r in R if r['hypothesis']]
    bad = [r for r in R if not (r['theorem_ok'] and r['criterion_ok'])]
    print('пар с p > n и непустым Omega: %d; из них с N > n^2: %d' % (len(R), len(hyp)))
    print('  по n:', dict(Counter(r['n'] for r in hyp)))
    print('расхождений с теоремой 1 / следствием 2: %d' % len(bad))
    for r in bad:
        print('  ', r)
    tr = [(r['n'], r['p']) for r in hyp if r['transitive']]
    print('транзитивное действие (в условиях теоремы), %d пар:' % len(tr))
    for n in sorted(set(a for a, _ in tr)):
        print('  n = %2d: p =' % n, ', '.join(str(b) for a, b in tr if a == n))
    print('время: %.1f с' % (time.time() - t0))

# ---------------------------------------------------------------- режим gl3

def omega(n, p):
    pw = [pow(t, n, p) for t in range(p)]
    roots = defaultdict(list)
    for z in range(1, p):
        roots[pw[z]].append(z)
    return [(x, y, z) for x in range(1, p) for y in range(1, p) for z in roots.get((pw[x] + pw[y]) % p, ())]

def det3(m, p):
    return (m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
            - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])) % p

def inv3(m, p):
    d = det3(m, p)
    di = pow(d, p - 2, p)
    c = [[(m[(j + 1) % 3][(i + 1) % 3] * m[(j + 2) % 3][(i + 2) % 3]
           - m[(j + 1) % 3][(i + 2) % 3] * m[(j + 2) % 3][(i + 1) % 3]) * di % p
          for j in range(3)] for i in range(3)]
    return c

def run_gl3(p, n=3):
    t0 = time.time()
    O = omega(n, p)
    S = set(O)
    m = len(O)
    print('p = %d: |Omega_{%d,%d}| = %d, троек образов: %d' % (p, n, p, m, m ** 3))
    # три линейно независимых вектора из Omega
    basis = next(([O[i], O[j], O[k]] for i in range(m) for j in range(i + 1, m) for k in range(j + 1, m)
                  if det3([[O[i][r], O[j][r], O[k][r]] for r in range(3)], p)), None)
    B = [[basis[c][r] for c in range(3)] for r in range(3)]      # столбцы — базисные векторы
    Binv = inv3(B, p)
    # координаты каждого v из Omega в базисе: v = B c_v
    coords = [tuple(sum(Binv[r][c] * v[c] for c in range(3)) % p for r in range(3)) for v in O]
    test = coords[: min(8, m)]                                   # быстрый фильтр по нескольким векторам
    found = []
    try:
        import numpy as np
        W = np.array(O, dtype=np.int64)
        C = np.array(coords, dtype=np.int64)
        T = np.zeros((p, p, p), dtype=bool)
        for v in O:
            T[v] = True
        Ct = C[: len(test)]
        for i1 in range(m):
            for i2 in range(m):
                part = (np.outer(Ct[:, 0], W[i1]) + np.outer(Ct[:, 1], W[i2])) % p         # (q,3)
                img = (part[:, None, :] + Ct[:, 2][:, None, None] * W[None, :, :]) % p     # (q,m,3)
                ok = T[img[..., 0], img[..., 1], img[..., 2]].all(axis=0)
                for i3 in np.nonzero(ok)[0]:
                    Wm = [[O[i1][r], O[i2][r], O[i3][r]] for r in range(3)]
                    M = [[sum(Wm[r][k] * Binv[k][c] for k in range(3)) % p for c in range(3)] for r in range(3)]
                    if det3(M, p) == 0:
                        continue
                    im = (np.array(M) @ W.T) % p
                    if T[im[0], im[1], im[2]].all():
                        found.append(M)
    except ImportError:
        for i1 in range(m):
            for i2 in range(m):
                for i3 in range(m):
                    w = (O[i1], O[i2], O[i3])
                    if any(tuple(sum(c[k] * w[k][r] for k in range(3)) % p for r in range(3)) not in S for c in test):
                        continue
                    Wm = [[w[0][r], w[1][r], w[2][r]] for r in range(3)]
                    M = [[sum(Wm[r][k] * Binv[k][c] for k in range(3)) % p for c in range(3)] for r in range(3)]
                    if det3(M, p) == 0:
                        continue
                    if all(tuple(sum(M[r][c] * v[c] for c in range(3)) % p for r in range(3)) in S for v in O):
                        found.append(M)
    mono = sum(1 for M in found if all(sum(1 for x in row if x) == 1 for row in M)
               and all(sum(1 for r in range(3) if M[r][c]) == 1 for c in range(3)))
    d = math.gcd(n, p - 1)
    print('автоморфизмов: %d, из них мономиальных: %d, немономиальных: %d; 6d^2(p-1) = %d'
          % (len(found), mono, len(found) - mono, 6 * d * d * (p - 1)))
    print('время: %.1f с' % (time.time() - t0))

# ---------------------------------------------------------------- запуск

if __name__ == '__main__':
    if len(sys.argv) >= 2 and sys.argv[1] == 'orbits':
        run_orbits(int(sys.argv[2]) if len(sys.argv) > 2 else 15, int(sys.argv[3]) if len(sys.argv) > 3 else 1500)
    elif len(sys.argv) >= 3 and sys.argv[1] == 'gl3':
        run_gl3(int(sys.argv[2]))
    else:
        print(__doc__)
