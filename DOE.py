"""
Tarea 3 - DOE Opcion A (factorial 3 temperaturas x 3 tamanos, 10 replicas)
Usa el Recocido Simulado con 2-opt de la Tarea 2 (Ej. 4), esquema exponencial.
Colocar junto a common.py (mismo directorio que la Tarea 2).
"""
import math
import random
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import statsmodels.api as sm
from statsmodels.formula.api import ols
from scipy import stats
from common import *          # matriz_distancias, instancia_aleatoria, costo_tour

# ----------------------------------------------------------------------
# Configuracion
# ----------------------------------------------------------------------
MATRICULA = 2043917           # <-- CAMBIA por tu matricula
TAMANOS = [5, 10, 15]
REPLICAS = 10
T_MIN, T_MAX = 0.01, 1.0      # rango de T0 (escala de tu instancia en [0,1]^2)
TF = 1e-3                     # temperatura final (igual que en Tarea 2)
N_ITER = 500                  # presupuesto fijo; calibrar ANTES del experimento final
N_BOOT = 10_000
ALFA = 0.05
ORDEN = ["bajo", "medio", "alto"]

random.seed(MATRICULA)
r = random.random()           # se genera UNA sola vez
FRAC = {"bajo": 0.3 * r, "medio": 0.3 * r + 0.3, "alto": 0.3 * r + 0.6}
NIVELES = {k: T_MIN + f * (T_MAX - T_MIN) for k, f in FRAC.items()}
print(f"r = {r:.4f}")
for k in ORDEN:
    print(f"  {k}: fraccion={FRAC[k]:.3f}  T0={NIVELES[k]:.4f}")


# ----------------------------------------------------------------------
# Tu recocido de la Tarea 2, con T0 como parametro (esquema exponencial)
# ----------------------------------------------------------------------
def temperatura_exp(k, N, T0):
    return T0 * (TF / T0) ** (k / N)


def recocido(D, T0, N, rnd):
    n = len(D)
    tour = list(range(n)); rnd.shuffle(tour)
    costo = costo_tour(D, tour); mejor = costo
    for k in range(N):
        T = temperatura_exp(k, N, T0)
        i, j = sorted(rnd.sample(range(n), 2))
        if i == 0 and j == n - 1:
            continue
        a, b, c, d = tour[i - 1], tour[i], tour[j], tour[(j + 1) % n]
        delta = D[a, c] + D[b, d] - D[a, b] - D[c, d]
        if delta < 0 or rnd.random() < math.exp(-delta / T):
            tour[i:j + 1] = reversed(tour[i:j + 1]); costo += delta
            if costo < mejor - 1e-12:
                mejor = costo
    return mejor


# ----------------------------------------------------------------------
# Experimento: 3 x 3 x 10 = 90 corridas
# ----------------------------------------------------------------------
instancias = {n: matriz_distancias(instancia_aleatoria(n, np.random.default_rng(100 + n)))
              for n in TAMANOS}
filas = []
for n in TAMANOS:
    D = instancias[n]
    for idx, nivel in enumerate(ORDEN):
        for rep in range(REPLICAS):
            semilla = 100_000 * n + 1_000 * idx + rep      # semillas distintas
            c = recocido(D, NIVELES[nivel], N_ITER, random.Random(semilla))
            filas.append({"tamano": str(n), "temperatura": nivel,
                          "T0": NIVELES[nivel], "rep": rep, "costo": round(c, 9)})
df = pd.DataFrame(filas)
df.to_csv("results/t3_resultados_doe.csv", index=False)

# ----------------------------------------------------------------------
# ANOVA de dos factores
# ----------------------------------------------------------------------
modelo = ols("costo ~ C(temperatura) * C(tamano)", data=df).fit()
anova = sm.stats.anova_lm(modelo, typ=2)
print("\n=== ANOVA de dos factores ===\n", anova)
anova.to_csv("results/t3_anova.csv")
p_temp = anova.loc["C(temperatura)", "PR(>F)"]
p_int = anova.loc["C(temperatura):C(tamano)", "PR(>F)"]
print(f"Temperatura significativa: {p_temp < ALFA} (p={p_temp:.4g})")
print(f"Interaccion significativa: {p_int < ALFA} (p={p_int:.4g})")

# Apoyo: efecto de la temperatura dentro de cada tamano + Levene (homocedasticidad)
print("\n=== ANOVA de un factor por tamano ===")
p_tam = {}
for n in TAMANOS:
    g = [df[(df.tamano == str(n)) & (df.temperatura == k)].costo.values for k in ORDEN]
    if all(np.ptp(x) == 0 for x in g) and len({x[0] for x in g}) == 1:
        p_tam[n] = np.nan
        print(f"n={n}: sin variabilidad (todas las corridas iguales)")
        continue
    F, p = stats.f_oneway(*g)
    p_tam[n] = p
    print(f"n={n}: F={F:.3f}  p={p:.4g}")
lev = stats.levene(*[df[df.tamano == str(n)].costo.values for n in TAMANOS])
print(f"Levene entre tamanos: p={lev.pvalue:.4g} (p<0.05 => varianzas distintas)")


# ----------------------------------------------------------------------
# Bootstrap de la mediana
# ----------------------------------------------------------------------
def boot_med(x, seed):
    rng = np.random.default_rng(seed)
    x = np.asarray(x)
    m = np.median(rng.choice(x, size=(N_BOOT, len(x)), replace=True), axis=1)
    return np.median(m), *np.percentile(m, [100 * ALFA / 2, 100 * (1 - ALFA / 2)])


fb = []
for n in TAMANOS:
    for k in ORDEN:
        x = df[(df.tamano == str(n)) & (df.temperatura == k)].costo.values
        med, lo, hi = boot_med(x, n * 10 + ORDEN.index(k))
        fb.append({"tamano": n, "nivel": k, "T0": NIVELES[k], "media": x.mean(),
                   "mediana_boot": med, "IC95_inf": lo, "IC95_sup": hi})
boot = pd.DataFrame(fb)
boot.to_csv("results/t3_bootstrap.csv", index=False)
print("\n=== Bootstrap por punto de diseno ===\n", boot.round(4).to_string(index=False))

# Mejor temperatura por tamano (empates por mediana se rompen con la media)
mejores = (boot.sort_values(["tamano", "mediana_boot", "media"])
           .groupby("tamano").head(1).copy())
mejores["p_temp_en_n"] = mejores["tamano"].map(p_tam)
mejores["significativo"] = mejores["p_temp_en_n"] < ALFA
print("\n=== Temperatura optima por tamano ===\n", mejores.round(4).to_string(index=False))
mejores.to_csv("results/t3_optimos.csv", index=False)

# ----------------------------------------------------------------------
# Graficas
# ----------------------------------------------------------------------
fig, ax = plt.subplots(1, 2, figsize=(10, 4))
ax[0].plot(ORDEN, df.groupby("temperatura").costo.mean().reindex(ORDEN).values, "o-")
ax[0].set(title="Efecto principal: temperatura", xlabel="Nivel", ylabel="Costo medio")
ax[1].plot([str(n) for n in TAMANOS],
           df.groupby("tamano").costo.mean().reindex([str(n) for n in TAMANOS]).values,
           "o-", color="tab:orange")
ax[1].set(title="Efecto principal: tamano", xlabel="n", ylabel="Costo medio")
for a in ax: a.grid(alpha=.3)
plt.tight_layout(); plt.savefig("figures/t3_efectos_principales.png", dpi=150); plt.close()

med = df.groupby(["tamano", "temperatura"]).costo.mean().unstack().reindex(columns=ORDEN)
fig, ax = plt.subplots(1, 2, figsize=(11, 4))
for n in TAMANOS:
    y = med.loc[str(n)]
    ax[0].plot(ORDEN, y.values, "o-", label=f"n={n}")
    ax[1].plot(ORDEN, (y / y.mean()).values, "o-", label=f"n={n}")
ax[0].set(title="Interaccion T x n (costo medio)", xlabel="Temperatura", ylabel="Costo medio")
ax[1].set(title="Interaccion T x n (normalizado por n)", xlabel="Temperatura",
          ylabel="Costo / media del tamano")
for a in ax: a.grid(alpha=.3); a.legend()
plt.tight_layout(); plt.savefig("figures/t3_interaccion.png", dpi=150); plt.close()

fig, ax = plt.subplots(1, 3, figsize=(13, 4))
for a, n in zip(ax, TAMANOS):
    b = boot[boot.tamano == n]
    a.errorbar(b.nivel, b.mediana_boot,
               yerr=[b.mediana_boot - b.IC95_inf, b.IC95_sup - b.mediana_boot],
               fmt="o-", capsize=4)
    a.set(title=f"n = {n}", xlabel="Temperatura", ylabel="Mediana bootstrap (IC 95%)")
    a.grid(alpha=.3)
plt.tight_layout(); plt.savefig("figures/t3_bootstrap_ic.png", dpi=150); plt.close()
print("\nListo.")