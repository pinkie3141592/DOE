# Tarea 3 – Ajuste de temperatura del Recocido Simulado mediante DOE

Ajuste del parámetro de temperatura inicial **T₀** del Recocido Simulado (movimiento 2-opt, enfriamiento exponencial) sobre el TSP, usando un diseño factorial completo (**Opción A**), ANOVA de dos factores y bootstrap.

- **Autora:** Rebeca Sanchez Rodriguez (matrícula 2043917)
- **Materia:** Algoritmia – FCFM
- **Reporte:** `reporte_tarea3.pdf`

## Diseño experimental

| Elemento | Valor |
|---|---|
| Algoritmo | Recocido Simulado + 2-opt (Tarea 2, Ej. 4), enfriamiento exponencial de T₀ a T_F = 10⁻³ |
| Iteraciones por corrida | 500 |
| Factor | Temperatura inicial T₀ (3 niveles) |
| Bloque | Tamaño del TSP, n ∈ {5, 10, 15} (una instancia por tamaño, puntos aleatorios en [0,1]²) |
| Réplicas | 10 por punto de diseño, con semillas distintas (90 corridas en total) |
| Niveles | r = 0.6553 (semilla = matrícula), rango [0.01, 1.0] → T₀ = 0.2046, 0.5016, 0.7986 |
| Análisis | ANOVA de dos factores (α = 0.05), ANOVA por tamaño, prueba de Levene, bootstrap de la mediana (10 000 remuestreos, IC 95 %) |

## Resultados principales

- Solo el **tamaño** del problema es significativo (p ≈ 6×10⁻⁷⁵). La temperatura (p = 0.726) y la interacción (p = 0.862) **no** lo son.
- En n = 5 y n = 10 todas las corridas alcanzan el mismo costo, por lo que cualquier temperatura es equivalente.
- En n = 15 el nivel bajo (T₀ = 0.2046) tiene la menor mediana, pero la diferencia no es significativa (p = 0.727) y los IC se traslapan.
- No hay una tendencia sustentada entre n y la temperatura óptima.

| n | Mejor T₀ | Mediana bootstrap | IC 95 % | ¿Significativo? |
|---|---|---|---|---|
| 5 | Cualquiera (empate) | 2.2668 | [2.2668, 2.2668] | No |
| 10 | Cualquiera (empate) | 2.2695 | [2.2695, 2.2695] | No |
| 15 | 0.2046 (nominal) | 4.0752 | [3.9037, 4.2141] | No (p = 0.727) |

## Estructura del repositorio

```
.
├── DOE.py                 # experimento completo: SA, ANOVA, bootstrap y gráficas
├── graficas_t3.py         # regenera las figuras a partir de los CSV
├── common.py              # utilidades de la Tarea 2 (instancias, distancias, costo)
├── reporte_tarea3.pdf     # reporte final
├── results/               # CSV generados
│   ├── t3_resultados_doe.csv   # las 90 corridas
│   ├── t3_anova.csv            # tabla ANOVA de dos factores
│   ├── t3_bootstrap.csv        # mediana bootstrap e IC por punto de diseño
│   └── t3_optimos.csv          # mejor temperatura por tamaño
└── figures/               # gráficas del reporte
    ├── t3_efectos_principales.png
    ├── t3_interaccion.png
    └── t3_bootstrap_ic.png
```

## Requisitos

- Python 3.10 o superior (probado con 3.14.9)
- Librerías: `numpy`, `pandas`, `scipy`, `statsmodels`, `matplotlib`

```bash
pip install numpy pandas scipy statsmodels matplotlib
```

## Cómo reproducir

```bash
# 1. Crear las carpetas de salida si no existen
mkdir -p results figures

# 2. Correr el experimento (90 corridas, ANOVA, bootstrap y gráficas)
python DOE.py

# 3. (Opcional) Regenerar solo las gráficas desde los CSV
python graficas_t3.py
```

El experimento es determinista: `r` se genera una sola vez a partir de la matrícula y cada corrida usa una semilla fija, así que los resultados se reproducen exactamente.

## Notas

- `common.py` debe estar en la misma carpeta que `DOE.py`.
- Para usar otra matrícula, cambia `MATRICULA` al inicio de `DOE.py`; los niveles de temperatura se recalculan automáticamente.
- La homocedasticidad del ANOVA global se viola (Levene, p < 0.001) porque n = 5 y n = 10 no tienen variabilidad; los detalles se discuten en el reporte.
