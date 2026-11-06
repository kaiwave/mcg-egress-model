# mcg-egress-model

This repository models post-match heavy rail egress and bottleneck dynamics around the Melbourne Cricket Ground (MCG) precinct over a 90-minute departure window. Using a coupled ordinary differential equation (ODE) system across Jolimont ($J$), Richmond ($R$), and Flinders Street ($F$) railway stations, the project investigates how AFL club demographic geographic distributions, walking lags, dynamic gate-throttling feedbacks, and density-dependent boarding friction govern passenger clearance efficiency. It was originally created as my semester Build Your Own Model project for [EVSC20007 Modelling the Real World](https://handbook.unimelb.edu.au/2026/subjects/evsc20007) at the University of Melbourne, but has been tweaked and updated since.

> **Status:** Complete. The entire report has been written up and published on [cloud.kaiwave.dev](https://cloud.kaiwave.dev/s/qAYrEqGSegxdxYR) - with a card on [kaiwave.dev/projects](https://kaiwave.dev/projects/). 

---

## 1. Mathematical Formulation

Let the corridors be indexed by $i \in \{J, R, F\}$:

* **1**: Jolimont ($J$)

* **2**: Richmond ($R$)

* **3**: Flinders Street ($F$)

For each corridor $i$, the model tracks three state variables over time $t \in [0, 90]$ minutes:

* $M_i(t)$: Remaining stadium crowd surge bound for station $i$ (pax)

* $Q_i(t)$: Platform and immediate concourse queue backlog at station $i$ (pax)

* $C_i(t)$: Cumulative passenger count cleared by departing trains through station $i$ by time $t$ (pax)


### Governing ODEs (Mass Conservation)

The mass-balance system governing the crowd evolution is:

$$\frac{dM_i}{dt} = -I_i(t)$$

$$\frac{dQ_i}{dt} = I_i(t) - O_i(t)$$

$$\frac{dC_i}{dt} = O_i(t)$$

where $I_i(t)$ and $O_i(t)$ denote instantaneous passenger inflow and outflow rates (pax/min). Total stadium egress achieved across all corridors at the conclusion of the operational window is:

$$\text{Total Cleared} = \sum_{i \in \{J, R, F\}} C_i(90)$$

---

## 2. Dynamic Feedback Mechanisms

```
Total Attendance (N_total)
         │
         ▼ (Demographic Share α_i)
Initial Station Crowd M_i(0)
         │
         ▼
[ Inflow Rate I_i(t) ] ◄── Arrival Hazard h_i(t) (Gamma)
         │             ◄── Gate Throttle Ψ_i(Q_i) (Balancing Loop)
         ▼
[ Station Queue Q_i(t) ] ──► Platform Density D_i(t)
         │                          │
         │                          ▼
         │                   Boarding Efficiency η(D_i) (log-logistic)
         │                          │
         ▼                          ▼
[ Outflow Rate O_i(t) ] ◄── Nominal Base Capacity C_{base,i}
         │
         ▼
[ Cleared Passengers C_i(t) ]

```

### Inflow & Surge Dynamics

Total rail-bound stadium attendance $N_{\text{total}}$ is partitioned into corridor allocations by demographic weights $\alpha_i$, such that $\sum_i \alpha_i = 1$ and $M_i(0) = \alpha_i \cdot N_{\text{total}}$.

Unimpeded pedestrian arrival follows a unimodal Gamma distribution discharge surge:

$$f_i(t) = \frac{t^{k-1} e^{-t / \theta_i}}{\theta_i^k \, \Gamma(k)}, \quad t \ge 0$$

where $k$ is the shape parameter, $\tau_i$ is the corridor-specific walking travel time, and $\theta_i = \tau_i / k$ is the scale parameter, yielding a peak arrival at $t_{\text{peak}, i} = (k - 1)\theta_i = \tau_i \frac{k - 1}{k}$.

Surge entry into the station concourse is regulated by a non-linear gate-throttling balancing feedback loop representing turnstile saturation and fan departure deferrals:

$$I_i(t) = M_i(t) \cdot h_i(t) \cdot \Psi_i(Q_i), \qquad h_i(t) = \frac{f_i(t)}{1 - F_i(t)}$$

where $F_i$ is the Gamma cumulative distribution and $h_i$ its hazard rate, the arrival rate per fan still outside. With the gate open, this reproduces the unthrottled surge $M_i(0) \cdot f_i(t)$ exactly; when the gate chokes, held-back fans stay in $M_i$ and keep arriving as the queue eases rather than dropping out.

$$\Psi_i(Q_i) = \frac{1}{1 + \exp\left(\beta_q (Q_i - Q_{i,\text{crit}})\right)}$$

* When $Q_i \ll Q_{i,\text{crit}}$, concourse entry is unconstrained ($\Psi_i \to 1.0$).


* When $Q_i = Q_{i,\text{crit}}$, inflow capacity is halved ($\Psi_i = 0.5$).


* When $Q_i \gg Q_{i,\text{crit}}$, the gate chokes inflow toward zero ($\Psi_i \to 0$).


### Platform Congestion & Outflow Friction

Live platform crowding density is computed from the physical station footprint $A_i$ as $D_i(t) = Q_i(t) / A_i$ (pax/m²). Boarding efficiency degrades non-linearly according to a log-logistic curve (a four-parameter logistic, 4PL, dose-response model with its upper asymptote fixed at 1, leaving three free parameters):

$$\eta(D_i) = \eta_{\min} + \frac{1 - \eta_{\min}}{1 + \left( \frac{D_i(t)}{D_{\text{crit}}} \right)^\gamma}$$

Outflow is the degraded train clearing capacity, a rate in pax/min:

$$O_i(t) = C_{\text{base}, i} \cdot \eta(D_i)$$

where nominal clearance capacity is defined by dispatch rate and train capacity: $C_{\text{base}, i} = f_{\text{train}, i} \cdot Cap_{\text{train}}$. The system is integrated with forward Euler at $\Delta t = 0.1$ min (results change by under 0.05 percentage points at $\Delta t = 0.01$). The queue cannot be over-drained within a step: the discrete scheme caps each step's departures at the passengers actually present.

---

## 3. Parameter Specifications

| Parameter / Variable | Symbol | Model Value | Plausible Range | Units | System Role & Calibration |
| --- | --- | --- | --- | --- | --- |
| **Total Rail Attendance** | $N_{\text{total}}$ | $60,000$ | $80,000$ for rivalries (unexercised) | pax | Total crowd departing via train; scaled from match attendance. |
| **Corridor Demographics** | $\alpha_i$ | Per fixture ($\sum \alpha_i = 1$) | – | dimensionless | Splits crowd based on playing team home territory (Section 4). |
| **Mean Walk / Lag Time** | $\tau_i$ | $J: 7.0; \, R: 9.0; \, F: 18.0$ | – | min | Estimated from maps of walking routes through Yarra Park / Birrarung Marr. |
| **Gamma Shape Factor** | $k$ | $4.0$ | $3.0\text{--}5.0$ | dimensionless | Regulates arrival skewness and unimodal arrival peak. |
| **Gamma Scale Factor** | $\theta_i$ | $\tau_i / k$ | – | min | Dispersion parameter for arrival distribution. |
| **Queue Throttle Sensitivity** | $\beta_q$ | $0.002$ | $0.001\text{--}0.003$ | $\text{pax}^{-1}$ | Controls logistic gate-choking sensitivity. |
| **Concourse Gate Capacity** | $Q_{i,\text{crit}}$ | $J: 4,000; \, R: 18,000; \, F: 25,000$ | – | pax | Physical gate throughput and turnstile holding saturation. |
| **Platform Surface Area** | $A_i$ | $J: 1,800; \, R: 9,000; \, F: 15,000$ | – | $\text{m}^2$ | Footprint across platforms (J: 2, R: 10, F: 13 platforms). |
| **Critical Boarding Density** | $D_{\text{crit}}$ | $2.25$ | $2.0\text{--}2.5$ | $\text{pax/m}^2$ | Boarding dwell-delay inflection point (Fruin LOS E/F). |
| **Friction Exponent** | $\gamma$ | $4.0$ | $3.5\text{--}4.5$ | dimensionless | Curvature of boarding efficiency breakdown. |
| **Minimum Efficiency Floor** | $\eta_{\min}$ | $0.35$ | $0.30\text{--}0.40$ | dimensionless | Enforced train departure efficiency threshold. |
| **Train Dispatch Rate** | $f_{\text{train}, i}$ | $J: 0.4; \, R: 0.8; \, F: 1.0$ | $J: 0.3\text{--}0.5; \, R: 0.8\text{--}1.2; \, F: 0.6\text{--}1.0$ | trains/min | Post-match clearing timetable dispatch frequency. |
| **Train Crush Capacity** | $Cap_{\text{train}}$ | $1,000$ | – | pax/train | Nominal capacity for standard 6-car EMU sets. |
| **Nominal Outflow Rate** | $C_{\text{base}, i}$ | $f_{\text{train}, i} \cdot Cap_{\text{train}}$ | – | pax/min | Baseline clearance rate prior to congestion degradation. |

---

## 4. Fanbase Corridor Demographics

Attendance assumes 30,000 fans per team (plus 10,000 fans per side during major rivalry fixtures). For each team, 18% of supporters are uniformly distributed across the three corridors (6% each), while the remaining 82% follow traditional home heartland demographic lines:

| Team | Proportion Jolimont ($\alpha_J$) | Proportion Richmond ($\alpha_R$) | Proportion Flinders ($\alpha_F$) | Catchment / Routing Notes |
| --- | --- | --- | --- | --- |
| **Collingwood** | 0.70 | 0.15 | 0.15 | Heavy Hurstbridge / Mernda line dominance via Jolimont. |
| **Hawthorn** | 0.10 | 0.70 | 0.20 | Eastern suburban rail corridor via Richmond. |
| **Melbourne (Demons)** | 0.20 | 0.30 | 0.50 | Broad distribution with central city alignment. |
| **Richmond** | 0.22 | 0.72 | 0.06 | Strong Eastern / South-Eastern rail dependence. |
| **Essendon** | 0.06 | 0.68 | 0.26 | Broad northern/cross-city split. |
| **Carlton** | 0.16 | 0.16 | 0.68 | Predominantly city-facing and northern transit via Flinders. |
| **St Kilda** | 0.06 | 0.68 | 0.26 | Bayside/Frankston flow channeling through Richmond. |
| **Geelong** | 0.06 | 0.33 | 0.61 | V/Line through Southern Cross/Flinders; Sandringham connection at Richmond. |

---

## 5. Headline Results

Across all 28 pairings of the eight parameterised teams, at 60,000 rail-bound fans:

* 21 fixtures clear their whole crowd within 90 minutes, taking 33 to 61 minutes to clear 95%; Richmond's share sets how long ($r = 0.85$).
* Only the seven Collingwood fixtures fail to clear (84% to 95%), with the shortfall queued at or waiting outside Jolimont, the smallest corridor (400 pax/min).
* Raising Jolimont's dispatch rate from 0.4 to 0.45–0.55 trains/min for Collingwood fixtures alone clears every crowd; widening Jolimont's gates instead makes clearance worse, since the bottleneck is the trains.

The full analysis, figures and limitations are in `main.ipynb` (Discussion section).

---

## 6. Repository Layout

* `main.ipynb`: the executed report (Quarto source, outputs included).
* `model.py`: model parameters, flow functions, the forward Euler simulation, and plotting.
* `assets/`: the stock-and-flow diagram and the figures exported by `main.ipynb`.
* `requirements.txt`: Python dependencies.

### Running it

From the repository root:

```bash
pip install -r requirements.txt
quarto render main.ipynb
```

Rendering uses the stored outputs; add `--execute` to rerun the notebook (which also re-exports the figures to `assets/`). The PDF needs a LaTeX install with `xelatex` (`quarto install tinytex`).

---

## References

1. Fruin, J. (1971). *Pedestrian Planning and Design*. Metropolitan Association of Urban Designers and Environmental Planners.

2. Ritz, C., Baty, F., Streibig, J. C., & Gerhard, D. (2015). *Dose-Response Analysis Using R*. *PLOS ONE*, 10(12), e0146021. [journals.plos.org](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0146021)
