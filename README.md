# mcg-egress-model

This repository models post-match heavy rail egress and bottleneck dynamics around the Melbourne Cricket Ground (MCG) precinct over a 90-minute departure window. Using a coupled ordinary differential equation (ODE) system across Jolimont ($J$), Richmond ($R$), and Flinders Street ($F$) railway stations, the project investigates how AFL club demographic geographic distributions, walking lags, dynamic gate-throttling feedbacks, and density-dependent boarding friction govern passenger clearance efficiency. It was created as my semester Build Your own Model project for EVSC20007.

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
[ Inflow Rate I_i(t) ] ◄── Incoming Arrival Surge f_i(t) (Gamma PDF)
         │             ◄── Gate Throttle Ψ_i(Q_i) (Balancing Loop)
         ▼
[ Station Queue Q_i(t) ] ──► Platform Density D_i(t)
         │                          │
         │                          ▼
         │                   Boarding Efficiency η(D_i) (4PL)
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

$$I_i(t) = M_i(t) \cdot f_i(t) \cdot \Psi_i(Q_i)$$

$$\Psi_i(Q_i) = \frac{1}{1 + \exp\left(\beta_q (Q_i - Q_{i, \text{gate\_crit}})\right)}$$

* When $Q_i \ll Q_{i, \text{gate\_crit}}$, concourse entry is unconstrained ($\Psi_i \to 1.0$).


* When $Q_i = Q_{i, \text{gate\_crit}}$, inflow capacity is halved ($\Psi_i = 0.5$).


* When $Q_i \gg Q_{i, \text{gate\_crit}}$, the gate chokes inflow toward zero ($\Psi_i \to 0$).


### Platform Congestion & Outflow Friction

Live platform crowding density is computed from the physical station footprint $A_i$ as $D_i(t) = Q_i(t) / A_i$ (pax/m²). Boarding efficiency degrades non-linearly according to a Four-Parameter Logistic (4PL) regression model:

$$\eta(D_i) = \eta_{\min} + \frac{1 - \eta_{\min}}{1 + \left( \frac{D_i(t)}{D_{\text{crit}}} \right)^\gamma}$$

Outflow is constrained by current queue availability and degraded train clearing capacity:

$$O_i(t) = \min\left(Q_i(t), \, C_{\text{base}, i} \cdot \eta(D_i)\right)$$

where nominal clearance capacity is defined by dispatch rate and train capacity: $C_{\text{base}, i} = f_{\text{train}, i} \cdot Cap_{\text{train}}$.

---

## 3. Parameter Specifications

| Parameter / Variable | Symbol | Value / Baseline Range | Units | System Role & Calibration |
| --- | --- | --- | --- | --- |
| **Total Rail Attendance** | $N_{\text{total}}$ | $60,000$ | pax | Total crowd departing via train; scaled from match attendance.

 |
| **Corridor Demographics** | $\alpha_i$ | Dynamic ($\sum \alpha_i = 1$) | dimensionless | Splits crowd based on playing team home territory.

 |
| **Mean Walk / Lag Time** | $\tau_i$ | $J: 7.0, \, R: 9.0, \, F: 18.0$ | min | Calibrated to pedestrian walking time through Yarra Park / Birrarung Marr.

 |
| **Gamma Shape Factor** | $k$ | $3.0\text{--}5.0$ | dimensionless | Regulates arrival skewness and unimodal arrival peak.

 |
| **Gamma Scale Factor** | $\theta_i$ | $\tau_i / k$ | min | Dispersion parameter for arrival distribution.

 |
| **Queue Throttle Sensitivity** | $\beta_q$ | $0.001\text{--}0.003$ | $\text{pax}^{-1}$ | Controls logistic gate-choking sensitivity.

 |
| **Concourse Gate Capacity** | $Q_{i,\text{gate\_crit}}$ | $J: 4,000; \, R: 18,000; \, F: 25,000$ | pax | Physical gate throughput and turnstile holding saturation.

 |
| **Platform Surface Area** | $A_i$ | $J: 1,800; \, R: 9,000; \, F: 15,000$ | $\text{m}^2$ | Footprint across platforms (J: 2, R: 10, F: 15 platforms).

 |
| **Critical Boarding Density** | $D_{\text{crit}}$ | $2.0\text{--}2.5$ | $\text{pax/m}^2$ | Boarding dwell-delay inflection point (Fruin LOS E/F).

 |
| **Friction Exponent** | $\gamma$ | $3.5\text{--}4.5$ | dimensionless | Curvature of boarding efficiency breakdown.

 |
| **Minimum Efficiency Floor** | $\eta_{\min}$ | $0.30\text{--}0.40$ | dimensionless | Enforced train departure efficiency threshold.

 |
| **Train Dispatch Rate** | $f_{\text{train}, i}$ | $J: 0.3\text{--}0.5; \, R: 0.8\text{--}1.2; \, F: 0.6\text{--}1.0$ | trains/min | Post-match clearing timetable dispatch frequency.

 |
| **Train Crush Capacity** | $Cap_{\text{train}}$ | $1,000$ | pax/train | Nominal capacity for standard 6-car EMU sets.

 |
| **Nominal Outflow Rate** | $C_{\text{base}, i}$ | $f_{\text{train}, i} \cdot Cap_{\text{train}}$ | pax/min | Baseline clearance rate prior to congestion degradation.

 |

---

## 4. Fanbase Corridor Demographics

Attendance assumes 30,000 fans per team (increasing by +10,000 fans per side during major rivalry fixtures). For each team, 18% of supporters are uniformly distributed across the three corridors (6% each), while the remaining 82% follow traditional home heartland demographic lines:

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

## 5. Implementation Architecture

The core mechanics verified in Python:

```python
from math import exp, gamma as gamma_function, isfinite

def _require_positive(value: float, name: str) -> None:
    if not isfinite(value) or value <= 0:
        raise ValueError(f"{name} must be a finite positive number")

def gamma_density(time: float, shape: float, scale: float) -> float:
    _require_positive(time, "time")
    _require_positive(shape, "shape")
    _require_positive(scale, "scale")
    return (time ** (shape - 1) * exp(-time / scale)) / (scale**shape * gamma_function(shape))

def queue_throttle(queue: float, gate_capacity: float, sensitivity: float) -> float:
    _require_positive(gate_capacity, "gate_capacity")
    _require_positive(sensitivity, "sensitivity")
    if queue < 0:
        raise ValueError("queue must not be negative")
    arg = sensitivity * (queue - gate_capacity)
    return exp(-arg) / (1.0 + exp(-arg)) if arg >= 0 else 1.0 / (1.0 + exp(arg))

def boarding_efficiency(density: float, d_crit: float, gamma_exp: float, eta_min: float) -> float:
    return eta_min + (1.0 - eta_min) / (1.0 + (density / d_crit) ** gamma_exp)

```
---

## References

1. Hohman, F. (2018). *CSE 6730 Modeling and Simulation Project Report*. [FRED HOHMAN CSE 6730](https://fredhohman.com/assets/cse-6730-project-report.pdf)

2. Ritz, C., et al. (2019). *Dose-Response Analysis Using R*. *The R Journal*, 11(2). [RJ-2019-003](https://journal.r-project.org/articles/RJ-2019-003/RJ-2019-003.pdf)