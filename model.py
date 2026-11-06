import numpy as np
import matplotlib.pyplot as plt
from math import exp, gamma as gamma_function, isfinite
from scipy.stats import gamma as gamma_dist

# ---------------------------------------
# CONSTANTS & PARAMETERS
# ---------------------------------------

STATION_NAMES = ["J", "R", "F"]

STATIONS = {
    "J": {"tau": 7.0,  "gate_crit": 4000.0,  "area": 1800.0,  "f_train": 0.4},
    "R": {"tau": 9.0,  "gate_crit": 18000.0, "area": 9000.0,  "f_train": 0.8},
    "F": {"tau": 18.0, "gate_crit": 25000.0, "area": 15000.0, "f_train": 1.0},
}

MODEL = {
    "shape": 4.0,
    "beta_q": 0.002,
    "d_crit": 2.25,
    "gamma_exp": 4.0,
    "eta_min": 0.35,
    "train_capacity": 1000.0,
}

# corridor shares (alpha_J, alpha_R, alpha_F), essendon corrected
TEAM_SHARES = {
    "Collingwood": (0.70, 0.15, 0.15),
    "Hawthorn":    (0.10, 0.70, 0.20),
    "Melbourne":   (0.20, 0.30, 0.50),
    "Richmond":    (0.22, 0.72, 0.06),
    "Essendon":    (0.06, 0.68, 0.26),
    "Carlton":     (0.16, 0.16, 0.68),
    "St Kilda":    (0.06, 0.68, 0.26),
    "Geelong":     (0.06, 0.33, 0.61),
}

COLORS = {"J": "tab:blue", "R": "tab:orange", "F": "tab:green"}
LABELS = {"J": "Jolimont", "R": "Richmond", "F": "Flinders"}

# ---------------------------------------
# FLOW FUNCTIONS
# ---------------------------------------

def _require_positive(value: float, name: str) -> None:
    if not isfinite(value) or value <= 0:
        raise ValueError(f"{name} must be a finite positive number")

def gamma_density(time: float, shape: float, scale: float) -> float:
    _require_positive(shape, "shape")
    _require_positive(scale, "scale")
    if time < 0:  # nobody arrives before the final siren
        return 0.0
    return time ** (shape - 1) * exp(-time / scale) / (scale**shape * gamma_function(shape))

def gamma_hazard(time: float, shape: float, scale: float) -> float:
    # f / (1 - F): arrival rate per person still outside, so throttled fans keep trying
    _require_positive(shape, "shape")
    _require_positive(scale, "scale")
    if time <= 0:
        return 0.0
    return exp(gamma_dist.logpdf(time, shape, scale=scale) - gamma_dist.logsf(time, shape, scale=scale))

def queue_throttle(queue: float, gate_capacity: float, sensitivity: float) -> float:
    _require_positive(gate_capacity, "gate_capacity")
    _require_positive(sensitivity, "sensitivity")
    if queue < 0:  # cannot use _require_positive since queue starts AT zero
        raise ValueError("queue must not be negative")
    argument = sensitivity * (queue - gate_capacity)  # split exp forms to dodge overflow
    if argument >= 0:
        return exp(-argument) / (1.0 + exp(-argument))
    return 1.0 / (1.0 + exp(argument))

def boarding_efficiency(density: float, d_crit: float, gamma_exp: float, eta_min: float) -> float:
    _require_positive(d_crit, "d_crit")
    _require_positive(gamma_exp, "gamma_exp")
    if not 0 < eta_min <= 1:
        raise ValueError("eta_min must be greater than 0 and at most 1")
    if density < 0:
        raise ValueError("density must not be negative")
    return eta_min + (1.0 - eta_min) / (1.0 + (density / d_crit) ** gamma_exp)

def inflow_rate(time, queue, remaining, shape, scale, gate_capacity, sensitivity):
    if remaining < 0:
        raise ValueError("remaining must not be negative")
    return remaining * gamma_hazard(time, shape, scale) * queue_throttle(queue, gate_capacity, sensitivity)

def outflow_rate(queue, platform_area, base_capacity, d_crit, gamma_exp, eta_min):
    _require_positive(platform_area, "platform_area")
    _require_positive(base_capacity, "base_capacity")
    if queue < 0:
        raise ValueError("queue must not be negative")
    density = queue / platform_area
    return base_capacity * boarding_efficiency(density, d_crit, gamma_exp, eta_min)

# ---------------------------------------
# SIMULATION
# ---------------------------------------

def match_allocation(home_team, away_team, attendance_each=30000.0, rivalry=False):
    _require_positive(attendance_each, "attendance_each")
    if home_team not in TEAM_SHARES or away_team not in TEAM_SHARES:
        raise ValueError("unrecognised team name")
    each = attendance_each + 10000.0 if rivalry else attendance_each  # rivalry fixtures draw extra fans
    total = 2.0 * each
    home_alpha = TEAM_SHARES[home_team]
    away_alpha = TEAM_SHARES[away_team]
    shares = {}
    for idx, name in enumerate(STATION_NAMES):
        shares[name] = (each * home_alpha[idx] + each * away_alpha[idx]) / total  # attendance-weighted blend
    return total, shares

def simulate_egress(attendance, shares, stations=STATIONS, model=MODEL, duration=90.0, dt=0.1):
    _require_positive(attendance, "attendance")
    _require_positive(duration, "duration")
    _require_positive(dt, "dt")
    if abs(sum(shares.values()) - 1.0) > 1e-9:
        raise ValueError("shares must sum to 1")

    n_steps = int(round(duration / dt))
    n_stations = len(STATION_NAMES)
    times = np.linspace(0.0, duration, n_steps + 1)
    M = np.zeros((n_steps + 1, n_stations))  # remaining allocation
    Q = np.zeros((n_steps + 1, n_stations))  # platform/gate queue
    C = np.zeros((n_steps + 1, n_stations))  # cleared

    for idx, name in enumerate(STATION_NAMES):
        M[0, idx] = attendance * shares[name]

    for step in range(n_steps):
        t = times[step]
        for idx, name in enumerate(STATION_NAMES):
            station = stations[name]
            scale = station["tau"] / model["shape"]
            base_capacity = station["f_train"] * model["train_capacity"]
            inflow = inflow_rate(t, Q[step, idx], M[step, idx], model["shape"], scale,
                                  station["gate_crit"], model["beta_q"])
            outflow = outflow_rate(Q[step, idx], station["area"], base_capacity,
                                    model["d_crit"], model["gamma_exp"], model["eta_min"])
            incoming = min(M[step, idx], inflow * dt)  # can't move more pax than exist
            departing = min(Q[step, idx] + incoming, outflow * dt)
            M[step + 1, idx] = M[step, idx] - incoming
            Q[step + 1, idx] = Q[step, idx] + incoming - departing
            C[step + 1, idx] = C[step, idx] + departing

    return times, M, Q, C

def fixture_clearance(home_team, away_team, station_overrides=None):
    # fraction cleared by t = 90 with optional {station: {param: value}} overrides
    stations = {name: {**params, **(station_overrides or {}).get(name, {})} for name, params in STATIONS.items()}
    attendance, shares = match_allocation(home_team, away_team)
    _, _, _, C = simulate_egress(attendance, shares, stations=stations)
    return C[-1].sum() / attendance

def required_frequency(home_team, away_team, station="J", target=0.999, step=0.05, max_frequency=2.0):
    # smallest dispatch rate (trains/min, on a `step` grid) at `station` that clears `target` of the crowd; nan if none
    frequency = STATIONS[station]["f_train"]
    while frequency <= max_frequency + 1e-9:
        if fixture_clearance(home_team, away_team, {station: {"f_train": frequency}}) >= target:
            return round(frequency, 2)
        frequency += step
    return np.nan

def time_to_clear(times, C, attendance, fraction=0.95):
    # first time the cumulative cleared total reaches `fraction` of attendance; nan if never
    reached = np.nonzero(C.sum(axis=1) >= fraction * attendance)[0]
    return times[reached[0]] if reached.size else np.nan

# ---------------------------------------
# VISUALIZATION
# ---------------------------------------

def plot_state_trajectories(times, M, Q, C, fixture):
    # one figure per state so the report can break them across pages
    titles = ["Remaining allocation $M_i(t)$", "Platform/gate queue $Q_i(t)$", "Cumulative cleared $C_i(t)$"]
    figs = []
    for arr, title in zip([M, Q, C], titles):
        fig, ax = plt.subplots(figsize=(7.5, 3))
        for idx, name in enumerate(STATION_NAMES):
            ax.plot(times, arr[:, idx], color=COLORS[name], label=LABELS[name])
        ax.set_xlabel("time since final siren (min)")
        ax.set_ylabel("pax")
        ax.set_title(f"{title}, {fixture}")
        ax.grid(alpha=0.3)
        ax.legend(loc="upper right")
        fig.tight_layout()
        figs.append(fig)
    return figs

def plot_clearance_vs_jolimont(j_shares, fractions):
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    ax.scatter(j_shares * 100, fractions * 100, color="tab:red")
    ax.set_xlabel(r"Jolimont demographic share $\alpha_J$ (%)")
    ax.set_ylabel("total cleared by t = 90 min (%)")
    ax.set_title(f"Clearance fraction vs. Jolimont corridor loading ({len(fractions)} fixtures)")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    return fig
