"""MCG post-match rail egress model.

The model follows the three conservation states described in ``README.md``:
remaining passengers (``M``), station/platform queue (``Q``), and cleared
passengers (``C``).  All rates are expressed per minute and all passenger
counts are people.
"""

from math import ceil, exp, gamma as gamma_function, isfinite


DEFAULT_MODEL = {
	"shape": 4.0,
	"queue_sensitivity": 0.002,
	"critical_density": 2.25,
	"density_exponent": 4.0,
	"minimum_efficiency": 0.35,
}


def _require_positive(value: float, name: str) -> None:
	if not isfinite(value) or value <= 0:
		raise ValueError(f"{name} must be a finite positive number")


def _validate_model_parameters(parameters: dict[str, float]) -> None:
	_require_positive(parameters["shape"], "shape")
	_require_positive(parameters["queue_sensitivity"], "queue_sensitivity")
	_require_positive(parameters["critical_density"], "critical_density")
	_require_positive(parameters["density_exponent"], "density_exponent")
	if not 0 < parameters["minimum_efficiency"] <= 1:
		raise ValueError("minimum_efficiency must be greater than 0 and at most 1")


def validate_station(station: dict[str, float]) -> None:
	"""Validate physical and service parameters for a station."""

	for value, name in (
		(station["walk_time"], "walk_time"),
		(station["gate_capacity"], "gate_capacity"),
		(station["platform_area"], "platform_area"),
		(station["train_frequency"], "train_frequency"),
		(station["train_capacity"], "train_capacity"),
	):
		_require_positive(value, name)


def gamma_density(time: float, shape: float, scale: float) -> float:
	"""Evaluate the Gamma surge density at ``time``.

	The density is zero before the match ends.  ``scale`` is normally
	``walk_time / shape``, as specified in the README.
	"""

	_require_positive(shape, "shape")
	_require_positive(scale, "scale")
	if time < 0:
		return 0.0
	return time ** (shape - 1) * exp(-time / scale) / (
		scale**shape * gamma_function(shape)
	)


def queue_throttle(queue: float, gate_capacity: float, sensitivity: float) -> float:
	"""Return the logistic inflow throttle ``Psi(Q)``."""

	if queue < 0:
		raise ValueError("queue must not be negative")
	_require_positive(gate_capacity, "gate_capacity")
	_require_positive(sensitivity, "sensitivity")
	argument = sensitivity * (queue - gate_capacity)
	if argument >= 0:
		return exp(-argument) / (1.0 + exp(-argument))
	return 1.0 / (1.0 + exp(argument))


def boarding_efficiency(
	queue: float,
	platform_area: float,
	critical_density: float,
	density_exponent: float,
	minimum_efficiency: float,
) -> float:
	"""Return boarding efficiency as a function of platform density."""

	if queue < 0:
		raise ValueError("queue must not be negative")
	_require_positive(platform_area, "platform_area")
	_require_positive(critical_density, "critical_density")
	_require_positive(density_exponent, "density_exponent")
	if not 0 < minimum_efficiency <= 1:
		raise ValueError("minimum_efficiency must be greater than 0 and at most 1")
	density = queue / platform_area
	return minimum_efficiency + (1.0 - minimum_efficiency) / (
		1.0 + (density / critical_density) ** density_exponent
	)


def nominal_outflow(station: dict[str, float]) -> float:
	"""Return ``f_train * Capacity_train`` for a station."""

	validate_station(station)
	return station["train_frequency"] * station["train_capacity"]


def station_rates(
	time: float,
	remaining: float,
	queue: float,
	station: dict[str, float],
	model: dict[str, float] | None = None,
) -> tuple[float, float]:
	"""Calculate station inflow and outflow at one point in time."""

	if remaining < 0 or queue < 0:
		raise ValueError("remaining and queue must not be negative")
	validate_station(station)
	model = DEFAULT_MODEL if model is None else model
	_validate_model_parameters(model)
	scale = station["walk_time"] / model["shape"]
	inflow = remaining * gamma_density(time, model["shape"], scale)
	inflow *= queue_throttle(queue, station["gate_capacity"], model["queue_sensitivity"])
	capacity = nominal_outflow(station) * boarding_efficiency(
		queue,
		station["platform_area"],
		model["critical_density"],
		model["density_exponent"],
		model["minimum_efficiency"],
	)
	return inflow, min(queue, capacity)


def state_derivatives(
	time: float,
	remaining: float,
	queue: float,
	station: dict[str, float],
	model: dict[str, float] | None = None,
) -> tuple[float, float, float]:
	"""Return ``(dM/dt, dQ/dt, dC/dt)`` for one station."""

	inflow, outflow = station_rates(time, remaining, queue, station, model)
	return -inflow, inflow - outflow, outflow


def simulate(
	attendance: float,
	stations: dict[str, dict[str, float]],
	shares: dict[str, float],
	model: dict[str, float] | None = None,
	duration: float = 90.0,
	step: float = 0.1,
) -> SimulationResult:
	"""Simulate post-match egress using bounded forward Euler integration.

	``shares`` gives each station's demographic share and must sum to one.
	``attendance`` is the total rail-bound attendance before applying those
	shares.
	"""

	if not isfinite(duration) or duration <= 0:
		raise ValueError("duration must be a finite positive number")
	if not isfinite(step) or step <= 0:
		raise ValueError("step must be a finite positive number")
	if set(stations) != set(shares):
		raise ValueError("stations and shares must contain the same station names")
	_require_positive(attendance, "attendance")
	model = DEFAULT_MODEL if model is None else model
	_validate_model_parameters(model)

	share_total = sum(shares.values())
	if not isfinite(share_total) or abs(share_total - 1.0) > 1e-9:
		raise ValueError("station shares must sum to 1")
	for name, station in stations.items():
		validate_station(station)
		if shares[name] < 0:
			raise ValueError("shares must not be negative")

	names = list(stations)
	times = [0.0]
	remaining = {name: [attendance * shares[name]] for name in names}
	queues = {name: [0.0] for name in names}
	cleared = {name: [0.0] for name in names}
	current_time = 0.0

	number_of_steps = ceil(duration / step)
	for index in range(number_of_steps):
		delta = min(step, duration - current_time)
		for name in names:
			station = stations[name]
			old_remaining = remaining[name][-1]
			old_queue = queues[name][-1]
			inflow, outflow = station_rates(
				current_time, old_remaining, old_queue, station, model
			)
			incoming = min(old_remaining, inflow * delta)
			departing = min(old_queue + incoming, outflow * delta)
			remaining[name].append(old_remaining - incoming)
			queues[name].append(old_queue + incoming - departing)
			cleared[name].append(cleared[name][-1] + departing)
		current_time = duration if index == number_of_steps - 1 else current_time + delta
		times.append(current_time)

	return {
		"times": times,
		"remaining": remaining,
		"queues": queues,
		"cleared": cleared,
		"total_cleared": sum(values[-1] for values in cleared.values()),
	}
