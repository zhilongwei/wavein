"""Generate validated PETSc YAML inputs for the WaveIn wave-tank apps."""

import argparse
import json
import math
import shlex
from collections.abc import Sequence
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path

GRAVITY = 9.80665
DEFAULT_CELLS_PER_WAVELENGTH = 20
DEFAULT_CELLS_PER_DEPTH = 10
DEFAULT_CFL = 0.5
DEFAULT_INLET_WAVELENGTHS = 1.0
DEFAULT_OUTLET_WAVELENGTHS = 2.0
DEFAULT_INLET_BUFFER_WAVELENGTHS = 1.0
DEFAULT_OUTLET_BUFFER_WAVELENGTHS = 1.0
DEFAULT_FRAMES_PER_PERIOD = 10
DEFAULT_OUTLET_RELAXATION_STRENGTH = 10.0
DEFAULT_RAMP_UP_PERIODS = Fraction(2)
DEFAULT_MINIMUM_FREQUENCY_RATIO = 1.0 / (2.0 * math.pi * 0.834)
DEFAULT_MAXIMUM_FREQUENCY_RATIO = 20.0 / (2.0 * math.pi * 0.834)


@dataclass(frozen=True)
class WaveProperties:
    wavelength: float
    celerity: float


@dataclass(frozen=True)
class Grid:
    requested_tank_length: float
    tank_wavelengths: int
    tank_length: float
    inlet_buffer_length: float
    outlet_buffer_length: float
    xmin: float
    xmax: float
    nx: int
    nz: int
    dx: float
    dz: float


@dataclass(frozen=True)
class FlowSchedule:
    simulation_periods: Fraction
    solution_record_periods: Fraction
    surface_elevation_record_periods: Fraction
    solution_record_start_time: float
    surface_elevation_record_start_time: float
    end_time: float
    interval: float
    dt: float
    steps_per_period: int
    num_steps: int
    output_stride: int
    solution_output_count: int
    surface_elevation_output_count: int


@dataclass(frozen=True)
class GeneratedInputs:
    app_name: str
    input_file: Path
    wave: WaveProperties
    grid: Grid
    schedule: FlowSchedule
    dt_limit: float
    options: tuple[tuple[str, object], ...]
    metadata: tuple[tuple[str, object], ...] = ()


def finite_positive(value: str) -> float:
    number = float(value)
    if not math.isfinite(number) or number <= 0.0:
        raise argparse.ArgumentTypeError(
            f"expected a positive finite number, got {value!r}"
        )
    return number


def finite_nonnegative(value: str) -> float:
    number = float(value)
    if not math.isfinite(number) or number < 0.0:
        raise argparse.ArgumentTypeError(
            f"expected a nonnegative finite number, got {value!r}"
        )
    return number


def finite_at_least_one(value: str) -> float:
    number = float(value)
    if not math.isfinite(number) or number < 1.0:
        raise argparse.ArgumentTypeError(
            f"expected a finite number at least 1, got {value!r}"
        )
    return number


def cfl_number(value: str) -> float:
    number = finite_positive(value)
    if number > 1.0:
        raise argparse.ArgumentTypeError(
            f"expected a CFL number no greater than 1, got {value!r}"
        )
    return number


def positive_integer(value: str) -> int:
    number = int(value)
    if number < 1:
        raise argparse.ArgumentTypeError(f"expected a positive integer, got {value!r}")
    return number


def nonnegative_integer(value: str) -> int:
    number = int(value)
    if number < 0:
        raise argparse.ArgumentTypeError(
            f"expected a nonnegative integer, got {value!r}"
        )
    return number


def nonnegative_multiple(value: str) -> Fraction:
    try:
        multiple = Fraction(value)
    except (ValueError, ZeroDivisionError) as error:
        raise argparse.ArgumentTypeError(
            f"expected a period multiple, got {value!r}"
        ) from error
    if multiple < 0:
        raise argparse.ArgumentTypeError(
            f"expected a nonnegative period multiple, got {value!r}"
        )
    return multiple


def positive_multiple(value: str) -> Fraction:
    multiple = nonnegative_multiple(value)
    if multiple == 0:
        raise argparse.ArgumentTypeError(
            f"expected a positive period multiple, got {value!r}"
        )
    return multiple


def linear_wave_properties(water_depth: float, period: float) -> WaveProperties:
    """Solve omega^2 = g k tanh(k h) using the same constants as AiryWave."""
    omega = 2.0 * math.pi / period
    nondimensional_frequency = omega * math.sqrt(water_depth / GRAVITY)
    kh = max(nondimensional_frequency, nondimensional_frequency**2)

    for _ in range(50):
        tanh_kh = math.tanh(kh)
        residual = kh * tanh_kh - nondimensional_frequency**2
        if abs(residual) <= 1.0e-12:
            break
        derivative = tanh_kh + kh * (1.0 - tanh_kh**2)
        kh -= residual / derivative
    else:
        raise ValueError("linear-wave dispersion solver did not converge")

    wavenumber = kh / water_depth
    wavelength = 2.0 * math.pi / wavenumber
    return WaveProperties(wavelength=wavelength, celerity=wavelength / period)


def compute_grid(
    water_depth: float,
    requested_tank_length: float,
    wave: WaveProperties,
    cells_per_wavelength: int,
    cells_per_depth: int,
    inlet_wavelengths: float,
    outlet_wavelengths: float,
    inlet_buffer_wavelengths: float,
    outlet_buffer_wavelengths: float,
) -> Grid:
    wavelength_ratio = requested_tank_length / wave.wavelength
    tank_wavelengths = max(
        1, math.ceil(wavelength_ratio - 1.0e-9 * max(1.0, wavelength_ratio))
    )
    tank_length = tank_wavelengths * wave.wavelength
    inlet_buffer_length = inlet_buffer_wavelengths * wave.wavelength
    outlet_buffer_length = outlet_buffer_wavelengths * wave.wavelength
    xmin = -(inlet_wavelengths + inlet_buffer_wavelengths) * wave.wavelength
    xmax = (
        tank_length + (outlet_buffer_wavelengths + outlet_wavelengths) * wave.wavelength
    )
    computational_wavelengths = (
        inlet_wavelengths
        + inlet_buffer_wavelengths
        + tank_wavelengths
        + outlet_buffer_wavelengths
        + outlet_wavelengths
    )
    target_nx = computational_wavelengths * cells_per_wavelength
    nx = max(2, math.ceil(target_nx - 1.0e-12 * max(1.0, target_nx)))
    dx = (xmax - xmin) / nx
    dz = water_depth / cells_per_depth
    return Grid(
        requested_tank_length=requested_tank_length,
        tank_wavelengths=tank_wavelengths,
        tank_length=tank_length,
        inlet_buffer_length=inlet_buffer_length,
        outlet_buffer_length=outlet_buffer_length,
        xmin=xmin,
        xmax=xmax,
        nx=nx,
        nz=cells_per_depth,
        dx=dx,
        dz=dz,
    )


def align_flow_schedule(
    period: float,
    simulation_periods: Fraction,
    solution_record_periods: Fraction,
    surface_elevation_record_periods: Fraction,
    frames_per_period: int,
    dt_limit: float,
) -> FlowSchedule:
    if solution_record_periods > simulation_periods:
        raise ValueError(
            "--solution-record-periods must not exceed --simulation-periods"
        )
    if surface_elevation_record_periods > simulation_periods:
        raise ValueError(
            "--surface-elevation-record-periods must not exceed --simulation-periods"
        )

    solution_output_intervals = solution_record_periods * frames_per_period
    if solution_output_intervals.denominator != 1:
        raise ValueError(
            "--solution-record-periods times --frames-per-period must be an integer"
        )
    surface_elevation_output_intervals = (
        surface_elevation_record_periods * frames_per_period
    )
    if surface_elevation_output_intervals.denominator != 1:
        raise ValueError(
            "--surface-elevation-record-periods times --frames-per-period must be an integer"
        )

    solution_start_periods = simulation_periods - solution_record_periods
    surface_elevation_start_periods = (
        simulation_periods - surface_elevation_record_periods
    )
    alignment = math.lcm(
        simulation_periods.denominator,
        solution_start_periods.denominator,
        surface_elevation_start_periods.denominator,
        frames_per_period,
    )
    minimum_steps_per_period = math.ceil(period / dt_limit)
    steps_per_period = (
        (minimum_steps_per_period + alignment - 1) // alignment
    ) * alignment

    num_steps_fraction = simulation_periods * steps_per_period
    solution_start_step_fraction = solution_start_periods * steps_per_period
    surface_elevation_start_step_fraction = (
        surface_elevation_start_periods * steps_per_period
    )
    if (
        num_steps_fraction.denominator != 1
        or solution_start_step_fraction.denominator != 1
        or surface_elevation_start_step_fraction.denominator != 1
    ):
        raise ValueError("could not construct an integer-step schedule")

    dt = period / steps_per_period
    output_stride = steps_per_period // frames_per_period
    return FlowSchedule(
        simulation_periods=simulation_periods,
        solution_record_periods=solution_record_periods,
        surface_elevation_record_periods=surface_elevation_record_periods,
        solution_record_start_time=float(solution_start_periods) * period,
        surface_elevation_record_start_time=float(surface_elevation_start_periods)
        * period,
        end_time=float(simulation_periods) * period,
        interval=period / frames_per_period,
        dt=dt,
        steps_per_period=steps_per_period,
        num_steps=num_steps_fraction.numerator,
        output_stride=output_stride,
        solution_output_count=solution_output_intervals.numerator + 1,
        surface_elevation_output_count=surface_elevation_output_intervals.numerator + 1,
    )


def validate_common_inputs(args: argparse.Namespace) -> None:
    if args.cells_per_depth < 2:
        raise ValueError("--cells-per-depth must be at least 2")
    if args.cells_per_wavelength < 2:
        raise ValueError("--cells-per-wavelength must be at least 2")
    output_path = Path(args.output)
    if not output_path.name:
        raise ValueError("--output must name an HDF5 file")
    if not output_path.parent.exists():
        raise ValueError(f"output directory does not exist: {output_path.parent}")
    if args.input_file.resolve() == output_path.resolve():
        raise ValueError("--input-file and --output must refer to different files")


def common_generated_inputs(
    args: argparse.Namespace, period: float
) -> tuple[WaveProperties, Grid, FlowSchedule, float, float, float]:
    validate_common_inputs(args)
    wave = linear_wave_properties(args.water_depth, period)
    grid = compute_grid(
        water_depth=args.water_depth,
        requested_tank_length=args.tank_length,
        wave=wave,
        cells_per_wavelength=args.cells_per_wavelength,
        cells_per_depth=args.cells_per_depth,
        inlet_wavelengths=args.inlet_wavelengths,
        outlet_wavelengths=args.outlet_wavelengths,
        inlet_buffer_wavelengths=args.inlet_buffer_wavelengths,
        outlet_buffer_wavelengths=args.outlet_buffer_wavelengths,
    )
    dt_limit = args.cfl * min(grid.dx, grid.dz) / wave.celerity
    solution_record_periods = (
        args.solution_record_periods
        if args.solution_record_periods is not None
        else args.simulation_periods
    )
    surface_elevation_record_periods = (
        args.surface_elevation_record_periods
        if args.surface_elevation_record_periods is not None
        else args.simulation_periods
    )
    schedule = align_flow_schedule(
        period=period,
        simulation_periods=args.simulation_periods,
        solution_record_periods=solution_record_periods,
        surface_elevation_record_periods=surface_elevation_record_periods,
        frames_per_period=args.frames_per_period,
        dt_limit=dt_limit,
    )
    relaxation_rate = args.outlet_relaxation_strength / period
    ramp_up_time = float(args.ramp_up_periods) * period
    return wave, grid, schedule, dt_limit, relaxation_rate, ramp_up_time


def common_options(
    args: argparse.Namespace,
    grid: Grid,
    schedule: FlowSchedule,
    relaxation_rate: float,
    ramp_up_time: float,
) -> list[tuple[str, object]]:
    return [
        ("sim_start_time", 0.0),
        ("sim_end_time", schedule.end_time),
        ("sim_dt", schedule.dt),
        ("flow_field_output_start_time", schedule.solution_record_start_time),
        ("flow_field_output_end_time", schedule.end_time),
        (
            "flow_surface_elevation_output_start_time",
            schedule.surface_elevation_record_start_time,
        ),
        ("flow_surface_elevation_output_end_time", schedule.end_time),
        ("flow_output_interval", schedule.interval),
        ("output", args.output),
        ("xmin", grid.xmin),
        ("xmax", grid.xmax),
        ("nx", grid.nx),
        ("nz", grid.nz),
        ("nin", args.inlet_wavelengths),
        ("nout", args.outlet_wavelengths),
        ("gamma", relaxation_rate),
        ("ramp_up_time", ramp_up_time),
    ]


def common_metadata(args: argparse.Namespace) -> tuple[tuple[str, object], ...]:
    return (
        ("cells_per_wavelength", args.cells_per_wavelength),
        ("cells_per_depth", args.cells_per_depth),
        ("cfl", args.cfl),
        ("inlet_wavelengths", args.inlet_wavelengths),
        ("outlet_wavelengths", args.outlet_wavelengths),
        ("inlet_buffer_wavelengths", args.inlet_buffer_wavelengths),
        ("outlet_buffer_wavelengths", args.outlet_buffer_wavelengths),
        ("frames_per_period", args.frames_per_period),
        ("outlet_relaxation_strength", args.outlet_relaxation_strength),
        ("ramp_up_periods", args.ramp_up_periods),
    )


def generate_regular_inputs(args: argparse.Namespace) -> GeneratedInputs:
    wave, grid, schedule, dt_limit, relaxation_rate, ramp_up_time = (
        common_generated_inputs(args, args.wave_period)
    )
    options = common_options(args, grid, schedule, relaxation_rate, ramp_up_time)
    options[9:9] = [
        ("wave_height", args.wave_height),
        ("wave_period", args.wave_period),
        ("water_depth", args.water_depth),
    ]
    return GeneratedInputs(
        app_name="regular_wave_tank",
        input_file=args.input_file,
        wave=wave,
        grid=grid,
        schedule=schedule,
        dt_limit=dt_limit,
        options=tuple(options),
        metadata=common_metadata(args),
    )


def generate_irregular_inputs(args: argparse.Namespace) -> GeneratedInputs:
    wave, grid, schedule, dt_limit, relaxation_rate, ramp_up_time = (
        common_generated_inputs(args, args.peak_period)
    )
    omega_peak = 2.0 * math.pi / args.peak_period
    omega_min = args.minimum_frequency_ratio * omega_peak
    omega_max = args.maximum_frequency_ratio * omega_peak
    if omega_max <= omega_min:
        raise ValueError(
            "--maximum-frequency-ratio must be greater than --minimum-frequency-ratio"
        )

    options = common_options(args, grid, schedule, relaxation_rate, ramp_up_time)
    options[9:9] = [
        ("significant_wave_height", args.significant_wave_height),
        ("peak_period", args.peak_period),
        ("peak_enhancement_factor", args.peak_enhancement_factor),
        ("water_depth", args.water_depth),
        ("omega_min", omega_min),
        ("omega_max", omega_max),
        ("component_count", args.component_count),
        ("random_seed", args.random_seed),
    ]
    return GeneratedInputs(
        app_name="irregular_wave_tank",
        input_file=args.input_file,
        wave=wave,
        grid=grid,
        schedule=schedule,
        dt_limit=dt_limit,
        options=tuple(options),
        metadata=(
            *common_metadata(args),
            ("minimum_frequency_ratio", args.minimum_frequency_ratio),
            ("maximum_frequency_ratio", args.maximum_frequency_ratio),
            ("omega_min", omega_min),
            ("omega_max", omega_max),
        ),
    )


def format_value(value: object) -> str:
    if isinstance(value, Fraction):
        return str(value.numerator) if value.denominator == 1 else f"{float(value):.8e}"
    if isinstance(value, float):
        return f"{value:.8e}"
    return str(value)


def yaml_scalar(value: object) -> str:
    if isinstance(value, str):
        return json.dumps(value)
    return format_value(value)


def render_yaml(generated: GeneratedInputs) -> str:
    metadata = [
        ("requested_tank_length", generated.grid.requested_tank_length),
        ("tank_wavelengths", generated.grid.tank_wavelengths),
        ("tank_length", generated.grid.tank_length),
        ("domain_of_interest_xmin", 0.0),
        ("domain_of_interest_xmax", generated.grid.tank_length),
        ("inlet_buffer_length", generated.grid.inlet_buffer_length),
        ("outlet_buffer_length", generated.grid.outlet_buffer_length),
        ("wavelength", generated.wave.wavelength),
        ("celerity", generated.wave.celerity),
        ("dx", generated.grid.dx),
        ("dz", generated.grid.dz),
        ("dt_cfl", generated.dt_limit),
        ("simulation_periods", generated.schedule.simulation_periods),
        ("solution_record_periods", generated.schedule.solution_record_periods),
        (
            "surface_elevation_record_periods",
            generated.schedule.surface_elevation_record_periods,
        ),
        ("steps_per_period", generated.schedule.steps_per_period),
        ("num_steps", generated.schedule.num_steps),
        ("output_stride", generated.schedule.output_stride),
        ("solution_output_count", generated.schedule.solution_output_count),
        (
            "surface_elevation_output_count",
            generated.schedule.surface_elevation_output_count,
        ),
        *generated.metadata,
    ]
    options = dict(generated.options)
    lines = [
        "# Generated by scripts/setup_wave_tank.py",
        f"# app: {generated.app_name}",
        *[f"# {name}: {format_value(value)}" for name, value in metadata],
        "",
        "sim:",
        f"  start_time: {yaml_scalar(options.pop('sim_start_time'))}",
        f"  end_time: {yaml_scalar(options.pop('sim_end_time'))}",
        f"  dt: {yaml_scalar(options.pop('sim_dt'))}",
        "",
        *[f"{name}: {yaml_scalar(value)}" for name, value in options.items()],
        "",
    ]
    return "\n".join(lines)


def write_input(generated: GeneratedInputs) -> None:
    if generated.input_file.suffix.lower() not in {".yaml", ".yml"}:
        raise ValueError("--input-file must end in .yaml or .yml")
    parent = generated.input_file.parent
    if not parent.exists():
        raise ValueError(f"input-file directory does not exist: {parent}")
    if not parent.is_dir():
        raise ValueError(f"input-file parent is not a directory: {parent}")
    if generated.input_file.is_dir():
        raise ValueError(f"input-file path is a directory: {generated.input_file}")
    generated.input_file.write_text(render_yaml(generated), encoding="utf-8")


def print_summary(generated: GeneratedInputs) -> None:
    print(f"Wrote {generated.input_file}")
    print(
        f"wave: wavelength={generated.wave.wavelength:.8e} m, "
        f"celerity={generated.wave.celerity:.8e} m/s"
    )
    print(
        f"domain of interest: requested={generated.grid.requested_tank_length:.8e} m, "
        f"actual={generated.grid.tank_length:.8e} m "
        f"({generated.grid.tank_wavelengths} wavelengths)"
    )
    print(
        f"domain: xmin={generated.grid.xmin:.8e} m, xmax={generated.grid.xmax:.8e} m, "
        f"nx={generated.grid.nx}, nz={generated.grid.nz}"
    )
    print(
        f"buffers: inlet={generated.grid.inlet_buffer_length:.8e} m, "
        f"outlet={generated.grid.outlet_buffer_length:.8e} m"
    )
    print(f"grid: dx={generated.grid.dx:.8e} m, dz={generated.grid.dz:.8e} m")
    print(
        f"schedule: periods={format_value(generated.schedule.simulation_periods)}, "
        f"steps/period={generated.schedule.steps_per_period}, "
        f"dt={generated.schedule.dt:.8e} s (limit={generated.dt_limit:.8e} s)"
    )
    print(f"output interval: {generated.schedule.interval:.8e} s")
    print(
        "solution output: final "
        f"{format_value(generated.schedule.solution_record_periods)} periods, "
        f"frames={generated.schedule.solution_output_count}"
    )
    print(
        "surface-elevation output: final "
        f"{format_value(generated.schedule.surface_elevation_record_periods)} periods, "
        f"frames={generated.schedule.surface_elevation_output_count}"
    )
    print(
        "run: "
        f"./build/release/src/apps/{generated.app_name} "
        f"-options_file_yaml {shlex.quote(str(generated.input_file))}"
    )


def add_common_arguments(parser: argparse.ArgumentParser) -> None:
    physics = parser.add_argument_group("tank and grid")
    physics.add_argument(
        "--water-depth", required=True, type=finite_positive, help="water depth h [m]"
    )
    physics.add_argument(
        "--tank-length",
        "--domain-length-of-interest",
        dest="tank_length",
        required=True,
        type=finite_positive,
        help="minimum domain length of interest [m], rounded up to a whole wavelength",
    )
    physics.add_argument(
        "--cells-per-wavelength",
        type=positive_integer,
        default=DEFAULT_CELLS_PER_WAVELENGTH,
        help=f"horizontal cells per wavelength (default: {DEFAULT_CELLS_PER_WAVELENGTH})",
    )
    physics.add_argument(
        "--cells-per-depth",
        type=positive_integer,
        default=DEFAULT_CELLS_PER_DEPTH,
        help=f"vertical cells per water depth (default: {DEFAULT_CELLS_PER_DEPTH})",
    )
    physics.add_argument(
        "--cfl",
        type=cfl_number,
        default=DEFAULT_CFL,
        help=f"CFL number (default: {format_value(DEFAULT_CFL)})",
    )
    physics.add_argument(
        "--inlet-wavelengths",
        type=finite_positive,
        default=DEFAULT_INLET_WAVELENGTHS,
        help=(
            "inlet-zone length in wavelengths "
            f"(default: {format_value(DEFAULT_INLET_WAVELENGTHS)})"
        ),
    )
    physics.add_argument(
        "--outlet-wavelengths",
        type=finite_positive,
        default=DEFAULT_OUTLET_WAVELENGTHS,
        help=(
            "outlet-zone length in wavelengths "
            f"(default: {format_value(DEFAULT_OUTLET_WAVELENGTHS)})"
        ),
    )
    physics.add_argument(
        "--inlet-buffer-wavelengths",
        type=finite_nonnegative,
        default=DEFAULT_INLET_BUFFER_WAVELENGTHS,
        help=(
            "distance from the inlet zone to the domain of interest in wavelengths "
            f"(default: {format_value(DEFAULT_INLET_BUFFER_WAVELENGTHS)})"
        ),
    )
    physics.add_argument(
        "--outlet-buffer-wavelengths",
        type=finite_nonnegative,
        default=DEFAULT_OUTLET_BUFFER_WAVELENGTHS,
        help=(
            "distance from the domain of interest to the outlet zone in wavelengths "
            f"(default: {format_value(DEFAULT_OUTLET_BUFFER_WAVELENGTHS)})"
        ),
    )

    schedule = parser.add_argument_group("dimensionless simulation and output schedule")
    schedule.add_argument(
        "--simulation-periods",
        required=True,
        type=positive_multiple,
        help="simulation duration as a multiple of T or Tp",
    )
    schedule.add_argument(
        "--solution-record-periods",
        type=nonnegative_multiple,
        help="number of final periods for solution output (default: entire simulation)",
    )
    schedule.add_argument(
        "--surface-elevation-record-periods",
        type=nonnegative_multiple,
        help=(
            "number of final periods for surface-elevation output "
            "(default: entire simulation)"
        ),
    )
    schedule.add_argument(
        "--frames-per-period",
        type=positive_integer,
        default=DEFAULT_FRAMES_PER_PERIOD,
        help=(
            "recording rate for both outputs in frames per period "
            f"(default: {DEFAULT_FRAMES_PER_PERIOD})"
        ),
    )

    forcing = parser.add_argument_group("dimensionless wavemaker settings")
    forcing.add_argument(
        "--outlet-relaxation-strength",
        type=finite_positive,
        default=DEFAULT_OUTLET_RELAXATION_STRENGTH,
        help=(
            "dimensionless gamma*T "
            f"(default: {format_value(DEFAULT_OUTLET_RELAXATION_STRENGTH)})"
        ),
    )
    forcing.add_argument(
        "--ramp-up-periods",
        type=positive_multiple,
        default=DEFAULT_RAMP_UP_PERIODS,
        help=(
            "linear ramp duration in periods "
            f"(default: {format_value(DEFAULT_RAMP_UP_PERIODS)})"
        ),
    )

    files = parser.add_argument_group("files")
    files.add_argument(
        "--input-file", required=True, type=Path, help="PETSc YAML input file to write"
    )
    files.add_argument(
        "--output", required=True, help="HDF5 output filename passed to the app"
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Validate wave-tank inputs and write a PETSc YAML file."
    )
    subparsers = parser.add_subparsers(dest="wave_type", required=True)

    regular = subparsers.add_parser("regular", help="configure the regular-wave tank")
    add_common_arguments(regular)
    regular.add_argument(
        "--wave-period", required=True, type=finite_positive, help="wave period T [s]"
    )
    regular.add_argument(
        "--wave-height", required=True, type=finite_positive, help="wave height H [m]"
    )
    regular.set_defaults(generator=generate_regular_inputs)

    irregular = subparsers.add_parser(
        "irregular", help="configure the JONSWAP irregular-wave tank"
    )
    add_common_arguments(irregular)
    irregular.add_argument(
        "--peak-period", required=True, type=finite_positive, help="peak period Tp [s]"
    )
    irregular.add_argument(
        "--significant-wave-height",
        required=True,
        type=finite_positive,
        help="significant wave height Hs [m]",
    )
    irregular.add_argument(
        "--peak-enhancement-factor",
        type=finite_at_least_one,
        default=3.3,
        help=f"JONSWAP peak-enhancement factor (default: {format_value(3.3)})",
    )
    irregular.add_argument(
        "--component-count",
        type=positive_integer,
        default=100,
        help="number of wave components (default: 100)",
    )
    irregular.add_argument(
        "--random-seed",
        type=nonnegative_integer,
        default=0,
        help="random phase seed (default: 0)",
    )
    irregular.add_argument(
        "--minimum-frequency-ratio",
        type=finite_positive,
        default=DEFAULT_MINIMUM_FREQUENCY_RATIO,
        help="minimum frequency divided by peak frequency",
    )
    irregular.add_argument(
        "--maximum-frequency-ratio",
        type=finite_positive,
        default=DEFAULT_MAXIMUM_FREQUENCY_RATIO,
        help="maximum frequency divided by peak frequency",
    )
    irregular.set_defaults(generator=generate_irregular_inputs)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        generated = args.generator(args)
        write_input(generated)
    except ValueError as error:
        parser.error(str(error))
    print_summary(generated)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
