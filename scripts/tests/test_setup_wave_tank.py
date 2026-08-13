from __future__ import annotations

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import setup_wave_tank as setup


class SetupWaveTankTest(unittest.TestCase):
    def parse(self, *arguments: str):
        return setup.build_parser().parse_args(arguments)

    def test_regular_case_uses_period_based_schedule(self) -> None:
        args = self.parse(
            "regular",
            "--water-depth",
            "0.7",
            "--wave-period",
            "1.5",
            "--wave-height",
            "0.07",
            "--tank-length",
            "21.0",
            "--cells-per-wavelength",
            "20",
            "--cells-per-depth",
            "10",
            "--cfl",
            "0.5",
            "--simulation-periods",
            "7",
            "--solution-record-periods",
            "1",
            "--surface-elevation-record-periods",
            "5",
            "--frames-per-period",
            "20",
            "--outlet-relaxation-strength",
            "10.02",
            "--input-file",
            "regular.yaml",
            "--output",
            "regular.h5",
        )

        generated = setup.generate_regular_inputs(args)

        self.assertAlmostEqual(generated.wave.wavelength, 3.1173251893, places=9)
        self.assertEqual(generated.grid.tank_wavelengths, 7)
        self.assertAlmostEqual(
            generated.grid.tank_length, 7 * generated.wave.wavelength
        )
        self.assertAlmostEqual(generated.grid.xmin, -2 * generated.wave.wavelength)
        self.assertAlmostEqual(generated.grid.xmax, 10 * generated.wave.wavelength)
        self.assertEqual(generated.grid.nx, 240)
        self.assertEqual(generated.grid.nz, 10)
        self.assertEqual(generated.schedule.steps_per_period, 100)
        self.assertEqual(generated.schedule.num_steps, 700)
        self.assertEqual(generated.schedule.output_stride, 5)
        self.assertEqual(generated.schedule.solution_output_count, 21)
        self.assertEqual(generated.schedule.surface_elevation_output_count, 101)
        self.assertAlmostEqual(generated.schedule.dt, 0.015)
        self.assertAlmostEqual(generated.schedule.interval, 0.075)
        self.assertAlmostEqual(generated.schedule.solution_record_start_time, 9.0)
        self.assertAlmostEqual(
            generated.schedule.surface_elevation_record_start_time, 3.0
        )
        self.assertAlmostEqual(generated.schedule.end_time, 10.5)
        self.assertLessEqual(generated.schedule.dt, generated.dt_limit)

    def test_irregular_case_rounds_tank_and_aligns_dt_with_frames(self) -> None:
        args = self.parse(
            "irregular",
            "--water-depth",
            "10.0",
            "--peak-period",
            "15.24",
            "--significant-wave-height",
            "1.0",
            "--tank-length",
            "700",
            "--cells-per-wavelength",
            "20",
            "--cells-per-depth",
            "10",
            "--cfl",
            "0.5",
            "--simulation-periods",
            "240",
            "--solution-record-periods",
            "40",
            "--surface-elevation-record-periods",
            "200",
            "--frames-per-period",
            "10",
            "--component-count",
            "100",
            "--random-seed",
            "20260805",
            "--input-file",
            "irregular.yaml",
            "--output",
            "irregular.h5",
        )

        generated = setup.generate_irregular_inputs(args)

        self.assertAlmostEqual(generated.wave.wavelength, 146.5489406424, places=9)
        self.assertEqual(generated.grid.tank_wavelengths, 5)
        self.assertAlmostEqual(
            generated.grid.tank_length, 5 * generated.wave.wavelength
        )
        self.assertAlmostEqual(generated.grid.xmin, -2 * generated.wave.wavelength)
        self.assertAlmostEqual(generated.grid.xmax, 8 * generated.wave.wavelength)
        self.assertEqual(generated.grid.nx, 200)
        self.assertEqual(generated.grid.nz, 10)
        self.assertEqual(generated.schedule.steps_per_period, 300)
        self.assertEqual(generated.schedule.num_steps, 72000)
        self.assertEqual(generated.schedule.output_stride, 30)
        self.assertEqual(generated.schedule.solution_output_count, 401)
        self.assertEqual(generated.schedule.surface_elevation_output_count, 2001)
        self.assertAlmostEqual(generated.schedule.dt, 15.24 / 300)
        self.assertAlmostEqual(generated.schedule.interval, 15.24 / 10)
        self.assertAlmostEqual(
            generated.schedule.solution_record_start_time, 200 * 15.24
        )
        self.assertAlmostEqual(
            generated.schedule.surface_elevation_record_start_time, 40 * 15.24
        )
        self.assertAlmostEqual(generated.schedule.end_time, 240 * 15.24)
        self.assertLessEqual(generated.schedule.dt, generated.dt_limit)

    def test_domain_length_is_rounded_up_to_a_whole_wavelength(self) -> None:
        wave = setup.linear_wave_properties(0.7, 1.5)

        grid = setup.compute_grid(
            water_depth=0.7,
            requested_tank_length=3.2 * wave.wavelength,
            wave=wave,
            cells_per_wavelength=20,
            cells_per_depth=10,
            inlet_wavelengths=1.0,
            outlet_wavelengths=2.0,
            inlet_buffer_wavelengths=1.0,
            outlet_buffer_wavelengths=1.0,
        )

        self.assertEqual(grid.tank_wavelengths, 4)
        self.assertAlmostEqual(grid.tank_length, 4 * wave.wavelength)
        self.assertAlmostEqual(grid.inlet_buffer_length, wave.wavelength)
        self.assertAlmostEqual(grid.outlet_buffer_length, wave.wavelength)
        self.assertEqual(grid.nx, 180)

    def test_buffer_lengths_separate_forcing_zones_from_domain_of_interest(
        self,
    ) -> None:
        wave = setup.linear_wave_properties(0.7, 1.5)

        grid = setup.compute_grid(
            water_depth=0.7,
            requested_tank_length=3.2 * wave.wavelength,
            wave=wave,
            cells_per_wavelength=20,
            cells_per_depth=10,
            inlet_wavelengths=1.0,
            outlet_wavelengths=2.0,
            inlet_buffer_wavelengths=0.5,
            outlet_buffer_wavelengths=1.5,
        )

        inlet_zone_end = grid.xmin + wave.wavelength
        outlet_zone_start = grid.xmax - 2.0 * wave.wavelength
        self.assertAlmostEqual(inlet_zone_end, -0.5 * wave.wavelength)
        self.assertAlmostEqual(
            outlet_zone_start, grid.tank_length + 1.5 * wave.wavelength
        )
        self.assertAlmostEqual(grid.inlet_buffer_length, 0.5 * wave.wavelength)
        self.assertAlmostEqual(grid.outlet_buffer_length, 1.5 * wave.wavelength)

    def test_incompatible_record_rate_is_rejected(self) -> None:
        with self.assertRaisesRegex(ValueError, "must be an integer"):
            setup.align_flow_schedule(
                period=2.0,
                simulation_periods=setup.Fraction("4"),
                solution_record_periods=setup.Fraction("1/3"),
                surface_elevation_record_periods=setup.Fraction("1"),
                frames_per_period=10,
                dt_limit=0.1,
            )

    def test_record_window_cannot_exceed_simulation(self) -> None:
        with self.assertRaisesRegex(ValueError, "must not exceed"):
            setup.align_flow_schedule(
                period=2.0,
                simulation_periods=setup.Fraction("4"),
                solution_record_periods=setup.Fraction("5"),
                surface_elevation_record_periods=setup.Fraction("4"),
                frames_per_period=10,
                dt_limit=0.1,
            )

        with self.assertRaisesRegex(ValueError, "surface-elevation-record-periods"):
            setup.align_flow_schedule(
                period=2.0,
                simulation_periods=setup.Fraction("4"),
                solution_record_periods=setup.Fraction("1"),
                surface_elevation_record_periods=setup.Fraction("5"),
                frames_per_period=10,
                dt_limit=0.1,
            )

    def test_irregular_frequency_ratios_must_be_ordered(self) -> None:
        args = self.parse(
            "irregular",
            "--water-depth",
            "10",
            "--peak-period",
            "15.24",
            "--significant-wave-height",
            "1",
            "--tank-length",
            "700",
            "--simulation-periods",
            "2",
            "--minimum-frequency-ratio",
            "2",
            "--maximum-frequency-ratio",
            "1",
            "--input-file",
            "irregular.yaml",
            "--output",
            "irregular.h5",
        )

        with self.assertRaisesRegex(ValueError, "must be greater"):
            setup.generate_irregular_inputs(args)

    def test_single_final_frame_is_supported(self) -> None:
        schedule = setup.align_flow_schedule(
            period=2.0,
            simulation_periods=setup.Fraction("5"),
            solution_record_periods=setup.Fraction("0"),
            surface_elevation_record_periods=setup.Fraction("0"),
            frames_per_period=10,
            dt_limit=0.1,
        )

        self.assertEqual(schedule.solution_output_count, 1)
        self.assertEqual(schedule.surface_elevation_output_count, 1)
        self.assertEqual(schedule.solution_record_start_time, 10.0)
        self.assertEqual(schedule.surface_elevation_record_start_time, 10.0)
        self.assertEqual(schedule.end_time, 10.0)

    def test_rendered_file_is_petsc_yaml(self) -> None:
        args = self.parse(
            "regular",
            "--water-depth",
            "1",
            "--wave-period",
            "2",
            "--wave-height",
            "0.1",
            "--tank-length",
            "10",
            "--simulation-periods",
            "10",
            "--input-file",
            "regular.yaml",
            "--output",
            "regular.h5",
        )

        rendered = setup.render_yaml(setup.generate_regular_inputs(args))

        self.assertIn("sim:\n  start_time: 0.00000000e+00\n  end_time:", rendered)
        self.assertIn("  dt:", rendered)
        for option in (
            "flow_field_output_start_time:",
            "flow_surface_elevation_output_start_time:",
            "flow_output_interval:",
            "wave_height:",
            "wave_period:",
            "water_depth:",
            "xmin:",
            "xmax:",
            "nx:",
            "nz:",
            "# inlet_buffer_wavelengths: 1.0",
            "# outlet_buffer_wavelengths: 1.0",
            "gamma:",
            "ramp_up_time:",
            'output: "regular.h5"',
        ):
            self.assertIn(option, rendered)
        self.assertNotIn("-sim_dt", rendered)

    def test_floats_use_scientific_notation_with_eight_decimal_places(self) -> None:
        self.assertEqual(setup.format_value(1.0), "1.00000000e+00")
        self.assertEqual(setup.format_value(setup.Fraction("1/2")), "5.00000000e-01")


if __name__ == "__main__":
    unittest.main()
