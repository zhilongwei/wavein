#include "wavein/airy_wave.h"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <petscsys.h>

namespace
{

constexpr PetscReal g = 9.80665;
constexpr PetscReal expected_h = 0.70;
constexpr PetscReal expected_L = 2.00;
constexpr PetscReal expected_k = 2.0 * PETSC_PI / expected_L;
constexpr PetscReal solver_tolerance = 1.0e-10;
constexpr PetscReal comparison_tolerance = 1.0e-8;
constexpr PetscReal velocity_test_h = 10.0;
constexpr PetscReal velocity_test_T = 6.0;
constexpr PetscReal dz = 1.0e-5;
constexpr PetscReal differential_equation_tolerance = 1.0e-8;
constexpr PetscReal first_relative_depth = -0.8;
constexpr PetscReal last_relative_depth = -0.2;
constexpr PetscInt depth_sample_count = 10;

const PetscComplex imaginary_unit{0.0, 1.0};

const PetscReal expected_omega =
    PetscSqrtReal(g * expected_k * PetscTanhReal(expected_k * expected_h));
const PetscReal expected_wave_period = 2.0 * PETSC_PI / expected_omega;
const PetscReal expected_celerity = expected_omega / expected_k;

} // namespace

TEST_CASE("Linear wave properties satisfy the dispersion relation", "[airy]")
{
    const wavein::AiryWave wave(PETSC_COMM_WORLD, expected_h, expected_wave_period, 2.0,
                                solver_tolerance);

    REQUIRE(wave.water_depth() == Catch::Approx(expected_h).epsilon(comparison_tolerance));
    REQUIRE(wave.wave_period() ==
            Catch::Approx(expected_wave_period).epsilon(comparison_tolerance));
    REQUIRE(wave.wave_frequency() == Catch::Approx(expected_omega).epsilon(comparison_tolerance));
    REQUIRE(wave.wavenumber() == Catch::Approx(expected_k).epsilon(comparison_tolerance));
    REQUIRE(wave.wavelength() == Catch::Approx(expected_L).epsilon(comparison_tolerance));
    REQUIRE(wave.celerity() == Catch::Approx(expected_celerity).epsilon(comparison_tolerance));

    const PetscReal dispersion_residual =
        wave.wave_frequency() * wave.wave_frequency() -
        g * wave.wavenumber() * PetscTanhReal(wave.wavenumber() * wave.water_depth());

    REQUIRE(PetscAbsReal(dispersion_residual) <= 2.0 * g / expected_h * solver_tolerance);
}

TEST_CASE("Linear wave velocities are incompressible and irrotational", "[airy]")
{
    const wavein::AiryWave wave(PETSC_COMM_WORLD, velocity_test_h, velocity_test_T, 2.0,
                                solver_tolerance);

    for (PetscInt index = 0; index < depth_sample_count; ++index)
    {
        const PetscReal interpolation_factor =
            static_cast<PetscReal>(index) / static_cast<PetscReal>(depth_sample_count - 1);
        const PetscReal relative_depth =
            first_relative_depth +
            interpolation_factor * (last_relative_depth - first_relative_depth);
        const PetscReal z = relative_depth * wave.water_depth();

        const PetscComplex du_dz = (wave.horizontal_velocity_transfer(z + dz) -
                                    wave.horizontal_velocity_transfer(z - dz)) /
                                   (2.0 * dz);
        const PetscComplex dw_dz =
            (wave.vertical_velocity_transfer(z + dz) - wave.vertical_velocity_transfer(z - dz)) /
            (2.0 * dz);

        const PetscComplex incompressibility_residual =
            -imaginary_unit * wave.wavenumber() * wave.horizontal_velocity_transfer(z) + dw_dz;
        const PetscComplex irrotationality_residual =
            -imaginary_unit * wave.wavenumber() * wave.vertical_velocity_transfer(z) - du_dz;

        REQUIRE(PetscAbsComplex(incompressibility_residual) ==
                Catch::Approx(0.0).margin(differential_equation_tolerance));
        REQUIRE(PetscAbsComplex(irrotationality_residual) ==
                Catch::Approx(0.0).margin(differential_equation_tolerance));
    }
}
