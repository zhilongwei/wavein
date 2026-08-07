#include "wavein/airy_wave.h"
#include "wavein/irregular_waves.h"
#include "wavein/jonswap.h"
#include "wavein/projection.h"
#include "wavein/wave_tank.h"
#include "wavein/wavemaker.h"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <petscdmstag.h>

#include <numeric>
#include <utility>
#include <vector>

namespace
{

constexpr PetscReal Hs = 2.00; // target significant wave height
constexpr PetscReal Tp = 10.0; // target peak period
constexpr PetscReal peak_enhancement_factor = 3.3;
constexpr PetscReal water_depth = 20.0;
constexpr PetscReal Tz = 0.834 * Tp; // zero-crossing period for JONSWAP spectrum with gamma = 3.3
constexpr PetscReal omega_min = 1.00 / Tz;        // minimum angular frequency for JONSWAP spectrum
constexpr PetscReal omega_max = 20.0 / Tz;        // maximum angular frequency for JONSWAP spectrum
constexpr PetscReal duration = 3.0 * 60.0 * 60.0; // 3 hours
constexpr PetscInt omega_point_count = 300;       // number of angular frequencies
constexpr PetscInt samples_per_zero_crossing_period = 40;
constexpr unsigned long random_seed = 20260805UL;
constexpr PetscReal mean_tolerance = 1.0e-3;
constexpr PetscReal variance_tolerance = 5.0e-3;

struct SpectralRealization
{
        std::vector<PetscReal> elevation;
        PetscReal spectral_variance;
};

[[nodiscard]] SpectralRealization realize_surface_elevation(const wavein::IrregularWaves &waves)
{
    PetscReal spectral_variance = 0.0;
    for (const wavein::WaveComponent &component : waves.components())
    {
        spectral_variance += 0.5 * PetscSqr(component.amplitude);
    }

    const auto zero_crossing_count = static_cast<PetscInt>(duration / Tz);
    const PetscInt sample_count = zero_crossing_count * samples_per_zero_crossing_period;
    const PetscReal delta_time = duration / static_cast<PetscReal>(sample_count - 1);
    std::vector<PetscReal> elevation(static_cast<std::size_t>(sample_count), 0.0);

    for (PetscInt sample = 0; sample < sample_count; ++sample)
    {
        const PetscReal time = static_cast<PetscReal>(sample) * delta_time;
        PetscReal eta = 0.0;
        for (const wavein::WaveComponent &component : waves.components())
        {
            eta += component.amplitude * PetscCosReal(component.omega * time + component.phase);
        }
        elevation[static_cast<std::size_t>(sample)] = eta;
    }

    return {std::move(elevation), spectral_variance};
}

[[nodiscard]] PetscReal mean(const std::vector<PetscReal> &values)
{
    return std::accumulate(values.begin(), values.end(), PetscReal{0.0}) /
           static_cast<PetscReal>(values.size());
}

[[nodiscard]] PetscReal variance(const std::vector<PetscReal> &values, PetscReal values_mean)
{
    PetscReal squared_deviation_sum = 0.0;
    for (const PetscReal value : values)
    {
        squared_deviation_sum += PetscSqr(value - values_mean);
    }
    return squared_deviation_sum / static_cast<PetscReal>(values.size());
}

} // namespace

TEST_CASE("Irregular-wave surface elevation recovers the spectral variance", "[irregular_waves]")
{
    const wavein::Jonswap spectrum(Hs, Tp, peak_enhancement_factor);
    const wavein::IrregularWaves waves(PETSC_COMM_WORLD, spectrum, water_depth, omega_min,
                                       omega_max, omega_point_count - 1, random_seed);
    const SpectralRealization realization = realize_surface_elevation(waves);
    const PetscReal elevation_mean = mean(realization.elevation);
    const PetscReal elevation_variance = variance(realization.elevation, elevation_mean);

    REQUIRE(realization.spectral_variance > 0.0);
    REQUIRE(elevation_mean == Catch::Approx(0.0).margin(mean_tolerance));
    REQUIRE(elevation_variance ==
            Catch::Approx(realization.spectral_variance).epsilon(variance_tolerance));
}

TEST_CASE("Irregular-wave tank advances a multi-component field", "[irregular_waves][wave_tank]")
{
    constexpr PetscInt nx = 16;
    constexpr PetscInt nz = 4;
    constexpr PetscInt component_count = 8;
    constexpr PetscReal dt = 0.01;
    constexpr PetscReal time = 0.37;

    const MPI_Comm comm = PETSC_COMM_WORLD;
    const wavein::Jonswap spectrum(Hs, Tp, peak_enhancement_factor);
    const wavein::IrregularWaves waves(comm, spectrum, water_depth, omega_min, omega_max,
                                       component_count, random_seed);
    const wavein::AiryWave peak_wave(comm, water_depth, Tp, Hs);

    const PetscReal wavelength = peak_wave.wavelength();
    const PetscReal xmin = -wavelength;
    const PetscReal xmax = 3.0 * wavelength;
    const PetscReal zmin = -water_depth;
    const PetscReal zmax = 0.0;

    DM dm = nullptr;
    PetscCallAbort(comm, DMStagCreate2d(comm, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, nx, nz,
                                        PETSC_DECIDE, PETSC_DECIDE, 0, 1, 1, DMSTAG_STENCIL_BOX, 1,
                                        nullptr, nullptr, &dm));
    PetscCallAbort(comm, DMSetUp(dm));
    PetscCallAbort(comm, DMStagSetUniformCoordinatesProduct(dm, xmin, xmax, zmin, zmax, 0.0, 0.0));

    const PetscReal dx = (xmax - xmin) / static_cast<PetscReal>(nx);
    const PetscReal dz = (zmax - zmin) / static_cast<PetscReal>(nz);

    Vec sol = nullptr;
    Vec eta = nullptr;
    PetscCallAbort(comm, DMCreateGlobalVector(dm, &sol));
    PetscCallAbort(comm, DMCreateGlobalVector(dm, &eta));
    PetscCallAbort(comm, VecZeroEntries(sol));
    PetscCallAbort(comm, VecZeroEntries(eta));

    PetscReal sol_norm = 0.0;
    PetscReal eta_norm = 0.0;
    PetscReal divergence_norm = 0.0;
    {
        wavein::Wavemaker wavemaker(comm, dm, wavelength, xmin, xmax, 1.0, 1.0, 1.0);
        wavein::Projection projection(comm, dm, dx, dz);
        wavein::WaveTank wave_tank(comm, dm, waves, wavemaker, projection);

        PetscCallAbort(comm, wave_tank.update(sol, eta, time, dt, 1.0));
        PetscCallAbort(comm, VecNorm(sol, NORM_2, &sol_norm));
        PetscCallAbort(comm, VecNorm(eta, NORM_2, &eta_norm));
        PetscCallAbort(comm, projection.divergence_norm(sol, &divergence_norm));
    }

    PetscCallAbort(comm, VecDestroy(&eta));
    PetscCallAbort(comm, VecDestroy(&sol));
    PetscCallAbort(comm, DMDestroy(&dm));

    REQUIRE_FALSE(PetscIsInfOrNanReal(sol_norm));
    REQUIRE_FALSE(PetscIsInfOrNanReal(eta_norm));
    REQUIRE_FALSE(PetscIsInfOrNanReal(divergence_norm));
    REQUIRE(sol_norm > 0.0);
    REQUIRE(eta_norm > 0.0);
    REQUIRE(divergence_norm <= 1.0e-10 * sol_norm);
}
