#include "wavein/jonswap.h"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <numeric>
#include <random>
#include <utility>
#include <vector>

namespace
{

constexpr PetscReal Hs = 2.00; // target significant wave height
constexpr PetscReal Tp = 10.0; // target peak period
constexpr PetscReal peak_enhancement_factor = 3.3;
constexpr PetscReal Tz = 0.834 * Tp; // zero-crossing period for JONSWAP spectrum with gamma = 3.3
constexpr PetscReal omega_min = 1.00 / Tz;        // minimum angular frequency for JONSWAP spectrum
constexpr PetscReal omega_max = 20.0 / Tz;        // maximum angular frequency for JONSWAP spectrum
constexpr PetscReal duration = 3.0 * 60.0 * 60.0; // 3 hours
constexpr PetscInt omega_point_count = 300;       // number of angular frequencies
constexpr PetscInt samples_per_zero_crossing_period = 40;
constexpr std::mt19937::result_type random_seed = 20260805U;
constexpr PetscReal mean_tolerance = 1.0e-3;
constexpr PetscReal variance_tolerance = 5.0e-3;

struct WaveComponent
{
        PetscReal omega;
        PetscReal amplitude;
        PetscReal phase;
};

struct SpectralRealization
{
        std::vector<PetscReal> elevation;
        PetscReal spectral_variance;
};

[[nodiscard]] SpectralRealization realize_surface_elevation(const wavein::WaveSpectrum &spectrum)
{
    const PetscReal delta_omega =
        (omega_max - omega_min) / static_cast<PetscReal>(omega_point_count - 1);
    std::mt19937 generator(random_seed);
    std::uniform_real_distribution<PetscReal> phase_distribution(0.0, 2.0 * PETSC_PI);
    std::vector<WaveComponent> components;
    components.reserve(static_cast<std::size_t>(omega_point_count - 1));
    PetscReal spectral_variance = 0.0;

    for (PetscInt strip = 0; strip < omega_point_count - 1; ++strip)
    {
        const PetscReal omega = omega_min + (static_cast<PetscReal>(strip) + 0.5) * delta_omega;
        const PetscReal component_variance = spectrum.spectrum(omega) * delta_omega;
        const PetscReal amplitude = PetscSqrtReal(2.0 * component_variance);
        spectral_variance += component_variance;
        components.push_back({omega, amplitude, phase_distribution(generator)});
    }

    const auto zero_crossing_count = static_cast<PetscInt>(duration / Tz);
    const PetscInt sample_count = zero_crossing_count * samples_per_zero_crossing_period;
    const PetscReal delta_time = duration / static_cast<PetscReal>(sample_count - 1);
    std::vector<PetscReal> elevation(static_cast<std::size_t>(sample_count), 0.0);

    for (PetscInt sample = 0; sample < sample_count; ++sample)
    {
        const PetscReal time = static_cast<PetscReal>(sample) * delta_time;
        PetscReal eta = 0.0;
        for (const WaveComponent &component : components)
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
    const SpectralRealization realization = realize_surface_elevation(spectrum);
    const PetscReal elevation_mean = mean(realization.elevation);
    const PetscReal elevation_variance = variance(realization.elevation, elevation_mean);

    REQUIRE(realization.spectral_variance > 0.0);
    REQUIRE(elevation_mean == Catch::Approx(0.0).margin(mean_tolerance));
    REQUIRE(elevation_variance ==
            Catch::Approx(realization.spectral_variance).epsilon(variance_tolerance));
}
