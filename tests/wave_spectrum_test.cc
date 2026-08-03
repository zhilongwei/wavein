#include "wavein/jonswap.h"
#include "wavein/pierson_moskowitz.h"
#include "wavein/wave_spectrum.h"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <petscsys.h>

namespace
{

struct AngularFrequencyGrid
{
        PetscReal minimum;
        PetscReal maximum;
        PetscInt point_count;
};

[[nodiscard]] PetscReal zeroth_moment(const wavein::WaveSpectrum &spectrum,
                                      const AngularFrequencyGrid &grid)
{
    const PetscReal spacing =
        (grid.maximum - grid.minimum) / static_cast<PetscReal>(grid.point_count - 1);
    PetscReal integral = 0.5 * (spectrum.spectrum(grid.minimum) + spectrum.spectrum(grid.maximum));

    for (PetscInt index = 1; index < grid.point_count - 1; ++index)
    {
        const PetscReal omega = grid.minimum + static_cast<PetscReal>(index) * spacing;
        integral += spectrum.spectrum(omega);
    }

    return integral * spacing;
}

constexpr PetscReal Hs = 2.0;
constexpr PetscReal Tp = 8.0;
constexpr PetscReal peak_enhancement_factor = 3.3;
constexpr AngularFrequencyGrid frequency_grid{0.2 / Tp, 200.0 / Tp, 20000};

} // namespace

TEST_CASE("Pierson-Moskowitz spectrum recovers the target significant wave height", "[spectrum]")
{
    const wavein::PiersonMoskowitz spectrum(Hs, Tp);
    const PetscReal m0 = zeroth_moment(spectrum, frequency_grid);

    REQUIRE(4.0 * PetscSqrtReal(m0) == Catch::Approx(Hs).epsilon(1.0e-6));
}

TEST_CASE("JONSWAP spectrum recovers the target significant wave height", "[spectrum]")
{
    const wavein::Jonswap spectrum(Hs, Tp, peak_enhancement_factor);
    const PetscReal m0 = zeroth_moment(spectrum, frequency_grid);

    REQUIRE(4.0 * PetscSqrtReal(m0) == Catch::Approx(Hs).epsilon(1.5e-3));
}
