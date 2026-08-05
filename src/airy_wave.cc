#include <petscsys.h>

#include "wavein/airy_wave.h"
#include "wavein/constants.h"

namespace wavein
{

// NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
AiryWave::AiryWave(MPI_Comm comm, PetscReal h, PetscReal T, PetscReal H, PetscReal atol)
    : h_(h), T_(T), H_(H), atol_(atol), omega_(2.0 * PETSC_PI / T),
      kh_(nondimensional_dispersion_relation(comm, omega_ * PetscSqrtReal(h_ / kGA))), k_(kh_ / h_),
      L_(2.0 * PETSC_PI / k_)
{
    // This is an intentional empty constructor body.
}

PetscComplex AiryWave::horizontal_velocity_transfer(PetscReal z) const noexcept
{
    return PetscComplex{kGA * k_ / omega_ * PetscCoshReal(k_ * (z + h_)) / PetscCoshReal(kh_), 0.0};
}

PetscComplex AiryWave::vertical_velocity_transfer(PetscReal z) const noexcept
{
    return PetscComplex{0.0, kGA * k_ / omega_ * PetscSinhReal(k_ * (z + h_)) / PetscCoshReal(kh_)};
}

// NOLINTBEGIN(bugprone-easily-swappable-parameters)
PetscReal AiryWave::nondimensional_dispersion_relation(
    MPI_Comm comm, PetscReal nondimensional_frequency) const noexcept
{
    constexpr PetscInt max_iterations = 50;

    PetscReal kh =
        PetscMax(nondimensional_frequency, nondimensional_frequency * nondimensional_frequency);

    for (PetscInt i = 0; i != max_iterations; ++i)
    {
        const PetscReal tanh_kh = PetscTanhReal(kh);
        const PetscReal residual =
            kh * tanh_kh - nondimensional_frequency * nondimensional_frequency;

        if (PetscAbsReal(residual) <= atol_)
        {
            return kh;
        }

        const PetscReal derivative = tanh_kh + kh * (1.0 - tanh_kh * tanh_kh);
        kh -= residual / derivative;
    }

    const PetscReal final_residual =
        kh * PetscTanhReal(kh) - nondimensional_frequency * nondimensional_frequency;

    PetscCheckAbort(PetscAbsReal(final_residual) <= atol_, comm, PETSC_ERR_CONV_FAILED,
                    "Linear wave dispersion solver did not converge within %" PetscInt_FMT
                    " iterations",
                    max_iterations);

    return kh;
}
// NOLINTEND(bugprone-easily-swappable-parameters)

} // namespace wavein