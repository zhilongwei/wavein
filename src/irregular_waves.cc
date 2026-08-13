#include "wavein/irregular_waves.h"
#include "wavein/airy_wave.h"

#include <random>

namespace wavein
{

IrregularWaves::IrregularWaves(MPI_Comm comm, const WaveSpectrum &spectrum, PetscReal h)
    : comm_(comm), spectrum_(spectrum), h_(h),
      wavelength_at_peak_period_(AiryWave(comm_, h_, spectrum_.peak_period()).wavelength())
{
}

// NOLINTBEGIN(bugprone-easily-swappable-parameters)
IrregularWaves::IrregularWaves(MPI_Comm comm, const WaveSpectrum &spectrum, PetscReal h,
                               PetscReal omega_min, PetscReal omega_max, PetscInt component_count,
                               unsigned long random_seed)
    : IrregularWaves(comm, spectrum, h)
{
    PetscFunctionBeginUser;

    const PetscReal delta_omega = (omega_max - omega_min) / static_cast<PetscReal>(component_count);
    std::mt19937 generator(static_cast<std::mt19937::result_type>(random_seed));
    std::uniform_real_distribution<PetscReal> phase_distribution(0.0, 2.0 * PETSC_PI);

    components_.reserve(static_cast<std::size_t>(component_count));
    for (PetscInt component = 0; component != component_count; ++component)
    {
        const PetscReal omega = omega_min + (static_cast<PetscReal>(component) + 0.5) * delta_omega;
        const PetscReal amplitude = PetscSqrtReal(2.0 * spectrum_.spectrum(omega) * delta_omega);
        components_.push_back({omega, amplitude, phase_distribution(generator)});
    }

    PetscFunctionReturnVoid();
}
// NOLINTEND(bugprone-easily-swappable-parameters)

PetscErrorCode IrregularWaves::velocity_spectra(Vec omegas, PetscReal z,
                                                Vec horizontal_velocity_spectrum,
                                                Vec vertical_velocity_spectrum) const
{
    PetscFunctionBeginUser;

    PetscInt n;
    PetscCall(VecGetLocalSize(omegas, &n));

    const PetscReal *c_arr_omegas;
    PetscReal *c_arr_horizontal_velocity, *c_arr_vertical_velocity;
    PetscCall(VecGetArrayRead(omegas, &c_arr_omegas));
    PetscCall(VecGetArray(horizontal_velocity_spectrum, &c_arr_horizontal_velocity));
    PetscCall(VecGetArray(vertical_velocity_spectrum, &c_arr_vertical_velocity));

    for (PetscInt i = 0; i != n; ++i)
    {
        PetscReal omega = c_arr_omegas[i];
        const AiryWave wave(comm_, h_, 2.0 * PETSC_PI / omega);
        const PetscReal spectral_density = spectrum_.spectrum(omega);
        c_arr_horizontal_velocity[i] =
            PetscSqr(PetscAbsComplex(wave.horizontal_velocity_transfer(z))) * spectral_density;
        c_arr_vertical_velocity[i] =
            PetscSqr(PetscAbsComplex(wave.vertical_velocity_transfer(z))) * spectral_density;
    }

    PetscCall(VecRestoreArrayRead(omegas, &c_arr_omegas));
    PetscCall(VecRestoreArray(horizontal_velocity_spectrum, &c_arr_horizontal_velocity));
    PetscCall(VecRestoreArray(vertical_velocity_spectrum, &c_arr_vertical_velocity));

    PetscFunctionReturn(PETSC_SUCCESS);
}

} // namespace wavein
