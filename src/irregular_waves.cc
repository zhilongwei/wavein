#include "wavein/irregular_waves.h"
#include "wavein/airy_wave.h"

PetscErrorCode wavein::IrregularWaves::velocity_spectra(Vec omegas, PetscReal z,
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