#include "wavein/jonswap.h"
#include "wavein/pierson_moskowitz.h"

[[nodiscard]] PetscReal wavein::Jonswap::spectrum(PetscReal omega) const
{
    const PetscReal omegap = 2.0 * PETSC_PI / Tp_;
    const PetscReal sigma = omega > omegap ? sigma_b_ : sigma_a_;

    const PetscReal A = 1.0 - 0.287 * PetscLogReal(gamma_);
    const PetscReal beta = (omega - omegap) / (sigma * omegap);
    const PiersonMoskowitz pm(Hs_, Tp_);

    return A * pm.spectrum(omega) * PetscPowReal(gamma_, PetscExpReal(-0.5 * beta * beta));
}