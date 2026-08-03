#pragma once

#include <petscsys.h>

namespace wavein
{

class AiryWave
{
    public:
        AiryWave() = delete;
        AiryWave(MPI_Comm comm, PetscReal h, PetscReal T, PetscReal atol = 1.0e-8);
        AiryWave(const AiryWave &) = delete;
        AiryWave &operator=(const AiryWave &) = delete;

        [[nodiscard]] PetscComplex horizontal_velocity_transfer(PetscReal z) const noexcept;
        [[nodiscard]] PetscComplex vertical_velocity_transfer(PetscReal z) const noexcept;

        [[nodiscard]] PetscReal water_depth() const
        {
            return h_;
        }

        [[nodiscard]] PetscReal wave_period() const
        {
            return T_;
        }

        [[nodiscard]] PetscReal wave_frequency() const
        {
            return omega_;
        }

        [[nodiscard]] PetscReal wavenumber() const
        {
            return k_;
        }

        [[nodiscard]] PetscReal wavelength() const
        {
            return L_;
        }

        [[nodiscard]] PetscReal celerity() const
        {
            return omega_ / k_;
        }

    private:
        [[nodiscard]] PetscReal nondimensional_dispersion_relation(
            MPI_Comm comm, PetscReal nondimensional_frequency) const noexcept;
        const PetscReal h_;     // water depth
        const PetscReal T_;     // wave period
        const PetscReal atol_;  // tolerance for solving the dispersion relation
        const PetscReal omega_; // wave frequency
        const PetscReal kh_;    // wavenumber * water depth
        const PetscReal k_;     // wavenumber
        const PetscReal L_;     // wavelength
};

} // namespace wavein