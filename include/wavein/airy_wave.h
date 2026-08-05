#pragma once

#include <petscsys.h>

namespace wavein
{

class AiryWave
{
    public:
        AiryWave() = delete;
        AiryWave(MPI_Comm comm, PetscReal h, PetscReal T, PetscReal H = 2.0,
                 PetscReal atol = 1.0e-8);
        AiryWave(const AiryWave &) = delete;
        AiryWave &operator=(const AiryWave &) = delete;

        [[nodiscard]] PetscComplex horizontal_velocity_transfer(PetscReal z) const noexcept;
        [[nodiscard]] PetscComplex vertical_velocity_transfer(PetscReal z) const noexcept;

        [[nodiscard]] PetscReal surface_elevation(PetscReal x, PetscReal t,
                                                  PetscReal phase = 0.0) const noexcept
        {
            return 0.5 * H_ *
                   PetscRealPartComplex(PetscExpComplex(PETSC_i * (omega_ * t - k_ * x + phase)));
        }

        [[nodiscard]] PetscReal horizontal_velocity(PetscReal x, PetscReal z, PetscReal t,
                                                    PetscReal phase = 0.0) const noexcept
        {
            return PetscRealPartComplex(0.5 * H_ * horizontal_velocity_transfer(z) *
                                        PetscExpComplex(PETSC_i * (omega_ * t - k_ * x + phase)));
        }

        [[nodiscard]] PetscReal vertical_velocity(PetscReal x, PetscReal z, PetscReal t,
                                                  PetscReal phase = 0.0) const noexcept
        {
            return PetscRealPartComplex(0.5 * H_ * vertical_velocity_transfer(z) *
                                        PetscExpComplex(PETSC_i * (omega_ * t - k_ * x + phase)));
        }

        [[nodiscard]] PetscReal wave_height() const
        {
            return H_;
        }

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
        const PetscReal H_;     // wave height (optional)
        const PetscReal atol_;  // tolerance for solving the dispersion relation
        const PetscReal omega_; // wave frequency
        const PetscReal kh_;    // wavenumber * water depth
        const PetscReal k_;     // wavenumber
        const PetscReal L_;     // wavelength
};

} // namespace wavein