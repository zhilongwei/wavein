#pragma once

#include <petscdmstag.h>
#include <petscvec.h>

namespace wavein
{

class Wavemaker
{
    public:
        Wavemaker() = delete;

        Wavemaker(MPI_Comm comm, DM dm, PetscReal wavelength, PetscReal xmin, PetscReal xmax,
                  PetscReal nin, PetscReal nout, PetscReal gamma);

        Wavemaker(const Wavemaker &) = delete;
        Wavemaker &operator=(const Wavemaker &) = delete;
        ~Wavemaker() noexcept;

        PetscErrorCode force(Vec ref_vel, Vec ref_eta, PetscReal dt, PetscReal factor, Vec vel,
                             Vec eta) noexcept;

    private:
        [[nodiscard]] PetscReal exp_blender_func(PetscReal x) const
        {
            return (PetscExpReal(x * x) - 1.0) / (PetscExpReal(1.0) - 1.0);
        }
        MPI_Comm comm_ = MPI_COMM_NULL;
        Vec blender_ = nullptr;
        Vec delta_vel_ = nullptr;
        Vec delta_eta_ = nullptr;
};

} // namespace wavein