#pragma once

#include <petscdmstag.h>
#include <petscvec.h>

namespace wavein
{

struct ForcingZoneGeometry
{
        PetscReal xmin;
        PetscReal xmax;
        PetscReal inlet_length;
        PetscReal outlet_length;

        [[nodiscard]] PetscBool is_in_inlet_forcing_zone(PetscReal x) const noexcept
        {
            return x >= xmin && x <= xmin + inlet_length;
        }

        [[nodiscard]] PetscBool is_in_outlet_forcing_zone(PetscReal x) const noexcept
        {
            return x >= xmax - outlet_length && x <= xmax;
        }
};

class Wavemaker
{
    public:
        Wavemaker() = delete;

        Wavemaker(MPI_Comm comm, DM dm, PetscReal wavelength, PetscReal xmin, PetscReal xmax,
                  PetscReal nin, PetscReal nout, PetscReal gamma);

        Wavemaker(const Wavemaker &) = delete;
        Wavemaker &operator=(const Wavemaker &) = delete;
        ~Wavemaker() noexcept;

        [[nodiscard]] ForcingZoneGeometry forcing_zone_geometry() const noexcept
        {
            return geometry_;
        }

        PetscErrorCode force(Vec ref_vel, Vec ref_eta, PetscReal dt, PetscReal ramp, Vec vel,
                             Vec eta) noexcept;

    private:
        [[nodiscard]] PetscReal exp_blender_func(PetscReal x) const
        {
            return (PetscExpReal(x * x) - 1.0) / (PetscExpReal(1.0) - 1.0);
        }
        MPI_Comm comm_ = MPI_COMM_NULL;
        ForcingZoneGeometry geometry_;
        PetscReal gamma_ = 0.0; // Outlet source-relaxation rate in 1/s.
        Vec inlet_blender_ = nullptr;
        Vec outlet_blender_ = nullptr;
        Vec delta_vel_ = nullptr;
        Vec delta_eta_ = nullptr;
};

} // namespace wavein
