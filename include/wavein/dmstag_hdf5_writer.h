#pragma once

#include <petscdmstag.h>

#include "wavein/airy_wave.h"
#include "wavein/irregular_waves.h"

namespace wavein
{

class DMStagHDF5Writer
{
    public:
        DMStagHDF5Writer() = delete;

        DMStagHDF5Writer(MPI_Comm comm, DM dm, PetscViewer viewer);

        DMStagHDF5Writer(const DMStagHDF5Writer &) = delete;
        DMStagHDF5Writer &operator=(const DMStagHDF5Writer &) = delete;
        ~DMStagHDF5Writer()
        {
            PetscCallAbort(comm_, destroy());
        }

        PetscErrorCode push_group();
        PetscErrorCode pop_group();

        // Writes the metadata for the AiryWave object. This includes the wave height, water depth,
        // wave period, and wavelength.
        PetscErrorCode write_airy_wave(const AiryWave &wave);

        PetscErrorCode write_irregular_waves(const IrregularWaves &waves);

        PetscErrorCode write_domain(PetscReal xmin, PetscReal xmax, PetscReal physical_xmin,
                                    PetscReal physical_xmax);

        PetscErrorCode write_simulation(PetscReal start_time, PetscReal end_time, PetscReal dt,
                                        PetscInt nsteps);

        // Writes the packed DMSTAG vector. Construction writes entry-aligned x, z, location,
        // i, and j datasets for decomposition-independent post-processing.
        PetscErrorCode write_solution(Vec sol, PetscReal t, const char *info = nullptr);

        // Writes the top-row DMSTAG_UP values as a 1D vector with one value per x element.
        PetscErrorCode write_surface_elevation(Vec eta, PetscReal t, const char *info = nullptr);

    private:
        struct EtaField
        {
                Vec values = nullptr;
                Vec time = nullptr;
                VecScatter scatter = nullptr;
                PetscInt timestep = 0;
        };

        PetscErrorCode create_eta_field(DM dm, Vec prototype);
        PetscErrorCode write_metadata(DM dm, Vec prototype);
        PetscErrorCode destroy();

        MPI_Comm comm_ = MPI_COMM_NULL;
        PetscViewer viewer_ = nullptr;
        Vec sol_time_ = nullptr;
        PetscInt sol_timestep_ = 0;
        EtaField eta_;
};

} // namespace wavein
