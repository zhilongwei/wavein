#include "wavein/petsc_info.h"

#include <petscsys.h>

int main(int argc, char **argv)
{
    const auto initialize_error = PetscInitialize(&argc, &argv, nullptr, nullptr);
    if (initialize_error != PETSC_SUCCESS)
    {
        return static_cast<int>(initialize_error);
    }

    PetscPrintf(PETSC_COMM_WORLD, "%s\n", wavein::petsc_version().c_str());

    return static_cast<int>(PetscFinalize());
}
