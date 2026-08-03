#include "wavein/petsc_info.h"

#include <catch2/catch_session.hpp>
#include <catch2/catch_test_macros.hpp>
#include <petscsys.h>

int main(int argc, char *argv[])
{
    // Keep Catch2's discovery flags out of PETSc's options database.
    const auto initialize_error = PetscInitializeNoArguments();
    if (initialize_error != PETSC_SUCCESS)
    {
        return static_cast<int>(initialize_error);
    }

    const int test_result = Catch::Session().run(argc, argv);
    const auto finalize_error = PetscFinalize();
    return test_result != 0 ? test_result : static_cast<int>(finalize_error);
}

TEST_CASE("PETSc reports an installed version", "[petsc]")
{
    PetscInt major = 0;
    PetscInt minor = 0;
    PetscInt subminor = 0;
    PetscInt release = 0;

    REQUIRE(PetscGetVersionNumber(&major, &minor, &subminor, &release) == PETSC_SUCCESS);
    REQUIRE(major >= 3);
    REQUIRE_FALSE(wavein::petsc_version().empty());
}
