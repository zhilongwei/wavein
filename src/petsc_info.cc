#include "wavein/petsc_info.h"

#include <petscsys.h>

#include <array>

namespace wavein
{

std::string petsc_version()
{
    std::array<char, 256> version{};
    PetscGetVersion(version.data(), version.size());
    return version.data();
}

} // namespace wavein
