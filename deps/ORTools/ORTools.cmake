# Google OR-Tools — CP-SAT constraint solver for Magma tube assignment.
# Built from source with all LP solvers disabled (only CP-SAT needed).
# BUILD_DEPS=ON lets OR-Tools self-bootstrap abseil, protobuf, re2, etc.

if(BUILD_SHARED_LIBS)
    set(_ortools_shared ON)
else()
    set(_ortools_shared OFF)
endif()

orcaslicer_add_cmake_project(ORTools
    URL "https://github.com/google/or-tools/archive/refs/tags/v9.15.tar.gz"
    URL_HASH SHA256=6395a00a97ff30af878ee8d7fd5ad0ab1c7844f7219182c6d71acbee1b5f3026
    CMAKE_ARGS
        -DCMAKE_POSITION_INDEPENDENT_CODE=ON
        -DBUILD_SHARED_LIBS=${_ortools_shared}
        -DBUILD_DEPS=ON
        -DBUILD_CXX=ON
        -DBUILD_PYTHON=OFF
        -DBUILD_JAVA=OFF
        -DBUILD_DOTNET=OFF
        -DBUILD_SAMPLES=OFF
        -DBUILD_EXAMPLES=OFF
        -DBUILD_TESTING=OFF
        -DBUILD_DOC=OFF
        -DUSE_SCIP=OFF
        -DUSE_COINOR=OFF
        -DUSE_GLPK=OFF
        -DUSE_HIGHS=OFF
        -DUSE_PDLP=OFF
)

if (MSVC)
    add_debug_dep(dep_ORTools)
endif ()
