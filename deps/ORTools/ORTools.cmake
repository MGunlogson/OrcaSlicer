# Google OR-Tools — CP-SAT constraint solver for Magma tube assignment.
# Built from source with all LP solvers disabled (only CP-SAT needed).
# BUILD_DEPS=ON lets OR-Tools self-bootstrap abseil, protobuf, re2, etc.

if(BUILD_SHARED_LIBS)
    set(_ortools_shared ON)
else()
    set(_ortools_shared OFF)
endif()

# Magma uses ONLY CP-SAT from OR-Tools — strip everything else.
# CP-SAT's irreducible deps are: abseil, protobuf, re2, eigen3, zlib, bzip2.
#
# Flatpak-friendly handling:
# OR-Tools' BUILD_DEPS=ON uses FetchContent to pull abseil/protobuf/re2 from
# GitHub at build time. Flatpak's sandbox blocks all network during build,
# so we point FetchContent at directories pre-populated by the Flatpak
# manifest's `type: git` sources. zlib/bzip2/eigen3 are provided by the
# GNOME SDK runtime, so we tell OR-Tools to use those via find_package
# instead of bundling — saves three more flatpak sources.
#
# When FLATPAK is OFF, these args are unset and OR-Tools' BUILD_DEPS=ON
# downloads everything normally — zero impact on Linux/Win/Mac builds.
set(_ortools_flatpak_args "")
if (FLATPAK)
    # OR-Tools' BUILD_DEPS=ON force-overrides BUILD_ZLIB/BZip2/Eigen3 back to ON
    # via CMAKE_DEPENDENT_OPTION, so we can't tell it "use system" — we have to
    # pre-populate its FetchContent sources. The Flatpak manifest stages all six
    # subdeps as git checkouts in external-packages/; the FETCHCONTENT_SOURCE_DIR
    # vars below point OR-Tools at those local copies so no network is needed
    # during the sandboxed build.
    list(APPEND _ortools_flatpak_args
        -DFETCHCONTENT_SOURCE_DIR_ABSL=${DEP_DOWNLOAD_DIR}/abseil-cpp
        -DFETCHCONTENT_SOURCE_DIR_PROTOBUF=${DEP_DOWNLOAD_DIR}/protobuf
        -DFETCHCONTENT_SOURCE_DIR_RE2=${DEP_DOWNLOAD_DIR}/re2
        -DFETCHCONTENT_SOURCE_DIR_ZLIB=${DEP_DOWNLOAD_DIR}/zlib
        -DFETCHCONTENT_SOURCE_DIR_BZIP2=${DEP_DOWNLOAD_DIR}/bzip2
        -DFETCHCONTENT_SOURCE_DIR_EIGEN3=${DEP_DOWNLOAD_DIR}/eigen
    )
endif ()

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
        # Magma uses CP-SAT only — drop FlatZinc to shrink the build.
        # (We tried -DBUILD_MATH_OPT=OFF too, but OR-Tools' Gurobi target has
        # no clean OFF switch and still links math_opt_proto, so disabling
        # MathOpt breaks the configure step. Leaving it ON is harmless: we
        # just don't use the resulting code.)
        -DBUILD_FLATZINC=OFF
        -DUSE_SCIP=OFF
        -DUSE_COINOR=OFF
        -DUSE_GLPK=OFF
        -DUSE_HIGHS=OFF
        -DUSE_PDLP=OFF
        ${_ortools_flatpak_args}
    DEPENDS dep_Eigen
)

# OR-Tools' BUILD_DEPS=ON builds AND INSTALLS its own eigen3 (3.4.0) into the shared dep
# prefix, overwriting the 5.0.1 that dep_Eigen installed ~40 build steps earlier. Nothing
# reports this; the slicer just fails to configure later with
#   "Could not find a configuration file for package Eigen3 compatible with 5.0.1"
# which points nowhere near OR-Tools. It was invisible until the upstream merge moved
# OrcaSlicer from Eigen 3.4.0 to 5.0.1 -- before that, both copies were the same version.
#
# There is no switch to stop it: INSTALL_BUILD_DEPS is declared in 9.15 but never referenced,
# and BUILD_Eigen3 is force-flipped back ON by CMAKE_DEPENDENT_OPTION whenever BUILD_DEPS
# is ON. Installing OR-Tools into a private prefix was tried and rejected: its bundled
# protobuf and absl are shared libraries, and moving them out of the shared lib dir breaks
# the AppImage bundler's runtime-dependency resolution (build_linux_image.sh).
#
# So let it install, then re-assert ours. Eigen is header-only, so this is a file copy.
# DEPENDS dep_Eigen above guarantees ours is built before this step runs.
# NOTE: if OR-Tools ever starts bundling another dep we also build, add it here. Today the
# only overlap is Eigen -- dep_ZLIB is not built on this platform.
ExternalProject_Add_Step(dep_ORTools restore_shared_eigen
    DEPENDEES install
    COMMENT "Restoring Eigen overwritten by OR-Tools' bundled copy"
    COMMAND ${CMAKE_COMMAND} --install ${CMAKE_BINARY_DIR}/dep_Eigen-prefix/src/dep_Eigen-build --prefix ${DESTDIR}
)

if (MSVC)
    add_debug_dep(dep_ORTools)
endif ()
