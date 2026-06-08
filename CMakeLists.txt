# ==============================================================================
# Configuration Initiale
# ==============================================================================
cmake_minimum_required(VERSION 3.15)
project(PhaseFieldFerro 
    VERSION 1.0.0 
    LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Type de build par defaut : Release" FORCE)
endif()

# ==============================================================================
# Dépendances 
# ==============================================================================
find_package(OpenMP REQUIRED)
find_package(Eigen3 REQUIRED) 

# ==============================================================================
# Architecture : Bibliothèque Centrale (CoreLib)
# ==============================================================================
set(CORE_SOURCES
    core/src/Mesh.cpp
    core/src/ResultsExporter.cpp
    core/src/PhaseFieldFerro.cpp
    core/src/Datafile.cpp
    core/src/Simulation.cpp
    core/src/Validation.cpp
)

add_library(CoreLib STATIC ${CORE_SOURCES})

# PUBLIC permet aux exécutables qui utilisent CoreLib de trouver les headers
target_include_directories(CoreLib PUBLIC core/include)

# Lier OpenMP et Eigen à la bibliothèque
if(OpenMP_CXX_FOUND)
    target_link_libraries(CoreLib PUBLIC OpenMP::OpenMP_CXX)
endif()
target_link_libraries(CoreLib PUBLIC Eigen3::Eigen) # <-- NOUVEAU : Lien avec Eigen

# Flags HPC appliqués à la bibliothèque
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(CoreLib PRIVATE -O3 -march=native -ffast-math -Wall -Wextra -Wpedantic)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    target_compile_options(CoreLib PRIVATE /O2 /arch:AVX2 /fp:fast /W4)
endif()

# ==============================================================================
# Exécutables
# ==============================================================================

# 1. Application Principale
add_executable(PhaseFieldFerro core/src/main.cpp)
target_link_libraries(PhaseFieldFerro PRIVATE CoreLib)

# 2. Unité de Tests
add_executable(UnitTests core/tests/test_physics.cpp)
target_link_libraries(UnitTests PRIVATE CoreLib)

enable_testing()
add_test(NAME PhysicsTest COMMAND UnitTests)

# ==========================================
# Cible personnalisée : Nettoyage des datas
# ==========================================
add_custom_target(clean-vtk
    # Note : remove_directory supprime le dossier entier, on ne met pas de *.vtk à la fin.
    COMMAND ${CMAKE_COMMAND} -E remove_directory "${CMAKE_BINARY_DIR}/data/VTK"
    COMMAND ${CMAKE_COMMAND} -E remove_directory "${CMAKE_SOURCE_DIR}/data/VTK"
    COMMAND ${CMAKE_COMMAND} -E remove_directory "${CMAKE_BINARY_DIR}/VTK"
    COMMENT "Nettoyage des dossiers de resultats VTK termine !"
)