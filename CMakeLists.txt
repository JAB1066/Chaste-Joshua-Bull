# Example CMake build configuration file for a Chaste user project.
# See also https://chaste.cs.ox.ac.uk/trac/wiki/ChasteGuides/CmakeBuildGuide for full documentation.

# You can set which Chaste components (or other projects) your project depends on by editing the
# find_package() call for Chaste.  E.g.
#   find_package(Chaste COMPONENTS cell_based)
# for a project just using the cell_based component (and its dependencies), or
#   find_package(Chaste COMPONENTS heart lung)
# for a project combining heart & lung simulations.
#
# Note that the order in which components are specified does not matter.

# Here we just depend on core components (nothing application-specific).
find_package(Chaste COMPONENTS cell_based)

# Alternatively, to specify a Chaste installation directory use a line like that below.
# This is needed if your project is not contained in the projects folder within a Chaste source tree.
#find_package(Chaste COMPONENTS heart crypt PATHS /path/to/chaste-install NO_DEFAULT_PATH)

# Change the project name in the line below to match the folder this file is in,
# i.e. the name of your project.
chaste_do_project(JoshuaBull)

# Add some custom test targets so we can just run on JoshuaBull code
# Stolen from MicrovesselChaste documentation
if (Chaste_PROFILE_GPROF OR Chaste_PROFILE_GPERFTOOLS)
    if (Chaste_PROFILE_GPERFTOOLS)
        set(extension gif)
    else()
        set(extension gmon)
    endif()
    
    set(CTEST_COMMAND ctest)
    add_custom_target(profile_project_JoshuaBull
        COMMAND ${CTEST_COMMAND} "-L" "project_JoshuaBull" "--output-on-failure"
        COMMAND ${PYTHON_EXECUTABLE} "${Chaste_SOURCE_DIR}/cmake/process_profile.py" "${Chaste_PROFILE_OUTPUT_DIR}" ${extension}
        DEPENDS Profile 
        WORKING_DIRECTORY ${Chaste_BINARY_DIR}
        VERBATIM)
endif()
