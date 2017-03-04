elements_subdir(AstrometricCalibration 0.1)


# Declare dependencies
#elements_depends_on_subdirs(ElementsKernel)
#elements_depends_on_subdirs(EuclidDmBindings)
#elements_depends_on_subdirs(EuclidDmUtils)
#elements_depends_on_subdirs(EuclidMdbParameters)

#find_package(Boost)
#find_package(CCfits)
#find_package(Cfitsio)
#find_package(Eigen3)

# Declare the library
#elements_add_library(biascorrection src/lib/*.cpp LINK_LIBRARIES ElementsKernel EuclidDmBindings EuclidDmUtils EuclidMdbParameters Boost CCfits Cfitsio Eigen3 INCLUDE_DIRS ElementsKernel EuclidDmBindings EuclidDmUtils EuclidMdbParameters Boost CCfits Cfitsio Eigen3 PUBLIC_HEADERS biascorrection)

#elements_add_library(LINK_LIBRARIES ElementsKernel EuclidDmBindings EuclidDmUtils EuclidMdbParameters INCLUDE_DIRS ElementsKernel EuclidDmBindings EuclidDmUtils EuclidMdbParameters)
#elements_add_library(LINK_LIBRARIES ElementsKernel INCLUDE_DIRS ElementsKernel)

# Create the C++ program
#elements_add_executable(wsCppProgram src/program/WsProgram.cpp
#                        LINK_LIBRARIES WorkshopModule
#                        INCLUDE_DIRS WorkshopModule)

# Install the python modules
elements_install_python_modules()

elements_install_scripts()


# Install the configuration files
elements_install_conf_files()

# Create the Python program
#elements_add_python_program(AstrometricCalibration)

# C++ tests
#elements_add_unit_test(BoostPfAlgorithm_test                  
#                       tests/src/PfAlgorithm_test.cpp
#                       LINK_LIBRARIES WorkshopModule TYPE Boost )
