cmake_minimum_required(VERSION 2.8)
project(hephaestusFoam)
add_executable(hephaestusFoam hephaestusFoam.C)

include_directories(${FOAM_INCLUDE_DIRS})
include_directories($ENV{FOAM_SRC}/finiteVolume/lnInclude)
include_directories($ENV{FOAM_SRC}/OpenFOAM/lnInclude)
include_directories($ENV{FOAM_SRC}/OSspecific/POSIX/lnInclude)
include_directories($ENV{FOAM_SRC}/OpenFOAM/dimensionedTypes/dimensionedScalar)


add_definitions(-D$ENV{WM_ARCH} -DWM_$ENV{WM_PRECISION_OPTION})
        
link_directories($ENV{FOAM_LIBBIN})
target_link_libraries(hephaestusFoam 
    ${FOAM_LIBBIN}/libfiniteVolume.so 
    ${FOAM_LIBBIN}/libOpenFOAM.so
    ${FOAM_LIBBIN}/libfoamGeneral.so
    ${FOAM_LIBBIN}/libdb.so)

install(TARGETS hephaestusFoam 
        RUNTIME DESTINATION /work/09629/m_hland10/ls6/apps/hephaestusFoam)