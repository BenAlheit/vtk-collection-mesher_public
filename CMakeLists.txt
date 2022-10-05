##
#  CMake script for the step-1 tutorial program:
##

# Set the name of the project and target:
SET(TARGET "vtk-collection-mesher")

# Declare all source files the target consists of. Here, this is only
# the one step-X.cc file, but as you expand your project you may wish
# to add other source files as well. If your project becomes much larger,
# you may want to either replace the following statement by something like
#    FILE(GLOB_RECURSE TARGET_SRC  "source/*.cc")
#    FILE(GLOB_RECURSE TARGET_INC  "include/*.h")
#    SET(TARGET_SRC ${TARGET_SRC}  ${TARGET_INC})
# or switch altogether to the large project CMakeLists.txt file discussed
# in the "CMake in user projects" page accessible from the "User info"
# page of the documentation.
#FILE(GLOB_RECURSE TARGET_SRC "*.cpp" "*.hpp" "*.h" "*.cc")
#SET(TARGET_SRC ${TARGET_SRC} ${TARGET_INC})

set(VTK_SMP_IMPLEMENTATION_TYPE "STDThread")
SET(TARGET_SRC
        ${TARGET}.cc
        ${TARGET_SRC}
        #  ${TARGET_INC}
        )

project(${TARGET} C CXX)
# Usually, you will not need to modify anything beyond this point...

CMAKE_MINIMUM_REQUIRED(VERSION 3.1.0)

FIND_PACKAGE(deal.II 9.3.0
        HINTS ${deal.II_DIR} ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR}
        )
IF(NOT ${deal.II_FOUND})
    MESSAGE(FATAL_ERROR "\n"
            "*** Could not locate a (sufficiently recent) version of deal.II. ***\n\n"
            "You may want to either pass a flag -DDEAL_II_DIR=/path/to/deal.II to cmake\n"
            "or set an environment variable \"DEAL_II_DIR\" that contains this path."
            )
ENDIF()

#set(JSON_BuildTests OFF CACHE INTERNAL "")
#add_subdirectory(external-libs/json)

find_package(VTK COMPONENTS
        CommonColor
        CommonCore
        CommonDataModel
        CommonTransforms
        CommonSystem
        FiltersGeneral
        FiltersModeling
        FiltersSources
        IOGeometry
        IOLegacy
        IOPLY
        IOXML
        InteractionStyle
        RenderingContextOpenGL2
        RenderingCore
        RenderingFreeType
        RenderingGL2PSOpenGL2
        RenderingOpenGL2
        )


if (NOT VTK_FOUND)
    message(FATAL_ERROR "ReadAllPolyDataTypes: Unable to find the VTK build folder.")
endif()




DEAL_II_INITIALIZE_CACHED_VARIABLES()
PROJECT(${TARGET})
DEAL_II_INVOKE_AUTOPILOT()
#add_executable(${TARGET} MACOSX_BUNDLE ${TARGET}.cc )
target_link_libraries(${TARGET} PRIVATE ${VTK_LIBRARIES})
vtk_module_autoinit(
        TARGETS ${TARGET}
        MODULES ${VTK_LIBRARIES}
)