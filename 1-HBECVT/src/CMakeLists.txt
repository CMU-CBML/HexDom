cmake_minimum_required (VERSION 2.8)
project (Segmentation)

# add our path and macros to cmake module path
set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake)
include (Common)
##############################################################################
# Mathematical libraries
##############################################################################
find_package (Eigen3 REQUIRED 3.0.0) 
if (NOT EIGEN3_INCLUDE_DIR)
  message (FATAL_ERROR "Eigen3 not found!")
endif ()
include_directories(${EIGEN3_INCLUDE_DIR})
add_definitions(${EIGEN3_DEFINITIONS})

if (WIN32)
  add_definitions(-D_USE_MATH_DEFINES -DNOMINMAX)
endif ()

find_package(OpenMP)
if (OPENMP_FOUND)
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
	set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
	add_definitions(-DUSE_OPENMP)
else ()
	set (OPENMP_NOTFOUND 1)
endif ()

macro(source_group_by_dir source_files)
 if(MSVC)
  set(sgbd_cur_dir ${CMAKE_CURRENT_SOURCE_DIR})
  foreach(sgbd_file ${${source_files}})
   string(REGEX REPLACE ${sgbd_cur_dir}/\(.*\) \\1 sgbd_fpath ${sgbd_file})
   string(REGEX REPLACE "\(.*\)/.*" \\1 sgbd_group_name ${sgbd_fpath})
   string(COMPARE EQUAL ${sgbd_fpath} ${sgbd_group_name} sgbd_nogroup)
   string(REPLACE "/" "\\" sgbd_group_name ${sgbd_group_name})
   if(sgbd_nogroup)
    set(sgbd_group_name "\\")
   endif(sgbd_nogroup)
   source_group(${sgbd_group_name} FILES ${sgbd_file})
  endforeach(sgbd_file)
 endif(MSVC)
endmacro(source_group_by_dir)


AUX_SOURCE_DIRECTORY(cpp SOURCES)
SET(SOURCES ${SOURCES})
FILE(GLOB_RECURSE project_headers cpp/*.h)
FILE(GLOB_RECURSE project_cpps cpp/*.c*)
set(all_files ${project_headers} ${project_cpps})
source_group_by_dir(all_files)

###################################################################################################
# Add special directory of files which cannot be generally added by the previous method
###################################################################################################
FILE (GLOB PCHFiles "${CMAKE_SOURCE_DIR}/stdafx.*")
source_group ("Precompiled Headers" FILES ${PCHFiles})

##############################################################################
# Application specific preprocessor definitions
##############################################################################
IF(MSVC)
add_definitions(-wd4018 -wd4267 -wd4101 -wd4190 -wd4244 -wd4305 -wd4099) # Turn off warnings
ENDIF(MSVC)

add_definitions(-D_CRT_SECURE_NO_WARNINGS)
#INCLUDE_DIRECTORIES(Entrance Constructor MeshStructure SparseSolver ANN CoMISo NURBS Analysis ConstructInclude Wap MeshIO)
ADD_EXECUTABLE(Segmentation ${all_files})
#target_link_libraries(TTSP_3D ${SuiteSparse_LIBS} ${GoTools_LIBS} ${Boost_LIBRARIES} ${Matlab_LIBS} ) 
#target_link_libraries(conv3D ${Matlab_LIBS} ) 
IF(MSVC)
SET_TARGET_PROPERTIES(Segmentation PROPERTIES COMPILE_FLAGS "/MP")
ELSE(MSVC)
add_definitions(-boost_system -boost_filesystem)
ENDIF(MSVC)
 

