# - Find MKL
# Find the native MKL headers and libraries.
#
#  MKL_INCLUDE_DIR  - where to find MKL.h, etc.
#  MKL_LIBRARY_DIR  - List of library dirs when using MKL, e.g. mkl, ipp, icl, intel threading
#  MKL_LIBRARIES    - List of libraries when using MKL.
#  MKL_FOUND        - True if MKL found.
#
# Copyright (c) 2012 Luoting Fu <luoting.fu@me.com>

IF (MKL_INCLUDE_DIR)
	# Already in cache, be silent
	SET(MKL_FIND_QUIETLY 1)
ENDIF (MKL_INCLUDE_DIR)

# Look for the header file.
if (NOT MKL_INCLUDE_DIR)
FIND_PATH(MKL_INCLUDE_DIR NAMES mkl.h 
	PATHS "$ENV{MKLROOT}\\include"
	"C:\\Program Files (x86)\\Intel\\Composer XE 2013\\mkl\\include"
	"C:\\Program Files (x86)\\Intel\\Composer XE 2011\\mkl\\include"
	)
endif()

# Grab the MKL Runtime library
if(MKL_INCLUDE_DIR)
	set(MKL_INCLUDE_DIR ${MKL_INCLUDE_DIR} CACHE STRING "")
	if(CMAKE_SIZEOF_VOID_P EQUAL 4)   # Regular x86
		SET(MKL_FOUND 0)
		MESSAGE(FATAL_ERROR "FindMKL only support the X64 version of MKL.")
	else()
		# find_library(MKL_LIBRARIES
			# NAMES mkl_rt mkl_rt.lib
			# PATHS 
			# $ENV{MKLROOT}/lib/intel64/
			# "C:\\Program Files (x86)\\Intel\\Composer XE 2011 SP1\\mkl\\lib\\intel64"
			# NO_DEFAULT_PATH
			# )
			
		find_path(MKL_MKL_LIB_DIR
			NAMES mkl_rt mkl_rt.lib
			PATHS 
			${MKL_INCLUDE_DIR}/../lib/intel64/
			$ENV{MKLROOT}/lib/intel64/
			"C:\\Program Files (x86)\\Intel\\Composer XE 2013\\mkl\\lib\\intel64"
			)
		
		find_path(MKL_IPP_LIB_DIR
			NAMES ipps ipps.lib
			PATHS 
			${MKL_INCLUDE_DIR}/../../ipp/lib/intel64/
			$ENV{IPPROOT}/lib/intel64/
			"C:\\Program Files (x86)\\Intel\\Composer XE 2013\\ipp\\lib\\intel64"
			)		
			
		find_path(MKL_ICL_LIB_DIR
			NAMES libiomp5md libiomp5md.lib
			PATHS 
			${MKL_INCLUDE_DIR}/../../compiler/lib/intel64/
			$ENV{ICPP_COMPILER12}/compier/lib/intel64/
			"C:\\Program Files (x86)\\Intel\\Composer XE 2013\\compiler\\lib\\intel64"
			)	
					
		set(MKL_LIBRARY_DIR ${MKL_MKL_LIB_DIR} ${MKL_IPP_LIB_DIR} ${MKL_ICL_LIB_DIR})
	endif()
endif()

if(MKL_INCLUDE_DIR AND MKL_LIBRARY_DIR)
	set(MKL_LIBRARIES 
		general mkl_intel_lp64_dll.lib 			# interface layer
		debug mkl_sequential_dll.lib			# threading layer
		optimized mkl_intel_thread_dll.lib 		# threading layer
		general mkl_core_dll.lib				# computational layer
		optimized libiomp5md.lib 				# runtime layer
		CACHE STRING "")
	set(MKL_LIBRARY_DIR ${MKL_LIBRARY_DIR} CACHE STRING "")
	SET(MKL_FOUND 1)
else()
	SET(MKL_FOUND 0)
endif()	
		
# Report the results.
IF(NOT MKL_FOUND)
	MESSAGE(FATAL_ERROR "${MKL_DIR_MESSAGE}")
ELSE (NOT MKL_FOUND)
	MARK_AS_ADVANCED(MKL_LIBRARY_DIR)
	IF(NOT MKL_FIND_QUIETLY)
		MESSAGE(STATUS "Looking for MKL - found")
	ENDIF(NOT MKL_FIND_QUIETLY)
ENDIF(NOT MKL_FOUND)
