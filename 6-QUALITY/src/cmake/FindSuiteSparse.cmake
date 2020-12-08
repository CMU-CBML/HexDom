# - Try to find SuiteSparse
# Once done this will define
#  
#  SuiteSparse_FOUND        - system has SuiteSparse
#  SuiteSparse_INCLUDE_DIR  - the SuiteSparse include directory
#  SuiteSparse_LIBS			- SuiteSparse libs, debug and optimized

IF (SuiteSparse_INCLUDE_DIR)
	# Already in cache, be silent
	SET(SuiteSparse_FIND_QUIETLY 1)
ENDIF (SuiteSparse_INCLUDE_DIR)

if( WIN32 )
    FIND_PATH( SuiteSparse_INCLUDE_DIR cholmod.h
              PATHS 
			  "C:\\libs\\win32\\SuiteSparse\\Include" 
			  "C:\\Users\\luotingf\\Dropbox\\opt\\SuiteSparse\\include"
			  "C:\\Users\\Luoting\\Dropbox\\opt\\SuiteSparse\\include"
			  "C:\\Users\\milan\\Videos\\Projects\\lib\\SuiteSparse\\include"
			  "C:\\Users\\milano\\Videos\\lib\\SuiteSparse\\include"
			  "C:\\Users\\Lei\\Videos\\lib\\SuiteSparse\\include"
			  usr/local/include/
			  )
    
    FIND_PATH(SuiteSparse_LIBRARY_DIR
	        NAMES libamd.lib amd.lib
            PATHS 
			"C:\\libs\\win32\\SuiteSparse\\libs" 
			"C:\\Users\\milan\\Videos\\Projects\\lib\\SuiteSparse\\lib"
			"C:\\Users\\milano\\Videos\\lib\\SuiteSparse\\include"
			"C:\\Users\\Lei\\Videos\\lib\\SuiteSparse\\include"
			${SuiteSparse_INCLUDE_DIR}/../lib/
			)

    set(SuiteSparse_LIBS)
    file (GLOB SuiteSparse_DEBUG_LIBS "${SuiteSparse_LIBRARY_DIR}/*-gd.lib")
    file (GLOB SuiteSparse_RELEASE_LIBS "${SuiteSparse_LIBRARY_DIR}/*.lib")	
    list (REMOVE_ITEM SuiteSparse_RELEASE_LIBS ${SuiteSparse_DEBUG_LIBS})

    foreach(_debug_lib ${SuiteSparse_DEBUG_LIBS})
    	list(APPEND SuiteSparse_LIBS debug ${_debug_lib} )
    endforeach()
    foreach(_release_lib ${SuiteSparse_RELEASE_LIBS})
    	list(APPEND SuiteSparse_LIBS optimized ${_release_lib} )
    endforeach()
else( WIN32 )
    FIND_PATH( SuiteSparse_INCLUDE_DIR cholmod.h
              PATHS 
			  usr/local/include/
			  )
    
    FIND_PATH(SuiteSparse_LIBRARY_DIR libamd.a
            PATHS 
			${SuiteSparse_INCLUDE_DIR}/../lib/
			)

 	set(SuiteSparse_LIBS)
    	file (GLOB SuiteSparse_LIBS "${SuiteSparse_LIBRARY_DIR}/*.a")
ENDIF( WIN32 )
   
IF (SuiteSparse_INCLUDE_DIR AND SuiteSparse_LIBS)
	set(SuiteSparse_FOUND 1)
	set(SuiteSparse_LIBS ${SuiteSparse_LIBS} CACHE STRING "")
	if (NOT SuiteSparse_FIND_QUIETLY)
		MESSAGE(STATUS "Looking for SuiteSparse - found in " ${SuiteSparse_LIBRARY_DIR})
	endif()
ELSE ()
	SET( SuiteSparse_FOUND 0 )
ENDIF ()

