# - Try to find GoTools
# Once done this will define
#  
#  GoTools_FOUND        - system has GoTools
#  GoTools_INCLUDE_DIR  - the GoTools include directory
#  GoTools_LIBS			- GoTools libs, debug and optimized

IF (GoTools_INCLUDE_DIR)
	# Already in cache, be silent
	SET(GoTools_FIND_QUIETLY 1)
ENDIF (GoTools_INCLUDE_DIR)

if( WIN32 )
    FIND_PATH( GoTools_INCLUDE_DIR gotools_include_find.txt
              PATHS 
#			  "C:\\libs\\win32\\GoTools\\Include" 
#			  "C:\\Users\\luotingf\\Dropbox\\opt\\GoTools\\include"
#			  "C:\\Users\\Luoting\\Dropbox\\opt\\GoTools\\include"
#			  "C:\\Users\\milan\\Videos\\Projects\\lib\\GoTools\\include"
#			  "C:\\Users\\milano\\Videos\\lib\\GoTools_Rebuild\\include"
			  "C:\\Users\\weixd\\Documents\\lib\\GoTools_Rebuild\\include"
			  usr/local/include/
			  )
    
    FIND_PATH(GoTools_LIBRARY_DIR
	        NAMES GoToolsCore.lib
            PATHS 
#			"C:\\libs\\win32\\GoTools\\libs" 
#			"C:\\Users\\milan\\Videos\\Projects\\lib\\GoTools\\lib"
#			"C:\\Users\\milano\\Videos\\lib\\GoTools_mybuild\\lib"
			"C:\\Users\\weixd\\Documents\\lib\\GoTools_Rebuild\\lib"
			${GoTools_INCLUDE_DIR}/../lib/
			)

    set(GoTools_LIBS)
    file (GLOB GoTools_DEBUG_LIBS "${GoTools_LIBRARY_DIR}/*_debug.lib")
    file (GLOB GoTools_RELEASE_LIBS "${GoTools_LIBRARY_DIR}/*.lib")	
    list (REMOVE_ITEM GoTools_RELEASE_LIBS ${GoTools_DEBUG_LIBS})
	MESSAGE(STATUS "Looking for GoTools - found in " ${GoTools_INCLUDE_DIR})
	MESSAGE(STATUS "Looking for GoTools lib - found in " ${GoTools_LIBRARY_DIR})
    foreach(_debug_lib ${GoTools_DEBUG_LIBS})
    	list(APPEND GoTools_LIBS debug ${_debug_lib} )
    endforeach()
    foreach(_release_lib ${GoTools_RELEASE_LIBS})
    	list(APPEND GoTools_LIBS optimized ${_release_lib} )
    endforeach()
else( WIN32 )
    FIND_PATH( GoTools_INCLUDE_DIR cholmod.h
              PATHS 
		  usr/local/include/
			  )
    
    FIND_PATH(GoTools_LIBRARY_DIR libamd.a
            PATHS 
			${GoTools_INCLUDE_DIR}/../lib/
			)
    set(GoTools_LIBRARY_DIR "/usr/local/lib")
    set(GoTools_LIBS)
    file (GLOB GoTools_LIBS "${GoTools_LIBRARY_DIR}/libGo*") 	 

ENDIF( WIN32 )
   
IF (GoTools_INCLUDE_DIR AND GoTools_LIBRARY_DIR)
	set(GoTools_FOUND 1)
	set(GoTools_LIBS ${GoTools_LIBS} CACHE STRING "")
	if (NOT GoTools_FIND_QUIETLY)
		MESSAGE(STATUS "Looking for GoTools - found in " ${GoTools_LIBRARY_DIR})
	endif()
ELSE ()
	SET( GoTools_FOUND 0 )
ENDIF ()

