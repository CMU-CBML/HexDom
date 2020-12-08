# - Try to find Matlab
# Once done this will define
#  
#  Matlab_FOUND        - system has Matlab
#  Matlab_INCLUDE_DIR  - the Matlab include directory
#  Matlab_LIBS			- Matlab libs, debug and optimized

IF (Matlab_INCLUDE_DIR)
	# Already in cache, be silent
	SET(Matlab_FIND_QUIETLY 1)
ENDIF (Matlab_INCLUDE_DIR)

if( WIN32 )
    FIND_PATH( Matlab_INCLUDE_DIR libmatlbmx.mlib
              PATHS  
			  "C:\\Program Files\\MATLAB\\R2015a\\extern\\include" 
			  )
    
    FIND_PATH(Matlab_LIBRARY_DIRcd #? where do you use this variable?
	        NAMES libfixedpoint.lib
            PATHS 
			"C:\\Program Files\\MATLAB\\R2015a\\extern\\lib\\win64\\microsoft" 
			)
	FIND_PATH(Matlab_DLIBRARY_DIR #? where do you use this variable?
	        NAMES tovideodevice.dll
            PATHS 
			"C:\\Program Files\\MATLAB\\R2015a\\bin\\win64" 
			)
			
    set(Matlab_LIBS) #? what is the meaning of this set command? No value is set to Matlab_LIBS
    file (GLOB Matlab_LIBS "${Matlab_LIBRARY_DIR}/*.lib")	 
	MESSAGE(STATUS "Looking for Matlab - found in " ${Matlab_INCLUDE_DIR})
	MESSAGE(STATUS "Looking for Matlab lib - found in " ${Matlab_LIBRARY_DIR}) #? where did you get Matlab_LIBRARY_DIR
ELSE()#Linux
    #FIND_PATH( Matlab_INCLUDE_DIR libmatlbmx.mlib
    #          PATHS  
	#		  #"C:\\Program Files\\MATLAB\\R2015a\\extern\\include" 
    #          /opt/packages/matlab/R2016a/extern/include
	#		  )
    #
    #FIND_PATH(Matlab_LIBRARY_DIR #Matlab_LIBRARY_DIRcd
	#        NAMES libfixedpoint.lib
    #        PATHS 
	#		#"C:\\Program Files\\MATLAB\\R2015a\\extern\\lib\\win64\\microsoft" 
    #        /opt/packages/matlab/R2016a/extern/lib/glnxa64
	#		)
	#FIND_PATH(Matlab_DLIBRARY_DIR
	#        NAMES tovideodevice.so
    #        PATHS 
	#		#"C:\\Program Files\\MATLAB\\R2015a\\bin\\win64" 
    #        /opt/packages/matlab/R2016a/bin/glnxa64
	#		)
    
    set(Matlab_INCLUDE_DIR "/opt/packages/matlab/R2016a/extern/include")
    set(Matlab_LIBRARY_DIR "/opt/packages/matlab/R2016a/extern/lib/glnxa64")
    set(Matlab_DLIBRARY_DIR "/opt/packages/matlab/R2016a/bin/glnxa64")
			
    set(Matlab_LIBS) 
    file (GLOB Matlab_LIBS "${Matlab_LIBRARY_DIR}/*.lib")	 
	MESSAGE(STATUS "Looking for Matlab - found in " ${Matlab_INCLUDE_DIR})
	MESSAGE(STATUS "Looking for Matlab lib - found in " ${Matlab_LIBRARY_DIR})
ENDIF( WIN32 )
   
IF (Matlab_INCLUDE_DIR AND Matlab_LIBRARY_DIR)
	set(Matlab_FOUND 1)
	set(Matlab_LIBS ${Matlab_LIBS} CACHE STRING "")
	if (NOT Matlab_FIND_QUIETLY)
		MESSAGE(STATUS "Looking for Matlab - found in " ${Matlab_LIBRARY_DIR})
	endif()
ELSE ()
	SET( Matlab_FOUND 0 )
ENDIF ()

