# Source files
set(LIB_SOURCES ProjectData.cpp)

# Library
add_library(ApplicationsLib ${LIB_SOURCES})

target_link_libraries(ApplicationsLib
	MeshGeoToolsLib
	NumLib
	ProcessLib
	${VTK_LIBRARIES}
)
ADD_VTK_DEPENDENCY(ApplicationsLib)

if(TARGET Eigen)
	add_dependencies(ApplicationsLib Eigen)
endif()
if(TARGET Boost)
	add_dependencies(ApplicationsLib Boost)
endif()
