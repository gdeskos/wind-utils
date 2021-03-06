
add_library(utilscore
  PerfUtils.cpp
  CFDMesh.cpp)
target_link_libraries(utilscore
  ${Trilinos_TPL_LIBRARIES}
  ${Trilinos_LIBRARIES}
  ${YAML_CPP_LIBRARIES}
  ${MPI_LIBRARIES})

target_compile_options(utilscore
  PRIVATE $<$<COMPILE_LANGUAGE:CXX>:${WIND_UTILS_WARN_FLAGS}>)

if(MPI_COMPILE_FLAGS)
  set_target_properties(utilscore PROPERTIES
    COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif(MPI_COMPILE_FLAGS)

if(MPI_LINK_FLAGS)
  set_target_properties(utilscore PROPERTIES
    LINK_FLAGS "${MPI_LINK_FLAGS}")
endif(MPI_LINK_FLAGS)

install(TARGETS utilscore
  EXPORT "${PROJECT_NAME}Libraries"
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)

file(GLOB CORE_HEADERS *.h)
install(FILES ${CORE_HEADERS}
  DESTINATION include)
