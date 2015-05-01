find_package(HDF5 REQUIRED)

LIST(APPEND sources "main-sacla-api.cpp")

include_directories("include")
include_directories(${CHEETAH_INCLUDES})
include_directories(${HDF5_INCLUDE_DIR})
include_directories(${SACLA_API_INCLUDE_DIR})

add_executable(cheetah-sacla-api ${sources})

add_dependencies(cheetah-sacla-api cheetah)

target_link_libraries(cheetah-sacla-api ${CHEETAH_LIBRARY})
target_link_libraries(cheetah-sacla-api ${SACLA_API_LIB})
target_link_libraries(cheetah-sacla-api ${MYSQL_LIB})

install(TARGETS cheetah-sacla-api
  RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
  LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib${LIB_SUFFIX}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib${LIB_SUFFIX})
