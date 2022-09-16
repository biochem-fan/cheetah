find_package(HDF5 REQUIRED)

#LIST(APPEND sources "api-stub.cpp")
LIST(APPEND sources "main-sacla-online-api.cpp")


include_directories("include")
include_directories(${CHEETAH_INCLUDES})
include_directories(${HDF5_INCLUDE_DIR})
include_directories(${SACLA_ARRAY_API_INCLUDE_DIR})
include_directories(${SACLA_ONLINE_API_INCLUDE_DIR})

add_executable(cheetah-sacla-online-api ${sources})

add_dependencies(cheetah-sacla-online-api cheetah)

target_link_libraries(cheetah-sacla-online-api pthread)
target_link_libraries(cheetah-sacla-online-api ${CHEETAH_LIBRARY})
target_link_libraries(cheetah-sacla-online-api ${MYSQL_LIB})
target_link_libraries(cheetah-sacla-online-api ${SACLA_ARRAY_API_LIB})
target_link_libraries(cheetah-sacla-online-api ${SACLA_ONLINE_API_LIB})

install(TARGETS cheetah-sacla-online-api
  RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
  LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib${LIB_SUFFIX}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib${LIB_SUFFIX})
