find_package(HDF5 REQUIRED)
# find_package(Ldap REQUIRED)

LIST(APPEND sources "main-sacla-api2.cpp")
LIST(APPEND prepare-sources "prepare-cheetah-sacla-api2.cpp")

include_directories("include")
include_directories(${CHEETAH_INCLUDES})
include_directories(${HDF5_INCLUDE_DIR})
include_directories(${SACLA_API2_INCLUDE_DIR})

add_executable(cheetah-sacla-api2 ${sources})
add_dependencies(cheetah-sacla-api2 cheetah)

target_link_libraries(cheetah-sacla-api2 ldap) # ${LDAP_LIBRARIES})
target_link_libraries(cheetah-sacla-api2 ${CHEETAH_LIBRARY})
target_link_libraries(cheetah-sacla-api2 ${SACLA_API2_LIB})
target_link_libraries(cheetah-sacla-api2 ${MYSQL_LIB})

add_executable(prepare-cheetah-sacla-api2 ${prepare-sources})

target_link_libraries(prepare-cheetah-sacla-api2 ldap) # ${LDAP_LIBRARIES})
target_link_libraries(prepare-cheetah-sacla-api2 ${CHEETAH_LIBRARY})
target_link_libraries(prepare-cheetah-sacla-api2 ${SACLA_API2_LIB})
target_link_libraries(prepare-cheetah-sacla-api2 ${MYSQL_LIB})

install(TARGETS cheetah-sacla-api2
  RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
  LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib${LIB_SUFFIX}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib${LIB_SUFFIX})
install(TARGETS prepare-cheetah-sacla-api2
  RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
  LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib${LIB_SUFFIX}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib${LIB_SUFFIX})
install(PROGRAMS prepare-cheetah-sacla-api2.py
  RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
