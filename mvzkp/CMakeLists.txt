add_library(mvzkp INTERFACE)

target_include_directories(mvzkp
INTERFACE
  ${CMAKE_CURRENT_SOURCE_DIR}/..
)

target_link_libraries(mvzkp
INTERFACE
  mv-qsvole
)
