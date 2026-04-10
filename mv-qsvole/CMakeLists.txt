add_library(mv-qsvole INTERFACE)

target_include_directories(mv-qsvole
INTERFACE
  ${CMAKE_CURRENT_SOURCE_DIR}/..
)

target_link_libraries(mv-qsvole
INTERFACE
  emp-zk
  emp-agmpc
)
