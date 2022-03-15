add_subdirectory(tests)

add_library(minigia_libff INTERFACE)

target_include_directories(minigia_libff INTERFACE
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>
  $<INSTALL_INTERFACE:include/>  # <prefix>/include/mylib
)

