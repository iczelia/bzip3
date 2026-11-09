cmake_minimum_required(VERSION 3.13...3.31 FATAL_ERROR)

set(BZIP3_VERSION "" CACHE STRING "Version of bzip3 (default: taken from .tarball-version or git)")
set(BZIP3_FULL_VERSION "${BZIP3_VERSION}")
if(NOT BZIP3_FULL_VERSION AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/.tarball-version")
  file(STRINGS "${CMAKE_CURRENT_SOURCE_DIR}/.tarball-version" BZIP3_FULL_VERSION LIMIT_COUNT 1)
  string(STRIP "${BZIP3_FULL_VERSION}" BZIP3_FULL_VERSION)
endif()
if(NOT BZIP3_FULL_VERSION AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/.git")
  find_package(Git QUIET)
  if(GIT_FOUND)
    execute_process(
      COMMAND "${GIT_EXECUTABLE}" describe --tags --abbrev=7 HEAD
      WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
      RESULT_VARIABLE BZIP3_GIT_RESULT
      OUTPUT_VARIABLE BZIP3_GIT_DESCRIBE
      OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
    if(BZIP3_GIT_RESULT EQUAL 0 AND BZIP3_GIT_DESCRIBE MATCHES "^v?[0-9]")
      # 1.5.5-2-g2f82d62 => 1.5.5.r2-g2f82d62
      string(REGEX REPLACE "^v" "" BZIP3_FULL_VERSION "${BZIP3_GIT_DESCRIBE}")
      if(BZIP3_FULL_VERSION MATCHES "^([^-]*)-(.*)$")
        set(BZIP3_FULL_VERSION "${CMAKE_MATCH_1}.r${CMAKE_MATCH_2}")
      endif()
      execute_process(
        COMMAND "${GIT_EXECUTABLE}" update-index --refresh
        WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
        OUTPUT_QUIET ERROR_QUIET)
      execute_process(
        COMMAND "${GIT_EXECUTABLE}" diff-index --name-only HEAD
        WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
        OUTPUT_VARIABLE BZIP3_GIT_DIRTY
        OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
      if(BZIP3_GIT_DIRTY)
        string(APPEND BZIP3_FULL_VERSION "-dirty")
      endif()
    endif()
  endif()
  if(NOT BZIP3_FULL_VERSION)
    message(WARNING "Could not determine the version of bzip3 from git (are the tags missing?). "
                    "Set it with -DBZIP3_VERSION=x.y.z.")
    set(BZIP3_FULL_VERSION "UNKNOWN")
    set(BZIP3_NUMERIC_VERSION "0.0.0")
  endif()
endif()
# The leading numeric part is what CMake itself can work with.
if(BZIP3_FULL_VERSION MATCHES "^([0-9]+(\\.[0-9]+)?(\\.[0-9]+)?)")
  set(BZIP3_NUMERIC_VERSION "${CMAKE_MATCH_1}")
elseif(NOT BZIP3_NUMERIC_VERSION)
  message(
    FATAL_ERROR
      "No version information found. Please either:\n"
      "* identify the correct version with `echo $version > .tarball-version` or -DBZIP3_VERSION=$version,\n"
      "* or use the source packages instead of a repo archive,\n"
      "* or use a full Git clone (with tags).")
endif()
message(STATUS "bzip3 version: ${BZIP3_FULL_VERSION}")

project(
  bzip3
  VERSION ${BZIP3_NUMERIC_VERSION}
  DESCRIPTION "A better and stronger spiritual successor to BZip2"
  HOMEPAGE_URL "https://github.com/iczelia/bzip3"
  LANGUAGES C)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)

# bzip3 is several times slower without compiler optimisations, so don't default to an unoptimised build.
get_property(BZIP3_MULTI_CONFIG GLOBAL PROPERTY GENERATOR_IS_MULTI_CONFIG)
if(NOT BZIP3_MULTI_CONFIG AND NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug Release RelWithDebInfo MinSizeRel)
endif()

option(BUILD_SHARED_LIBS "Build libbz3 as a shared library" ON)
option(BZIP3_BUILD_APPS "Build bzip3 applications" ON)
if(WIN32 AND NOT (MINGW OR CYGWIN))
  set(BZIP3_ENABLE_PTHREAD_DEFAULT OFF)
else()
  set(BZIP3_ENABLE_PTHREAD_DEFAULT ON)
endif()
option(BZIP3_ENABLE_PTHREAD "Enable use of pthread library" "${BZIP3_ENABLE_PTHREAD_DEFAULT}")
option(BZIP3_ENABLE_ARCH_NATIVE "Enable CPU-specific optimizations" OFF)
option(BZIP3_ENABLE_STATIC_EXE "Enable static builds of the executable" OFF)

include(CheckCCompilerFlag)
include(GNUInstallDirs)
# Defines the BUILD_TESTING option (ON by default) and enables `ctest`.
include(CTest)

set(CMAKE_INSTALL_RPATH ${CMAKE_INSTALL_FULL_LIBDIR})
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# Variables substituted into bzip3.pc.
set(prefix ${CMAKE_INSTALL_PREFIX})
set(exec_prefix ${CMAKE_INSTALL_PREFIX})
set(bindir ${CMAKE_INSTALL_FULL_BINDIR})
set(libdir ${CMAKE_INSTALL_FULL_LIBDIR})
set(includedir ${CMAKE_INSTALL_FULL_INCLUDEDIR})
set(PACKAGE ${PROJECT_NAME})
set(PACKAGE_VERSION ${BZIP3_FULL_VERSION})
if(WIN32 AND BUILD_SHARED_LIBS)
  set(extra_cflags " -DBZIP3_DLL_IMPORT=1")
else()
  set(extra_cflags "")
endif()
configure_file(bzip3.pc.in ${CMAKE_CURRENT_BINARY_DIR}/bzip3.pc @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/bzip3.pc
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig)

if(BZIP3_ENABLE_PTHREAD)
  set(THREADS_PREFER_PTHREAD_FLAG TRUE)
  find_package(Threads REQUIRED)
endif()

# Definitions needed to compile the sources of bzip3 itself. They are deliberately not propagated to
# the users of the library.
set(BZIP3_PRIVATE_DEFINITIONS VERSION="${BZIP3_FULL_VERSION}")
if(BZIP3_ENABLE_PTHREAD)
  list(APPEND BZIP3_PRIVATE_DEFINITIONS PTHREAD)
endif()

set(BZIP3_NATIVE_FLAGS "")
if(BZIP3_ENABLE_ARCH_NATIVE)
  check_c_compiler_flag(-march=native CC_SUPPORT_MARCH_NATIVE_FLAG)
  check_c_compiler_flag(-mtune=native CC_SUPPORT_MTUNE_NATIVE_FLAG)
  if(CC_SUPPORT_MARCH_NATIVE_FLAG AND CC_SUPPORT_MTUNE_NATIVE_FLAG)
    set(BZIP3_NATIVE_FLAGS -march=native -mtune=native)
  else()
    message(
      FATAL_ERROR
        "Compiler does not support native optimizations, disable `BZIP3_ENABLE_ARCH_NATIVE`"
    )
  endif()
endif()

add_library(bz3 src/libbz3.c)
target_compile_definitions(bz3 PRIVATE ${BZIP3_PRIVATE_DEFINITIONS})
target_compile_options(bz3 PRIVATE ${BZIP3_NATIVE_FLAGS})
if(WIN32 AND BUILD_SHARED_LIBS)
  target_compile_definitions(bz3 PRIVATE BZIP3_DLL_EXPORT=1)
  target_compile_definitions(bz3 INTERFACE BZIP3_DLL_IMPORT=1)
endif()
target_include_directories(
  bz3
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
         $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src
          ${CMAKE_CURRENT_SOURCE_DIR}/3rdparty/libsais)
if(BZIP3_ENABLE_PTHREAD)
  target_link_libraries(bz3 PUBLIC Threads::Threads)
endif()
# Keep in sync with the libtool -version-info in Makefile.am.
set_target_properties(
  bz3
  PROPERTIES OUTPUT_NAME bzip3
             SOVERSION "1"
             VERSION "1.0.1"
             PUBLIC_HEADER include/libbz3.h
             C_VISIBILITY_PRESET hidden)
if(BUILD_SHARED_LIBS)
  set_target_properties(bz3 PROPERTIES POSITION_INDEPENDENT_CODE ON)
endif()
install(
  TARGETS bz3
  EXPORT ${PROJECT_NAME}-targets
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(
  EXPORT ${PROJECT_NAME}-targets
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}"
  NAMESPACE ${PROJECT_NAME}::)
configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/${PROJECT_NAME}-config.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}-config.cmake"
  @ONLY)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}-config.cmake"
        DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}")

if(BZIP3_BUILD_APPS)
  add_executable(bzip3 src/main.c)
  target_compile_definitions(bzip3 PRIVATE ${BZIP3_PRIVATE_DEFINITIONS})
  target_compile_options(bzip3 PRIVATE ${BZIP3_NATIVE_FLAGS})
  if(BZIP3_ENABLE_STATIC_EXE)
    if(BUILD_SHARED_LIBS)
      message(
        FATAL_ERROR
          "libbz3 is not built as a static library, disable `BUILD_SHARED_LIBS`"
      )
    endif()
    check_c_compiler_flag(-static CC_SUPPORT_STATIC_FLAG)
    if(CC_SUPPORT_STATIC_FLAG)
      target_link_options(bzip3 PRIVATE -static)
    else()
      message(
        FATAL_ERROR
          "Compiler does not support static linking, disable `BZIP3_ENABLE_STATIC_EXE`"
      )
    endif()
  endif()
  target_link_libraries(bzip3 PRIVATE bz3)
  install(TARGETS bzip3 RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

  set(BZIP3_APP_SCRIPTS bunzip3 bz3cat bz3grep bz3less bz3more bz3most)
  install(PROGRAMS ${BZIP3_APP_SCRIPTS} DESTINATION ${CMAKE_INSTALL_BINDIR})

  if(UNIX)
    set(BZIP3_MANS
        bunzip3.1
        bz3cat.1
        bz3grep.1
        bz3less.1
        bz3more.1
        bz3most.1
        bzip3.1)
    string(TIMESTAMP MAN_DATE "%d %B %Y" UTC)
    set(TRANSFORMED_PACKAGE_NAME ${PROJECT_NAME})
    set(VERSION ${BZIP3_FULL_VERSION})
    foreach(BZIP3_MAN ${BZIP3_MANS})
      if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${BZIP3_MAN}.in)
        configure_file(${CMAKE_CURRENT_SOURCE_DIR}/${BZIP3_MAN}.in
                       ${CMAKE_CURRENT_BINARY_DIR}/${BZIP3_MAN} @ONLY)
      else()
        configure_file(${CMAKE_CURRENT_SOURCE_DIR}/${BZIP3_MAN}
                       ${CMAKE_CURRENT_BINARY_DIR}/${BZIP3_MAN} COPYONLY)
      endif()
      install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${BZIP3_MAN}
              DESTINATION ${CMAKE_INSTALL_MANDIR}/man1)
    endforeach()
  endif()
endif()

# The test suite; run it with `ctest`.
if(BUILD_TESTING AND CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
  add_executable(bz3-test-api tests/api.c)
  target_link_libraries(bz3-test-api PRIVATE bz3)
  if(BZIP3_ENABLE_PTHREAD)
    target_compile_definitions(bz3-test-api PRIVATE PTHREAD)
  endif()
  add_test(NAME api COMMAND bz3-test-api)

  find_program(BZIP3_SH sh)
  if(BZIP3_SH AND NOT CMAKE_CROSSCOMPILING)
    if(BZIP3_BUILD_APPS)
      add_test(NAME cli
               COMMAND ${BZIP3_SH} ${CMAKE_CURRENT_SOURCE_DIR}/tests/cli.sh
                       $<TARGET_FILE:bzip3>
                       ${CMAKE_CURRENT_SOURCE_DIR}/tests/data)
    endif()
    add_test(NAME symbols
             COMMAND ${BZIP3_SH} ${CMAKE_CURRENT_SOURCE_DIR}/tests/symbols.sh
                     $<TARGET_FILE:bz3>)
    # Exit status 77: the test could not be run (no `nm`).
    set_tests_properties(symbols PROPERTIES SKIP_RETURN_CODE 77)
  endif()
endif()
