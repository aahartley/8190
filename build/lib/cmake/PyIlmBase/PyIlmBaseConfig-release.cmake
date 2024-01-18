#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "PyIlmBase::PyIex_Python3_8" for configuration "Release"
set_property(TARGET PyIlmBase::PyIex_Python3_8 APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(PyIlmBase::PyIex_Python3_8 PROPERTIES
  IMPORTED_LINK_DEPENDENT_LIBRARIES_RELEASE "Python3::Python"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libPyIex_Python3_8-2_5.so.25.0.7"
  IMPORTED_SONAME_RELEASE "libPyIex_Python3_8-2_5.so.25"
  )

list(APPEND _IMPORT_CHECK_TARGETS PyIlmBase::PyIex_Python3_8 )
list(APPEND _IMPORT_CHECK_FILES_FOR_PyIlmBase::PyIex_Python3_8 "${_IMPORT_PREFIX}/lib/libPyIex_Python3_8-2_5.so.25.0.7" )

# Import target "PyIlmBase::PyImath_Python3_8" for configuration "Release"
set_property(TARGET PyIlmBase::PyImath_Python3_8 APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(PyIlmBase::PyImath_Python3_8 PROPERTIES
  IMPORTED_LINK_DEPENDENT_LIBRARIES_RELEASE "PyIlmBase::PyIex_Python3_8;Python3::Python"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libPyImath_Python3_8-2_5.so.25.0.7"
  IMPORTED_SONAME_RELEASE "libPyImath_Python3_8-2_5.so.25"
  )

list(APPEND _IMPORT_CHECK_TARGETS PyIlmBase::PyImath_Python3_8 )
list(APPEND _IMPORT_CHECK_FILES_FOR_PyIlmBase::PyImath_Python3_8 "${_IMPORT_PREFIX}/lib/libPyImath_Python3_8-2_5.so.25.0.7" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
