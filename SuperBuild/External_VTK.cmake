if( NOT EXTERNAL_SOURCE_DIRECTORY )
  set( EXTERNAL_SOURCE_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/ExternalSources )
endif()

# Make sure this file is included only once by creating globally unique varibles
# based on the name of this included file.
get_filename_component(CMAKE_CURRENT_LIST_FILENAME ${CMAKE_CURRENT_LIST_FILE} NAME_WE)
if(${CMAKE_CURRENT_LIST_FILENAME}_FILE_INCLUDED)
  return()
endif()
set(${CMAKE_CURRENT_LIST_FILENAME}_FILE_INCLUDED 1)

## External_${extProjName}.cmake files can be recurisvely included,
## and cmake variables are global, so when including sub projects it
## is important make the extProjName and proj variables
## appear to stay constant in one of these files.
## Store global variables before overwriting (then restore at end of this file.)
ProjectDependancyPush(CACHED_extProjName ${extProjName})
ProjectDependancyPush(CACHED_proj ${proj})

# Make sure that the ExtProjName/IntProjName variables are unique globally
# even if other External_${ExtProjName}.cmake files are sourced by
# SlicerMacroCheckExternalProjectDependency
set(extProjName VTK) #The find_package known name
set(proj        VTK) #This local name

set(${extProjName}_REQUIRED_VERSION "7.0")

# Set dependency list
set(${proj}_DEPENDENCIES "")

# Include dependent projects if any
SlicerMacroCheckExternalProjectDependency(${proj})

if(NOT ( DEFINED "USE_SYSTEM_${extProjName}" AND "${USE_SYSTEM_${extProjName}}" ) )
  #message(STATUS "${__indent}Adding project ${proj}")

  # Set CMake OSX variable to pass down the external project
  set(CMAKE_OSX_EXTERNAL_PROJECT_ARGS)
  if(APPLE)
    list(APPEND CMAKE_OSX_EXTERNAL_PROJECT_ARGS
      -DCMAKE_OSX_ARCHITECTURES=${CMAKE_OSX_ARCHITECTURES}
      -DCMAKE_OSX_SYSROOT=${CMAKE_OSX_SYSROOT}
      -DCMAKE_OSX_DEPLOYMENT_TARGET=${CMAKE_OSX_DEPLOYMENT_TARGET}
      -DVTK_REQUIRED_OBJCXX_FLAGS:STRING="")
  endif()

  ### --- Project specific additions here
  # if(USE_VTKv6)
    list(APPEND EXTERNAL_PROJECT_OPTIONAL_ARGS
        -DVTK_Group_Qt:BOOL=ON
        -DModule_vtkTestingRendering:BOOL=ON
        )
  # endif()
  list(APPEND EXTERNAL_PROJECT_OPTIONAL_ARGS
      -DVTK_QT_VERSION:STRING=${ShapePopulationViewer_QT_VERSION}
      -DVTK_USE_QT:BOOL=ON
      -DVTK_USE_QVTK_QTOPENGL:BOOL=ON
      -DVTK_USE_GUISUPPORT:BOOL=ON
      -DQT_QMAKE_EXECUTABLE:FILEPATH=${QT_QMAKE_EXECUTABLE}
      )
  if(APPLE)
    list(APPEND EXTERNAL_PROJECT_OPTIONAL_ARGS
         -DVTK_USE_CARBON:BOOL=OFF
         -DVTK_USE_COCOA:BOOL=ON # Default to Cocoa, VTK/CMakeLists.txt will enable Carbon and disable cocoa if needed
         -DVTK_USE_X:BOOL=OFF
         -DVTK_REQUIRED_OBJCXX_FLAGS:STRING=
         #-DVTK_USE_RPATH:BOOL=ON # Unused
        )
  endif()
  if(UNIX AND NOT APPLE)
    find_package(FontConfig QUIET)
    if(FONTCONFIG_FOUND)
      list(APPEND EXTERNAL_PROJECT_OPTIONAL_ARGS
      -DModule_vtkRenderingFreeTypeFontConfig:BOOL=ON
    )
    endif()
  endif()

  set(${proj}_CMAKE_OPTIONS
      -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_CURRENT_BINARY_DIR}/${proj}-install
      -DBUILD_EXAMPLES:BOOL=OFF
      -DBUILD_TESTING:BOOL=OFF
      -DVTK_USE_PARALLEL:BOOL=ON
      -DVTK_USE_GL2PS:BOOL=ON
      -DVTK_DEBUG_LEAKS:BOOL=${${PROJECT_NAME}_USE_VTK_DEBUG_LEAKS}
      -DVTK_LEGACY_REMOVE:BOOL=OFF
      -DVTK_WRAP_TCL:BOOL=OFF
      #-DVTK_USE_RPATH:BOOL=ON # Unused
      -DVTK_WRAP_PYTHON:BOOL=OFF
      -DVTK_INSTALL_LIB_DIR:PATH=${${PROJECT_NAME}_INSTALL_LIB_DIR}
      ${EXTERNAL_PROJECT_OPTIONAL_ARGS}
    )
  ### --- End Project specific additions
  set(${proj}_GIT_TAG "v7.0.0")
  set(${proj}_REPOSITORY ${git_protocol}://vtk.org/VTK.git)
  ExternalProject_Add(${proj}
    GIT_REPOSITORY ${${proj}_REPOSITORY}
    GIT_TAG ${${proj}_GIT_TAG}
    SOURCE_DIR ${EXTERNAL_SOURCE_DIRECTORY}/${proj}
    BINARY_DIR ${proj}-build
    BUILD_COMMAND ${VTK_BUILD_STEP}
    LOG_CONFIGURE 0  # Wrap configure in script to ignore log output from dashboards
    LOG_BUILD     0  # Wrap build in script to to ignore log output from dashboards
    LOG_TEST      0  # Wrap test in script to to ignore log output from dashboards
    LOG_INSTALL   0  # Wrap install in script to to ignore log output from dashboards
    ${cmakeversion_external_update} "${cmakeversion_external_update_value}"
    CMAKE_GENERATOR ${gen}
    CMAKE_ARGS
      ${CMAKE_OSX_EXTERNAL_PROJECT_ARGS}
      ${COMMON_EXTERNAL_PROJECT_ARGS}
      ${${proj}_CMAKE_OPTIONS}
## We really do want to install in order to limit # of include paths INSTALL_COMMAND ""
    DEPENDS
      ${${proj}_DEPENDENCIES}
    )

  set(${extProjName}_DIR ${CMAKE_BINARY_DIR}/${proj}-install/lib/cmake/vtk-7.0)

else()
  if(${USE_SYSTEM_${extProjName}})
    find_package(${extProjName} ${${extProjName}_REQUIRED_VERSION} REQUIRED)
    message("USING the system ${extProjName}, set ${extProjName}_DIR=${${extProjName}_DIR}")
  endif()
  # The project is provided using ${extProjName}_DIR, nevertheless since other
  # project may depend on ${extProjName}, let's add an 'empty' one
  SlicerMacroEmptyExternalProject(${proj} "${${proj}_DEPENDENCIES}")
endif()

list(APPEND ${CMAKE_PROJECT_NAME}_SUPERBUILD_EP_VARS ${extProjName}_DIR:PATH)

ProjectDependancyPop(CACHED_extProjName extProjName)
ProjectDependancyPop(CACHED_proj proj)
