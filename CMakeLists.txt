cmake_minimum_required(VERSION 3.18)

# Find includes in corresponding build directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)

# Instruct CMake to run moc automatically when needed
set(CMAKE_AUTOMOC ON)

# Create code from a list of Qt designer ui files
set(CMAKE_AUTOUIC ON)

set(PROJECT_NAME tv9k)

project(${PROJECT_NAME})

# Not sure why exactly I need this, but I was getting a linking error otherwise
#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -luuid")

message(STATUS "${PROJECT_NAME}")

add_compile_options(-std=c++11)

find_package(Qt5Core REQUIRED)
find_package(Qt5Widgets REQUIRED)
find_package(Qt5OpenGL REQUIRED)
find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)
#find_package(OpenMP REQUIRED)
#FIND_PACKAGE(Boost  REQUIRED)



#if(VTK_DIR)
    #find_package(VTK 8.2 REQUIRED)
    #ADD_DEFINITIONS(-DVTK_DIR)
#endif(VTK_DIR)

#if(TTKBase_DIR)
    ## This should be added as a preprocessor definition automaticaly, but I had to do it manually.
    #ADD_DEFINITIONS(-DTTKBase_DIR)
    #find_package(TTKBase REQUIRED)
#endif(TTKBase_DIR)

set(SOURCES
    ./src/Data.cpp
    ./src/PlotWidget.cpp
    ./src/TracerVisualiserMain.cpp
    ./src/TracerVisualiserWidget.cpp
    ./src/TracerVisualiserWindow.cpp
    ./src/ArcBall/Ball.cpp
    ./src/ArcBall/BallAux.cpp
    ./src/ArcBall/BallMath.cpp
    ./src/utility/Geometry.cpp
    )

#INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIR} )

#if(VTK_DIR)
    #include(${VTK_USE_FILE})
#endif(VTK_DIR)

add_executable(${PROJECT_NAME} ${SOURCES} ${HEADERS})

target_link_libraries(${PROJECT_NAME}
    PRIVATE
    ${Qt5Core_LIBRARIES}
    ${Qt5Widgets_LIBRARIES}
    ${Qt5OpenGL_LIBRARIES}
    ${OPENGL_LIBRARIES}  # required so that include is "-framework OpenGL" on macos instad of "-lOpenGL"
    ${GLUT_LIBRARIES}
    #${netCDF_LIBRARIES}
    #${Boost_LIBRARIES}
    #OpenMP::OpenMP_CXX
    )

#if(APPLE)
    ## Leif: this was the only way I could get compilation with netcdf-cxx4 to work on macos
    #target_link_libraries(${PROJECT_NAME} PUBLIC /usr/local/lib/libnetcdf-cxx4.dylib)
#else()
    ##target_link_libraries(${PROJECT_NAME} PUBLIC netcdf ncdf_c++4 )
#endif(APPLE)

# Is this only needed for apple?
#target_include_directories(${PROJECT_NAME}
  #PUBLIC
  #${netCDF_INCLUDE_DIR}
  #${OPENGL_INCLUDE_DIR}
  #${Qt5Core_INCLUDE_DIRS}
#)

#if(VTK_DIR)
    #target_link_libraries(${PROJECT_NAME} ${VTK_LIBRARIES})
#endif(VTK_DIR)

#if(TTKBase_DIR)
    #target_link_libraries(${PROJECT_NAME} ttk::base::baseAll)
#endif(TTKBase_DIR)

# Is this only needed for apple?
#set_target_properties(${PROJECT_NAME}
    #PROPERTIES
    #AUTOMOC ON
    #AUTOUIC ON
#)

#add_custom_target(
        #clangformat
        #COMMAND /usr/bin/clang-format-10
        #-i
        #../src/*.cpp
        #../src/*.h
        #../src/*.hpp
        #../src/utility/*.cpp
        #../src/utility/*.h
        #../src/utility/*.hpp
#)


# for debugging cmake variables:
#get_cmake_property(_variableNames VARIABLES)
#list (SORT _variableNames)
#foreach (_variableName ${_variableNames})
    #message(STATUS "${_variableName}=${${_variableName}}")
#endforeach()
