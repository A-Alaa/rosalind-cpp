################################################################################
### Copyrights (c) 2015
### Marwan Abdellah <abdellah.marwan@gmail.com>
################################################################################

# Qt5 moccing
set( CMAKE_INCLUDE_CURRENT_DIR ON )
set( CMAKE_AUTOMOC ON )

# Find each component independently.
find_package( Qt5Core REQUIRED  )
find_package( Qt5Widgets REQUIRED )
find_package( Qt5OpenGL REQUIRED )
find_package( Qt5Gui REQUIRED )
#find_package( Qt5Network REQUIRED )



include_directories( ${Qt5Core_INCLUDES} )
include_directories( ${Qt5Widgets_INCLUDES} )
include_directories( ${Qt5OpenGL_INCLUDES} )
include_directories( ${Qt5Gui_INCLUDES} )
#include_directories( ${Qt5Network_INCLUDES} )

link_libraries( Qt5::Core Qt5::Widgets Qt5::OpenGL Qt5::Gui )
