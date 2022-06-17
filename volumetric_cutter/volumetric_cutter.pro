CINOLIB_PATH    = $$PWD/../lib/cinolib/
QT             += core opengl
TEMPLATE        = app
TARGET          = volumetric_cutter
CONFIG         += c++11
CONFIG         -= app_bundle
INCLUDEPATH    += $$CINOLIB_PATH/include
INCLUDEPATH    += $$CINOLIB_PATH/external/eigen
DEFINES        += CINOLIB_USES_QT
DEFINES        += CINOLIB_USES_GRAPH_CUT
DEFINES        += CINOLIB_USES_OPENGL
QMAKE_CXXFLAGS += -Wno-deprecated-declarations

# enable Tetgen (used in cinolib/tetgen_wrap.cpp)
DEFINES     += CINOLIB_USES_TETGEN
DEFINES     += TETLIBRARY
INCLUDEPATH *= /usr/local/include
LIBS        += -L/usr/local/lib -ltet

SOURCES += main.cpp
SOURCES += finalization.cpp
SOURCES += mesh_smoother.cpp
SOURCES += mesh_extractor.cpp
SOURCES += export_helper.cpp
SOURCES += subdivision_helper.cpp
SOURCES += batch.cpp
SOURCES += polyhedral_decomposition.cpp
SOURCES += cut.cpp
SOURCES += gui.cpp
SOURCES += loops.cpp

HEADERS += gui.h
HEADERS += finalization.h
HEADERS += mesh_smoother.h
HEADERS += mesh_extractor.h
HEADERS += export_helper.h
HEADERS += subdivision_helper.h
HEADERS += labeling.h
HEADERS += batch.h
HEADERS += polyhedral_decomposition.h
HEADERS += cut.h
HEADERS += definitions.h
HEADERS += state.h
HEADERS += undoredo.h
HEADERS += loops.h
