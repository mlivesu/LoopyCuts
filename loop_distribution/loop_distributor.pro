############################ PROJECT FILES ############################

include(libs.pri)


HEADERS       = glwidget.h \
                mesh_type.h \
                loop_stats.h \
                loop_mesher.h \
                loop_splitter.h \
                smooth_loops_mesh.h \
                tracing_field/inter_cross_validator.h \
                tracing_field/loop_finder.h \
                tracing_field/sharp_feature_sampler.h\
                tracing_field/sharp_feature_manager.h\
                tracing_field/loop_common_functions.h\
                tracing_field/triangle_mesh_functions.h\
                tracing_field/conflict_finder.h\
                tracing_field/stats_collector.h\
                tracing_field/artifact_removal.h\
                tracing_field/anisotropic_geodesic.h\
                tracing_field/graph/sink_geodesic_node.h\
                tracing_field/graph/singulatity_geodesic_node.h\
                tracing_field/graph/edge_geodesic_node.h\
                tracing_field/graph/base_geodesic_node.h\
                tracing_field/graph/anisotropic_graph.h\
                tracing_field/separatrix_parametrizer.h\
                tracing_field/graph/path_geodesic_node.h\
                tracing_field/graph/neigh_info.h\
                tracing_field/remesh/edge_splitter.h\
                tracing_field/remesh/edge_mesh_type.h\


SOURCES       = glwidget.cpp \
                main.cpp

QT           += opengl

############################ TARGET ############################

#App config
TARGET = loop_distributor

TEMPLATE = app
CONFIG += qt
CONFIG += c++11
CONFIG -= app_bundle

QT += core gui opengl xml widgets

#Debug/release optimization flags
CONFIG(debug, debug|release){
    DEFINES += DEBUG
}
CONFIG(release, debug|release){
    DEFINES -= DEBUG
    #just uncomment next line if you want to ignore asserts and got a more optimized binary
    CONFIG += FINAL_RELEASE
}

#Final release optimization flag
#FINAL_RELEASE {
#    unix:!macx{
#        QMAKE_CXXFLAGS_RELEASE -= -g -O2
#        QMAKE_CXXFLAGS += -O3 -DNDEBUG
#    }
#}

#macx {
#    QMAKE_MACOSX_DEPLOYMENT_TARGET = 10.13
#    QMAKE_MAC_SDK = macosx10.13
#}

############################ INCLUDES ############################

DEFINES += NO_PATCH_SIZING

INCLUDEPATH += $$VCGLIBDIR
INCLUDEPATH += $$ANTDIR/include
INCLUDEPATH += $$EIGENLIB

SOURCES += $$VCGLIBDIR/wrap/ply/plylib.cpp
SOURCES += $$VCGLIBDIR/wrap/gui/trackball.cpp
SOURCES += $$VCGLIBDIR/wrap/gui/trackmode.cpp
SOURCES += $$VCGLIBDIR/wrap/qt/anttweakbarMapperNew.cpp

# Awful problem with windows..
win32{
  DEFINES += NOMINMAX
  LIBS +=$$ANTDIR/lib/AntTweakBar.lib
}

mac{
# Mac specific Config required to avoid to make application bundles
# CONFIG -= app_bundle
  LIBS +=$$ANTDIR/lib/libAntTweakBar.dylib
  QMAKE_POST_LINK +="cp -P ../lib/AntTweakBar1.16/lib/libAntTweakBar.dylib . ; "
  QMAKE_POST_LINK +="install_name_tool -change ../lib/libAntTweakBar.dylib ./libAntTweakBar.dylib $$TARGET ; "
  DEPENDPATH += .
}

