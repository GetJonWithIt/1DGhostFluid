QT += core
QT -= gui

TARGET = 1DGhostFluidProject
CONFIG += console
CONFIG -= app_bundle

TEMPLATE = app

SOURCES += main.cpp \
    statevector.cpp \
    equationofstate.cpp \
    hllcsolver.cpp \
    wavespeeds.cpp \
    linearalgebra.cpp \
    solvers.cpp \
    godunovsolver.cpp \
    exactsolver.cpp \
    musclsolver.cpp \
    slopelimiter.cpp \
    ghostfluidmethod.cpp \
    multimaterialsystem.cpp \
    tests.cpp

HEADERS += \
    statevector.h \
    equationofstate.h \
    hllcsolver.h \
    wavespeeds.h \
    linearalgebra.h \
    solvers.h \
    godunovsolver.h \
    exactsolver.h \
    musclsolver.h \
    slopelimiter.h \
    ghostfluidmethod.h \
    multimaterialsystem.h \
    tests.h

