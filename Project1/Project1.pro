TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt
LIBS += -llapack -lblas -larmadillo
SOURCES += main.cpp \
    problems.cpp \
    lib.cpp

HEADERS += \
    problems.h \
    lib.h
