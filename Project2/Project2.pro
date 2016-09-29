TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

INCLUDEPATH += /usr/local/include\

INCLUDEPATH += /uio/hume/student-u68/trygvels/armadillo-7.400.3/include

LIBS += -llapack -lblas

LIBS += -L/usr/local/lib

SOURCES += main.cpp \
    jacobiMethod.cpp

HEADERS += \
    jacobimethod.h
