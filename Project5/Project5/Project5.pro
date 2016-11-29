CONFIG -= qt

INCLUDEPATH += /uio/hume/student-u68/trygvels/armadillo_install/usr/local/include
LIBS += -L /uio/hume/student-u68/trygvels/armadillo_install/usr/local/lib64
LIBS += -llapack -lblas -larmadillo

SOURCES += main.cpp \


