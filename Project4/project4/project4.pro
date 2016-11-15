
CONFIG -= qt

INCLUDEPATH += /uio/hume/student-u68/trygvels/armadillo_install/usr/local/include
LIBS += -L /uio/hume/student-u68/trygvels/armadillo_install/usr/local/lib64
LIBS += -llapack -lblas -larmadillo

SOURCES += main.cpp \
        ising.cpp \
    lib.cpp

HEADERS += ising.h \
    lib.h



INCLUDEPATH += /usr/lib64/openmpi/bin/mpirun

#Comment out at astro
#QMAKE_LFLAGS += -pthread -m64 -Wl,-rpath -Wl,/usr/lib64/openmpi/lib -Wl,--enable-new-dtags -L/usr/lib64/openmpi/lib -lmpi_cxx -lmpi

QMAKE_CFLAGS += $$system(mpicc --showme:compile)
QMAKE_LFLAGS += $$system(mpicxx --showme:link)
QMAKE_CXXFLAGS += $$system(mpicxx --showme:compile) -DMPICH_IGNORE_CXX_SEEK
QMAKE_CXXFLAGS_RELEASE += $$system(mpicxx --showme:compile) -DMPICH_IGNORE_CXX_SEEK
