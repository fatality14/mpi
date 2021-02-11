TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
        main.cpp

LIBS += -L$$PWD/libs/ -lmsmpi
LIBS += -L$$PWD/libs/ -lmsmpifec
LIBS += -L$$PWD/libs/ -lmsmpifmc

INCLUDEPATH += $$PWD/libs
DEPENDPATH += $$PWD/libs

HEADERS += \
    headers/mpi.h \
    headers/mpif.h \
    headers/mpifptr.h \
    headers/mpio.h \
    headers/mspms.h \
    headers/pmidbg.h

