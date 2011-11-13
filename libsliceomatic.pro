#-------------------------------------------------
#
# Project created by QtCreator 2011-11-13T00:07:06
#
#-------------------------------------------------

QT       += opengl svg

QT       -= gui

TARGET = libsliceomatic
TEMPLATE = lib

DEFINES += LIBSLICEOMATIC_LIBRARY

SOURCES += libsliceomatic.cpp \
    stl.cpp

HEADERS += libsliceomatic.h\
        libsliceomatic_global.h \
    stl.h \
    cgaldefs.h

symbian {
    #Symbian specific definitions
    MMP_RULES += EXPORTUNFROZEN
    TARGET.UID3 = 0xE76F15F9
    TARGET.CAPABILITY = 
    TARGET.EPOCALLOWDLLDATA = 1
    addFiles.sources = libsliceomatic.dll
    addFiles.path = !:/sys/bin
    DEPLOYMENT += addFiles
}

unix:!symbian {
    maemo5 {
        target.path = /opt/usr/lib
    } else {
        target.path = /usr/local/lib
    }
    INSTALLS += target
}
