QT += core gui \
      widgets

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets printsupport

TARGET = ClosedLoop
TEMPLATE = app

CONFIG += console

SOURCES += main.cpp\
           mainwindow.cpp \
           qcustomplot.cpp \
           plotwindow.cpp \
    hilbert.cpp \      
    sockstream.cpp \
    appconnect.cpp \    
    settings.cpp \
    filt.cpp \
    exludech.cpp

HEADERS  += mainwindow.h \
            qcustomplot.h \
            plotwindow.h \
    hilbert.h \        
    sockstream.h \
    appconnect.h \
    settings.h \
    filt.h \
    exludech.h

FORMS    += mainwindow.ui \
            plotwindow.ui \
    settings.ui \
    exludech.ui

LIBS += -lwsock32 -lws2_32 -mthreads

INCLUDEPATH += $$PWD/.
DEPENDPATH += $$PWD/.

RESOURCES += \
    resf.qrc
