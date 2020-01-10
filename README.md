# tACS-EEG closed-loop system 

the code for an adaptive (alpha-phase dependent) non-invasive brain stimulation system with alternating currents

the code is written in C++, Qt 5.6.0, compiled with MinGW 4.9.2 32bit

the system provides with an options for EEG alpha phase analysis and prediction 

both in offline mode for testing data sets / methods and in online mode for experimental studies

==== hardware requirements ====

for performing an actual stimulation National Instruments DAQ-card 6229/6343 and neuroConn DC-STIMULATOR PLUS

are requiered, for a detailed configuration see file systemconfig.pdf 

==== software dependencies ====

National Instruments software and drivers should be installed

www.ni.com/en-us/support/downloads.html

==== third-party code / libraries licenses: ====

qcustomplot.h / qcustomplot.cpp codes by Emanuel Eichhammer, GNU General Public License v. 3.0    

appconnect.h / appconnect.cpp codes by J. Adam Wilson, BCI2000, GNU General Public License v. 3.0

sockstream.h / sockstream.cpp codes by Juergen Mellinger, GNU General Public License v. 3.0

Eigen library, Mozilla Public License v. 2.0

filt.h / filt.cpp codes by Mike Perkins, Cardinal Peak, LLC

NIDAQmx.h / NIDAQmx.cpp / NIDAQmx.lib, National Instruments
