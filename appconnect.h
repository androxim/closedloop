/*
 * $BEGIN_BCI2000_LICENSE$
 * 
 * This file is part of BCI2000, a platform for real-time bio-signal research.
 * [ Copyright (C) 2000-2012: BCI2000 team and many external contributors ]
 * 
 * BCI2000 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 * 
 * BCI2000 is distributed in the hope that it will be useful, but
 *                         WITHOUT ANY WARRANTY
 * - without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * $END_BCI2000_LICENSE$
 */
#ifndef MAINWIDGET_H
#define MAINWIDGET_H

#include <string>
#include <map>

#include <QtWidgets>
#include <QThread>

#include "SockStream.h"
#include "plotwindow.h"

//forward declare our qt classes to decrease compile time
class QGridLayout;
class QGroupBox;
class QHBoxLayout;
class QLabel;
class QLineEdit;
class QComboBox;
class QPushButton;
class QVBoxLayout;
class QCheckBox;
class QButtonGroup;
class QForm;
class QTableWidget;

/*
Our main GUI class will be called mainUI. It inherits from the Qt
QWidget class, allowing it to use slots/signals, layouts, and
everything a form would need.
*/

class plotwindow;

class appconnect : public QWidget
{
  //this macro MUST be declared for any object which will use signals/slots
  //it provides the ability to use the private slots: class declaration

  Q_OBJECT

public:
  bool ready;
  int chnum, totalch;
  int count;
  double** datach;
  plotwindow* wd;  
  appconnect();
  ~appconnect();
  void connectButCallback();
  void clearstates();

private slots:
  void timerUpdate();

private:
  // QThread descendant
  friend struct ReceivingThread;
  struct ReceivingThread : QThread
  {
    ReceivingThread( appconnect& p ) : parent( p ) {}
    void run();
    appconnect& parent;
  } mReceivingThread;

private:
  //although Qt Designer allows us to graphically design our UI, it is helpful
  // to know what is going on in the actual code. Therefore, this class will manually
  // setup the forms rather than using the Qt Designer ui file.
  // layout ===========================================================
  QVBoxLayout *mainLayout;

  //connection objects
  QGroupBox *connectGroup;
  QVBoxLayout *connectLayout;
  QLabel *receiveLabel;
  QLineEdit *receiveBox;
  QLabel *sendLabel;
  QLineEdit *sendBox; 
  QLineEdit *statusBox;

  //table objects
  QGroupBox *statesGroup;
  QGridLayout *statesLayout;
  QTableWidget *stateTable;

  //layout setup functions
  void setupUI();
  void setupConnectionGroup();
  void setupStatesGroup();

  // data members ====================================================
  std::string address;
  //these are the objects used to connect to bci2000
  receiving_udpsocket recSocket;
  sockstream recConnection;
  sending_udpsocket sendSocket;
  sockstream sendConnection;
  //this map holds the state names and associated values
  std::map<std::string, std::string> states;
  QMutex mMutex;
  //a timer to update the display
  QTimer* timer;

};
#endif
