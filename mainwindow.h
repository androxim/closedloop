#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <plotwindow.h>
#include <hilbert.h>
#include "appconnect.h"

class plotwindow;
class appconnect;
class Settings;

namespace Ui {

class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:  
    hilbert* ht;
    QString daqport;
    bool cleardata;
    explicit MainWindow(QWidget *parent = 0);      
    void adddata(string s, QString fpath);
    void printdata(QString str);
    void savelogtofile(QString str);
    ~MainWindow();

private slots:    
    void on_pushButton_clicked();
    void on_pushButton_3_clicked();

    void on_pushButton_2_clicked();

    void on_pushButton_4_clicked();

    void OpenDataFile();

    void on_pushButton_5_clicked();

    void on_pushButton_6_clicked();

    void on_checkBox_clicked();

private:
    Ui::MainWindow *ui;
    plotwindow *plotw;
    appconnect *connectWin;

};

#endif // MAINWINDOW_H
