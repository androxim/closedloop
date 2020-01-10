#ifndef SETTINGS_H
#define SETTINGS_H

#include <QWidget>
#include <plotwindow.h>
#include <exludech.h>

class plotwindow;

namespace Ui {
class Settings;
}

class Settings : public QWidget
{
    Q_OBJECT

public:  
    plotwindow* pwd;
    exludech* exch;
    double esavg, esprop;
    int eslength, siml;
    explicit Settings(QWidget *parent = 0);   
    ~Settings();
    void init();
    void ssdready();

private slots:
    void on_pushButton_2_clicked();

    void on_pushButton_clicked();

    void on_spinBox_valueChanged(int arg1);

    void on_spinBox_3_valueChanged(int arg1);   

    void on_pushButton_3_clicked();

    void on_checkBox_2_clicked();

    void on_checkBox_3_clicked(); 

    void on_checkBox_clicked();

    void on_toolButton_clicked();

    void on_checkBox_8_clicked();

    void on_pushButton_4_clicked(); 

    void on_spinBox_10_editingFinished();

    void on_spinBox_9_editingFinished();

    void on_checkBox_5_clicked();

private:
    Ui::Settings *ui;

};

#endif // SETTINGS_H
