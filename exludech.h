#ifndef EXLUDECH_H
#define EXLUDECH_H

#include <QDialog>
#include "plotwindow.h"

namespace Ui {
class exludech;
}

class exludech : public QDialog
{
    Q_OBJECT

public:
    int totalchnumb;
    plotwindow *plw;
    explicit exludech(QWidget *parent = 0);
    ~exludech();
    void process();

private slots:
    void on_buttonBox_accepted();

    void on_pushButton_clicked();

    void on_pushButton_2_clicked();

    void setcheckboxs(bool fl);

private:
    Ui::exludech *ui;
};

#endif // EXLUDECH_H
