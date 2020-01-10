#ifndef PAINTFORM_H
#define PAINTFORM_H

#include <QWidget>
#include <QGraphicsScene>
#include <paintscene.h>
#include "plotwindow.h"

class paintScene;
class plotwindow;

namespace Ui {
class paintform;
}

class paintform : public QWidget
{
    Q_OBJECT

public:
    explicit paintform(QWidget *parent = 0);
    paintScene *scene;
    plotwindow *pw;
    QPixmap pm;
    int eegsize, pensize;
    bool erasepen;
    bool eventFilter(QObject *target, QEvent *event);
    void delay(int temp);
    ~paintform();


private slots:
    void on_pushButton_clicked();

    void on_checkBox_clicked();

    void on_pushButton_2_clicked();

    void on_spinBox_valueChanged(int arg1);

    void on_spinBox_2_valueChanged(int arg1);

    void on_pushButton_3_clicked();

    void on_pushButton_4_clicked();

    void on_checkBox_2_clicked();

    void on_radioButton_clicked();

    void on_radioButton_2_clicked();

    void on_radioButton_3_clicked();

private:
    Ui::paintform *ui;
    QPointF previousPoint;

    void resizeEvent(QResizeEvent *);
    void showEvent(QShowEvent *);
    void fitView();

};

#endif // PAINTFORM_H
