#ifndef PAINTSCENE_H
#define PAINTSCENE_H

#include <QGraphicsScene>
#include <QGraphicsSceneMouseEvent>
#include <QTimer>
#include <QDebug>
#include <plotwindow.h>
#include <QSound>
#include "paintform.h"

const int TMAX=1024;

class plotwindow;
class MainWindow;
class paintform;

class paintScene : public QGraphicsScene
{

    Q_OBJECT

public:

    plotwindow* pw;  
    MainWindow* mww;
    paintform* paintf;

    QTimer* tim;

    int rc, gc, bc, ac;
    int t0, t1, pos, counter, length, pointnum;
    int* data;
    int* angles;
    QPointF* linecoords;
    double sumd, meand;
    Complex t[TMAX];
    CArray cdata;
    QColor drcolor, fxcolor;
    bool randcolor, fixcolor, freqcolor, drawbpoints, startedline;
    QPointF previousPoint, p1, p2;

    void init(plotwindow* pww, MainWindow* mw);
    void getdata(int x);
    void drawline();
    void setcolor();
    void clearlast();
    bool eventFilter(QObject *target, QEvent *event);
    explicit paintScene(QObject *parent = 0);
    ~paintScene();

private:
    void drawBackground(QPainter *p, const QRectF &rect);

private:
    void mousePressEvent(QGraphicsSceneMouseEvent * event);
    void mouseMoveEvent(QGraphicsSceneMouseEvent *event);
    void mouseReleaseEvent(QGraphicsSceneMouseEvent *event);
    void mouseScrEvent(QGraphicsSceneMouseEvent *event);

private slots:
    void timerUpdate();

};

#endif // PAINTSCENE_H
