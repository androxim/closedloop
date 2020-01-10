#include "paintscene.h"
#include "random"
#include "QtAlgorithms"
#include <QSound>
#include "qmath.h"

// TO DO:
// 1. correct angle on 45 degree
// 2. add transparancy control

paintScene::paintScene(QObject *parent) : QGraphicsScene(parent)
{
    qDeleteAll(this->items());
    pos=0;
    t0=0;
    randcolor=false;
    drawbpoints=false;
    startedline = false;
    fxcolor=Qt::green;
    linecoords = new QPointF[500];
    angles = new int[500];
    pointnum = 0;
  //  qApp->installEventFilter(this);
}

bool paintScene::eventFilter(QObject *target, QEvent *event)
{

}

paintScene::~paintScene()
{

}

void paintScene::init(plotwindow* pww, MainWindow* mw)
{
    pw=pww;
    mww=mw;
    mww->psstart=true;
    counter=1; length=512; sumd=0; t0=0;
    data = new int[length];
    tim = new QTimer(this);
    tim->connect(tim, SIGNAL(timeout()), this, SLOT(timerUpdate()));
    rc=0; gc=0; bc=0; ac=255;
 //   tim->setInterval(10);
 //   tim->start();
}


void paintScene::getdata(int x)
{
    if (counter==length)
    {
        counter=1;
        sumd=0;
    }
    else
    {
       // data[counter]=x;
        sumd+=x; meand=sumd/counter;
        counter++;
    }
    //if (x>meand*2)
     //   x=meand;
    t0=(x-meand)*paintf->eegsize;
    //QCursor::setPos(400,t0);
   // qDebug()<<x<<" "<<meand;
}

void paintScene::drawBackground(QPainter *p, const QRectF &rect)
{
   // QPixmap back("D:/brain0.jpg");
   // p->drawPixmap(rect.left(), rect.top(), back.scaled(1490,824));
}

void paintScene::timerUpdate()
{
}

void paintScene::setcolor()
{
    if (pw->start)
    {
        rc=pw->alpha*6;
        gc=pw->beta*6;
        bc=pw->theta*6;
       // rc = 50 + pw->theta*5;
       // gc = 250 - pw->alpha*2;
       // bc = 150 - pw->beta*2;
      /*  if ((pw->theta > pw->alpha) && (pw->theta>pw->beta))
        {
            rc=0;
            gc=0;
            bc=255;
        } else
        if ((pw->alpha > pw->theta) && (pw->alpha>pw->beta))
        {
            rc=0;
            gc=255;
            bc=0;
        } else
        if ((pw->beta > pw->theta) && (pw->beta>pw->alpha))
        {
            rc=255;
            gc=0;
            bc=0;
        }*/
        ac=198-pw->delta*3;
        // qDebug()<<rc<<" "<<gc<<" "<<bc<<" "<<ac;
    }
    if (rc>255) rc=255;
    if (gc>255) gc=255;
    if (gc<0) gc=29;
    if (bc>255) bc=255;
    if (bc<0) bc=30;
    if (ac<0) ac=25;
    if (freqcolor)
        drcolor=QColor(rc,gc,bc,ac);
    else if (fixcolor)
        drcolor=QColor(fxcolor.red(),fxcolor.green(),fxcolor.blue(),ac);
}

void paintScene::clearlast()
{
//    double angle;
    qDebug()<<pointnum;
    if ((pointnum>0))
    {
      //  angle = qAtan2(linecoords[pointnum].y()-linecoords[0].y(),linecoords[pointnum].x()-linecoords[0].x());
        for (int i=0; i<pointnum-1; i++)
        {
            addLine(linecoords[i].x(), linecoords[i].y(), linecoords[i+1].x()+qSin(angles[i])*data[i], linecoords[i+1].y()+qCos(angles[i])*data[i], QPen(QColor(255,255,255,255), paintf->pensize, Qt::SolidLine, Qt::RoundCap));
            addEllipse(linecoords[i].x(), linecoords[i].y(), paintf->pensize, paintf->pensize, QPen(Qt::NoPen), QBrush(QColor(255,255,255,255)));
        }
        pointnum=0;
    }
}

void paintScene::mouseReleaseEvent(QGraphicsSceneMouseEvent *event)
{
 //   if ((!drawbpoints) && (startedline))
 //       startedline=false;
}

void paintScene::mousePressEvent(QGraphicsSceneMouseEvent *event)
{
    if (((freqcolor) || (fixcolor)) && (!paintf->erasepen))
    {
        setcolor();
        addEllipse(event->scenePos().x(), event->scenePos().y(), paintf->pensize, paintf->pensize, QPen(Qt::NoPen), QBrush(drcolor));
    }
    else
    if ((randcolor) && (!paintf->erasepen))
        addEllipse(event->scenePos().x(), event->scenePos().y(), paintf->pensize, paintf->pensize, QPen(Qt::NoPen), QBrush(QColor(qrand() % 256, qrand() % 256, qrand() % 256, qrand() % 256)));
    else
        addEllipse(event->scenePos().x(), event->scenePos().y(), paintf->pensize, paintf->pensize, QPen(Qt::NoPen), QBrush(QColor(255,255,255,255)));
    previousPoint = event->scenePos();
  /*  if ((event->button() == Qt::LeftButton))
    {
        if ((!drawbpoints) && (!startedline))
        {
            startedline=true;
            pointnum=0;
        } else
        if ((!drawbpoints) && (startedline))
        {
            startedline=false;
            pointnum=0;
        }
    }*/
}

void paintScene::mouseMoveEvent(QGraphicsSceneMouseEvent *event)
{
   /* if ((!drawbpoints) && (startedline))
    {
        data[pointnum] = t0;
        linecoords[pointnum].setX(event->scenePos().x());
        linecoords[pointnum].setY(event->scenePos().y());
        angles[pointnum] = qAtan2(event->scenePos().y()-previousPoint.y(),event->scenePos().x()-previousPoint.x());
        pointnum++;
        //  qDebug()<<pointnum;
    }*/
    double angle = qAtan2(event->scenePos().y()-previousPoint.y(),event->scenePos().x()-previousPoint.x());
   // qDebug()<<qRadiansToDegrees(angle);
    if (((freqcolor) || (fixcolor)) && (!paintf->erasepen))
    {
        setcolor();
        addLine(previousPoint.x(), previousPoint.y(), event->scenePos().x()+qSin(angle)*t0, event->scenePos().y()+qCos(angle)*t0, QPen(drcolor, paintf->pensize, Qt::SolidLine, Qt::RoundCap));
    }
    else
    if ((randcolor) && (!paintf->erasepen))
        addLine(previousPoint.x(), previousPoint.y(), event->scenePos().x()+qSin(angle)*t0, event->scenePos().y()+qCos(angle)*t0, QPen(QColor(qrand() % 256, qrand() % 256, qrand() % 256, qrand() % 256), paintf->pensize, Qt::SolidLine, Qt::RoundCap));
    else
        addLine(previousPoint.x(), previousPoint.y(), event->scenePos().x(), event->scenePos().y(), QPen(QColor(255,255,255,255), paintf->pensize, Qt::SolidLine, Qt::RoundCap));
    previousPoint = event->scenePos();
}

void paintScene::drawline()
{
    int x0, x1, y0, y1;
    x0=p1.x(); y0=p1.y(); x1=previousPoint.x(); y1=previousPoint.y();
    double angle = qAtan2(y1-y0,x1-x0);
  //  qDebug()<<x0<<" "<<y0;
  //  qDebug()<<x1<<" "<<y1;
  //  qDebug()<<angle;
    int t1 = x1-x0; int t2=abs(y1-y0);
    int p = t1 / t2;
    int k=0;
    while (k<t1)
    {
        if (k==p)
        {
            y0++;
            p=p+p;
        }
        setcolor();
        addLine(x0, y0, x0+1+qSin(angle)*t0, y0+qCos(angle)*t0, QPen(QColor(rc,gc,bc,ac), paintf->pensize, Qt::SolidLine, Qt::RoundCap));
        k++;
        paintf->delay(1);
    }
   /* for (int i=x0; i<x1; i++)
    {
        setcolor();
        addLine(i, y0, i+1+qSin(angle)*t0, y1+qCos(angle)*t0, QPen(QColor(rc,gc,bc,ac), paintf->pensize, Qt::SolidLine, Qt::RoundCap));
        paintf->delay(1);
    }*/
}
