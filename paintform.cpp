#include "paintform.h"
#include "ui_paintform.h"
#include "paintscene.h"
#include <QColorDialog>

paintform::paintform(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::paintform)
{
    ui->setupUi(this);
    scene = new paintScene();
    scene->paintf=this;
    scene->randcolor=false; scene->fixcolor=false; scene->freqcolor=true;
    ui->graphicsView->setFixedSize(1425,720);
    ui->graphicsView->setScene(scene);
    ui->pushButton_3->setGeometry(700,740,100,30);
    ui->pushButton->setGeometry(850,740,80,30);
    ui->pushButton_4->setGeometry(950,740,80,30);
    ui->pushButton_2->setGeometry(1070,740,125,30);
    ui->checkBox->setGeometry(550,735,110,20);
    ui->checkBox_2->setGeometry(550,755,110,20);
    ui->radioButton->setGeometry(550,735,120,20);
    ui->radioButton_2->setGeometry(550,755,120,20);
    ui->radioButton_3->setGeometry(550,775,145,20);
    ui->radioButton->setChecked(scene->randcolor);
    ui->radioButton_2->setChecked(scene->fixcolor);
    ui->radioButton_3->setChecked(scene->freqcolor);
    ui->label->setGeometry(380,748,100,20);
    ui->spinBox->setGeometry(475,748,40,20);
    ui->label_2->setGeometry(260,748,60,20);
    ui->spinBox_2->setGeometry(320,748,40,20);
    ui->checkBox->setChecked(scene->randcolor);
    ui->checkBox->setVisible(false);
    ui->checkBox_2->setChecked(scene->drawbpoints);
    ui->checkBox_2->setVisible(false);
    eegsize=1; pensize=3;
    ui->spinBox->setValue(eegsize);
    ui->spinBox_2->setValue(pensize);
    erasepen=false;
    ui->graphicsView->installEventFilter(this);
}

paintform::~paintform()
{
    delete ui;
}


void paintform::resizeEvent(QResizeEvent *)
{
    fitView();
}

void paintform::showEvent(QShowEvent *)
{
    fitView();
}

bool paintform::eventFilter(QObject *target, QEvent *event)
{
    QKeyEvent* keyEvent = static_cast<QKeyEvent*>(event);

    if ((target == ui->graphicsView) && (event->type() == QEvent::KeyPress))
    {
        if (keyEvent->key() == Qt::Key_Z)
        {
            if (eegsize>1)
                eegsize--;
            ui->spinBox->setValue(eegsize);
        }
        if (keyEvent->key() == Qt::Key_X)
        {
            if (eegsize<12)
                eegsize++;
            ui->spinBox->setValue(eegsize);
        }
        if (keyEvent->key() == Qt::Key_A)
        {
            if (pensize>1)
                pensize--;
            ui->spinBox_2->setValue(pensize);
        }
        if (keyEvent->key() == Qt::Key_S)
        {
            if (pensize<50)
                pensize++;
            ui->spinBox_2->setValue(pensize);
        }

    }

    if ((target == ui->graphicsView) && (event->type() == QEvent::MouseButtonPress))
    {
        QMouseEvent* mouseEvent = static_cast<QMouseEvent*>(event);
        if (mouseEvent->button() == Qt::MiddleButton)
        {
            erasepen=!erasepen;
        }
    }

    if ((target == ui->graphicsView) && (event->type() == QEvent::MouseMove))
    {
      //  qDebug()<<scene->pointnum;
      //  QMouseEvent* mouseEvent = static_cast<QMouseEvent*>(event);
      //  if ((scene->drawbpoints) && (scene->startedline))
      //  {
      //      scene->linecoords[ scene->pointnum].setX(mouseEvent->pos().x());
      //      scene->linecoords[ scene->pointnum].setY(mouseEvent->pos().y());
     //       scene->pointnum++;
     //       qDebug()<<scene->pointnum;
      //  }
    }

    return false;
}

void paintform::fitView()
{
    const QRectF rect = QRectF(1,1, ui->graphicsView->width(), ui->graphicsView->height());
    ui->graphicsView->fitInView(rect, Qt::KeepAspectRatio);
    ui->graphicsView->setSceneRect(rect);
}

void paintform::on_checkBox_clicked()
{
    if (!ui->checkBox->isChecked())
        scene->randcolor=false;
    else
        scene->randcolor=true;
}

void paintform::on_pushButton_clicked()
{
    if (!pm.isNull())
        scene->addPixmap(pm.scaled(1425,724));
    else
        scene->clear();
}

void paintform::on_pushButton_2_clicked()
{
    QString filename=QFileDialog::getOpenFileName(this,tr("Open File"),"D://","All files (*.*)");
    if (filename!="")
    {
        pm.load(filename);
        scene->addPixmap(pm.scaled(1425,724));
    }
}

void paintform::on_spinBox_valueChanged(int arg1)
{
    eegsize=arg1;
}

void paintform::on_spinBox_2_valueChanged(int arg1)
{
    pensize=arg1;
}

void paintform::delay(int temp)
{
    QTime dieTime = QTime::currentTime().addMSecs(temp);
    while (QTime::currentTime() < dieTime)
        QCoreApplication::processEvents(QEventLoop::AllEvents, 100);
}

void paintform::on_pushButton_3_clicked()
{
    scene->fxcolor = QColorDialog::getColor(Qt::green, this, "Select Color", QColorDialog::DontUseNativeDialog);

    //scene->clearlast();
    //scene->drawline();
   // int t=0;
  //  while (t<200)
   // {
   //     QCursor::setPos(500, 100 + t);
    //    delay(10);
   //     t++;
   // }
}

void paintform::on_pushButton_4_clicked()
{
    QString fileName=QFileDialog::getSaveFileName(this, "Save image", QCoreApplication::applicationDirPath(), "BMP Files (*.bmp);;JPEG (*.JPEG);;PNG (*.png)" );
    if (!fileName.isNull())
    {
        QPixmap pixMap = this->ui->graphicsView->grab();
        pixMap.save(fileName);
    }
}

void paintform::on_checkBox_2_clicked()
{
    if (!ui->checkBox_2->isChecked())
        scene->drawbpoints=false;
    else
        scene->drawbpoints=true;
}

void paintform::on_radioButton_clicked()
{
    scene->randcolor=true;
    scene->fixcolor=false;
    scene->freqcolor=false;

}

void paintform::on_radioButton_2_clicked()
{
    scene->randcolor=false;
    scene->fixcolor=true;
    scene->freqcolor=false;

}

void paintform::on_radioButton_3_clicked()
{
    scene->randcolor=false;
    scene->fixcolor=false;
    scene->freqcolor=true;

}
