#include "mainwindow.h"
#include "plotwindow.h"
#include "ui_mainwindow.h"
#include "ui_plotwindow.h"
#include <QtGui>
#include <iostream>
#include <QMainWindow>
#include <QScrollArea>
#include <qdebug.h>
#include <hilbert.h>
#include "appconnect.h"
#include <QFileDialog>
#include <settings.h>
#include "NIDAQmx.h"

QStringList strList1;
QStringListModel *strListM1;

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{     
    setWindowFlags(Qt::WindowTitleHint | Qt::WindowMinimizeButtonHint);
    ui->setupUi(this);
    plotw = new plotwindow();              
    plotw->setFixedSize(1500,815);
    plotw->move(QApplication::desktop()->screen()->rect().center() - plotw->rect().center()-QPoint(0,30));
    plotw->start=false;  
    connectWin = new appconnect();
    connectWin->wd=plotw;
    strListM1 = new QStringListModel();
    strListM1->setStringList(strList1);
    cleardata=false;
    ui->checkBox->setChecked(cleardata);
    ui->listView->setModel(strListM1);
    connect(ui->actionOpen_Data_File, SIGNAL(triggered()), this, SLOT(OpenDataFile()));
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::adddata(string s, QString spath)
{
    ui->label->setText(QString::fromStdString(s));
    plotw->hnt=ht;
    for (int i=0; i<ht->npt; i++)
    {
        plotw->arrc.xc[i] = i;
        plotw->arrc.amp1[i]=ht->x[i];
    }
    ui->label_2->setText("Data File:  " + spath);
    plotw->setWindowTitle("Data File:  " + spath + "    offline EEG channel: " + QString::fromStdString(ht->channel) + "    Interval: " + QString::number(0) + " - " + QString::number(ht->npt));
}

void MainWindow::printdata(QString str)
{   
    strList1.push_front(str);
    strListM1->setStringList(strList1);
}

void MainWindow::on_pushButton_clicked()
{
    this->close();    
    QApplication::quit();
}

void MainWindow::on_pushButton_3_clicked()
{
    if (plotw->start)
       plotw->move(QApplication::desktop()->screen()->rect().center() - plotw->rect().center()-QPoint(0,15));
    if (plotw->cleardata)
        plotw->clean();
    plotw->show();
    plotw->doplot();   
    plotw->appcn=connectWin;    
}

void MainWindow::on_pushButton_2_clicked() // add sample_block_length
{       
    bool ok1,ok2;
    stringstream str1, str2;
    str1 << "Enter total number of EEG channels: ";
    int chnumb = QInputDialog::getInt(this,"",str1.str().c_str(), 1, 1, 64, 1, &ok1);
    if (ok1)
    {
        plotw->chnums=chnumb;
        plotw->exch=chnumb;
        connectWin->totalch=chnumb;
        for (int i=0; i<chnumb; i++)
            plotw->exlchannels[i]=1;
    }
    str2 << "Enter EEG source channel number: ";
    int channel = QInputDialog::getInt(this,"",str2.str().c_str(), 1, 1, plotw->chnums, 1, &ok2);
    if (ok2)
    {
        connectWin->chnum=channel-1;
        plotw->setWindowTitle(plotw->windowTitle()+"    Online EEG source channel: "+QString::number(channel));
        plotw->sourcech=channel-1;
        connectWin->clearstates();
        connectWin->show();
        connectWin->connectButCallback();     
    }
}

void MainWindow::on_pushButton_4_clicked()
{
    // add DAQ status check
//   DAQmxCreateTask("",&plotw->taskHandle);
 //  DAQmxCreateAOVoltageChan(plotw->taskHandle,"Dev1/ao2","",-5.0,5.0,DAQmx_Val_Volts,"");
 //  DAQmxCfgSampClkTiming(plotw->taskHandle,"",500,DAQmx_Val_Falling,DAQmx_Val_FiniteSamps,1500);
 //  DAQmxStartTask(plotw->taskHandle);
}

void MainWindow::OpenDataFile()
{
    QString filename=QFileDialog::getOpenFileName(this,tr("Open File"),"E://EEGdat//","Data file (*.dat);;All files (*.*)");
    if (filename!="")
    {
        if (plotw->start)
            plotw->close();
        ht->firstinit(filename,NMAX);
        string s = "EEG channel: " + ht->channel + "   Interval: " + to_string(ht->pos) + " - " + to_string(ht->npt);
        strList1.clear();
        strListM1->setStringList(strList1);
      //  plotw = new plotwindow();
      //  plotw->setFixedSize(1500,815);
       // plotw->start=fal
        connectWin->wd=plotw;
        plotw->appcn=connectWin;
        plotw->clean();
        plotw->loadstimfromfile(filename);
        adddata(s,filename);
        if (plotw->start)
        {
            plotw->updatedata();
            plotw->clearstim();
            plotw->clearfiltered();            
        }
     //   ht->imlength=100;
        plotw->show();
        plotw->doplot();
    }
}

void MainWindow::savelogtofile(QString str)
{
    QFile outputFile(str);
    outputFile.open(QIODevice::WriteOnly);
    QTextStream fout(&outputFile);
    for (int i = 0; i < strList1.size(); ++i)
          fout << strList1.at(i) << '\n';
    fout << '\n';
    outputFile.close();
}

void MainWindow::on_pushButton_5_clicked()
{
    strList1.clear();
    strListM1->setStringList(strList1);
}

void MainWindow::on_pushButton_6_clicked()
{
    QString fname=QFileDialog::getSaveFileName(this,tr("Save File as"),"E://","Text file (*.txt);;All files (*.*)");
    if (fname!="")
        savelogtofile(fname);
}

void MainWindow::on_checkBox_clicked()
{
    if (!ui->checkBox->isChecked())
        plotw->cleardata=false;
    else
    {
        plotw->hronar=false; plotw->zerocr=false; plotw->optstim=false; plotw->fftpredict=false;
        ht->imlength=750; ht->imstprop = 2; ht->pos=100; ht->extinterval=0.15; ht->lfilt=50;
        plotw->cleardata=true;
    }
}
