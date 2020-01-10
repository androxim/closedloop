#include "plotwindow.h"
#include "mainwindow.h"
#include "ui_plotwindow.h"
#include <iostream>
#include <QShortcut>
#include <hilbert.h>
#include <QFileDialog>
#include <settings.h>
#include <math.h>
#include <cmath>
#include <complex>
#include <stdio.h>
#include <NIDAQmx.h>
#include <filt.h>
#include <vector>
#include <exception>
#include <algorithm>
#include <iterator>
#include <random>
#include <chrono>
#include <thread>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <QTime>
#include <bits/stdc++.h>

/*
    QMessageBox msgBox;
    msgBox.setText("Something wrong!");
    msgBox.exec();
 */

// DONE:

// 1. plot optimal phase for last 2 Periods in pre-stim intervals and for first 2 Periods in post-stim intervals
// 2. degree deviaton between optimal phases for:
//      1) last 2 Periods in pre-stim interval and first 2 Periods in stim interval - actual accuracy of phase pred in real exp with all delays
//      2) last 2 Periods in stim interval and first 2 Periods in post-stim interval - entrainment traces
// 3. compare 2.1) and 2.2) in closed and opened eyes states
//     plot distribution in Python, check if there are any trends in Phase diffs, PLV within blocks or through blocks, check PLV inside sorted blocks
//     try different number of periods as intervals
//     same for Synchrony analysis
// 4. phase preservation index (PPI)
// 5. common spatial patterns (CSP)
// 6. ANOVA
// 7. ZC- and AR-optimization

const complex<double> I(0.0,1.0);

using namespace Eigen;

typedef std::vector<int> vectori;
typedef std::vector<double> vectord;
typedef std::complex<double> Complex;



plotwindow::plotwindow(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::plot)
{   
   // setWindowFlags(Qt::WindowTitleHint | Qt::WindowMinimizeButtonHint);    
    //datamap = QCPDataMap();
    //datamap_ptr = &datamap;
    QMenuBar* menu = new QMenuBar(this);
    QMenu* mainMenu = new QMenu("Menu");
    QMenu* mOpt = new QMenu("Options");
    QAction* acOpt = mOpt->addAction("Settings");   
    strLstM1 = new QStringListModel();
    strLstM1->setStringList(strLst1);
    strLstM2 = new QStringListModel();        
    strLstM2->setStringList(strLst2);
    strLstM3 = new QStringListModel();
    strLstM3->setStringList(strLst3);
    connect(acOpt, SIGNAL(triggered()),this,SLOT(slSettings()));
    menu->addMenu(mainMenu);
    menu->addMenu(mOpt);
    alltimesc = 0;
    start = false;
    counter=0; stims=0; startpos=0; ampval=70; tscale=12; recparts=0; chnums=64, sampleblock=10;
    offlinedata=false; addmode=false; addmodeon=false; adaptampl=false; hideardata=true;
    estimpredstim=true; filtereddata=false; correctdelay=true; filterstims=false;
    maxentropy=false; leastsqr=true; regforstim=false; hronarrun=false;
    phasepr=true; hronar=true; zerocr=true; optstim=false; fftpredict=true; // methods for stimulation
    usefiltering=false; filterardata=false; carfilter=false;
    recurbutter=false; filtfir=false; zerobutter=true; // filters
    applyssd=false; extractingssd=false; imagingafterstim=false;
    skipdelayedata=false; artifactremove=false; readypredict=false; randomstims=false;
    compavgdeviation=false; estimatenoisy=false; estimateclear=false; cleardata=false;
    randnoisestim=false; randphasestim=false;
    clearoptphase = new int [500];
    for (int i=0; i<500; i++)
        clearoptphase[i]=0;
    nstimforopt=15;
    daqscalevalue=180; // 250 mV
    stimshift=0;
    offsetcomp=0;
    transdelay=42; // 80 ms
    receivedelay=0;
    plv=0; averopttime=0; iti1=333; iti2=333; // 1000 - for 1.5 sec ITI, 500 for 1 sec, 333 for 0.5
    flength=25; fsampr=500; flowpass=12;
    butterord=2; lcutoff=6; hcutoff=15; // l 6, h 15
    autoreginput=new double[1000];   
    autoregdegree=75;
    totalar_error=0;
    autoregcoeffs=new double[500];
    fftphase = 0;
    use_mean_trial=false;
    exlchannels=new int[64]; exch=64;
    indexes = new int [64];    
    bc = vector<double>(butterord*2+1);
    ac = vector<double>(butterord*2+1);
    for (int i=0; i<64; i++)
        indexes[i]=0;
    rawdata = new double *[64]; // max channels number
    for (int i=0; i<chnums; i++)
    {
        rawdata[i]=new double[5000]; // max interval length
        for (int j=0; j<5000; j++)
            rawdata[i][j]=0;
    }
    resblock = new double *[200]; // stims in seq number
    for (int i=0; i<50; i++)
    {
        resblock[i]=new double[1500]; // length of pre + post intervals
        for (int j=0; j<1500; j++)
            resblock[i][j]=0;
    }
    arrres = new double *[20];
    for (int i=0; i<20; i++)
    {
        arrres[i]=new double[16]; // in-phase: ADD*4, ARA*4, PLV*4, SYN*4
        for (int j=0; j<16; j++)
            arrres[i][j]=0;
    }
    phaseshiftdiffs = new double *[500];
    for (int i=0; i<500; i++)
    {
        phaseshiftdiffs[i]=new double[2];
        for (int j=0; j<2; j++)
            phaseshiftdiffs[i][j]=0;
    }

    prestimplvs = new double *[20];
    for (int i=0; i<20; i++)
        prestimplvs[i] = new double[10];

    poststimplvs = new double *[20];
    for (int i=0; i<20; i++)
        poststimplvs[i] = new double[10];

    prestimsynchs = new double *[20];
    for (int i=0; i<20; i++)
        prestimsynchs[i] = new double[10];

    poststimsynchs = new double *[20];
    for (int i=0; i<20; i++)
        poststimsynchs[i] = new double[10];

    prestim = new double [500];
    for (int i=0; i<500; i++)
        prestim[i]=0;
    poststim = new double [500];
    for (int i=0; i<500; i++)
        poststim[i]=0;

    synchprestim = new double [500];
    for (int i=0; i<500; i++)
        synchprestim[i]=0;
    synchpoststim = new double [500];
    for (int i=0; i<500; i++)
        synchpoststim[i]=0;

    relativeacc =  new double [500];
    for (int i=0; i<500; i++)
        relativeacc[i]=0;

    allopttimes = new double [11000];
    opttimecount = 0;

    NumC = new double [20];
    DenC = new double [20];

    currentsubj=0;

    indexesforperms = new int[100];
    for (int i=0; i<100; i++)
        indexesforperms[i]=i;

    rsinitarray = new double *[100];
    for (int i=0; i<100; i++)
    {
        rsinitarray[i]=new double[500];
        for (int j=0; j<500; j++)
            rsinitarray[i][j]=0;
    }

    postinstphasediffs = new double *[100];
    for (int i=0; i<100; i++)
    {
        postinstphasediffs[i]=new double[5];
        for (int j=0; j<5; j++)
            postinstphasediffs[i][j]=0;
    }

    preinstphasediffs = new double *[100];
    for (int i=0; i<100; i++)
    {
        preinstphasediffs[i]=new double[5];
        for (int j=0; j<5; j++)
            preinstphasediffs[i][j]=0;
    }

    postppi = new double[8];
    for (int i=0; i<8; i++)
        postppi[i]=0;

    preppi = new double[8];
    for (int i=0; i<8; i++)
        preppi[i]=0;

    perms = new double[2000];

    blockparts=5;
    prestimtrials = new double[blockparts];
    poststimtrials = new double[blockparts];
    for (int i=0; i<blockparts; i++)
    {
        prestimtrials[i]=0;
        poststimtrials[i]=0;
    }

    arrpowerpre = new double[500];
    arrpowerpost = new double [500];
    arrplvpre = new double [500];
    arrplvpost = new double [500];
    for (int i=0; i<500; i++)
    {
        arrpowerpre[i]=0;
        arrpowerpost[i]=0;
        arrplvpre[i]=0;
        arrplvpost[i]=0;
    }

    meanprepower = new double [10];
    for (int i=0; i<10; i++)
        meanprepower[i]=0;

    meanpostpower = new double [10];
    for (int i=0; i<10; i++)
        meanpostpower[i]=0;

    meanpre_in = new double[10];
    for (int i=0; i<10; i++)
        meanpre_in[i]=0;

    meanpost_in = new double[10];
    for (int i=0; i<10; i++)
        meanpost_in[i]=0;

    meanpre_anti = new double[10];
    for (int i=0; i<10; i++)
        meanpre_anti[i]=0;

    meanpost_anti = new double[10];
    for (int i=0; i<10; i++)
        meanpost_anti[i]=0;

    meanpretrial = new double[500];
    meanposttrial = new double[500];

    meantrEOin = new double[1000];
    meantrEOanti = new double[1000];
    meantrECin = new double[1000];
    meantrECanti = new double[1000];

    tx = new double[750];
    ty = new double[750];
    for (int i=0; i<750; i++)
    {
        tx[i]=0;
        ty[i]=0;
    }
    fftframpl = new double[512];
    spatial_fitering = false;

   // hnt->imstprop=3;
   // ui->doubleSpinBox_3->setValue(hnt->imstprop);

    // SSD parameters
    butterord=2;
    sigb.lf=8; sigb.hf=12;
    noibp.lf=6; noibp.hf=14;
    samplenums=30000; samplesize=27; ssdcomp=1; maxssdcomp=1; ssdtime=5; ssdt=0; ssdscale=50, ssdshift=-100;

    arrc.xc = QVector<double>(NMAX);
    arrc.amp0 = QVector<double>(NMAX);
    arrc.amp1 = QVector<double>(NMAX);
    arrc.amp2 = QVector<double>(NMAX);
    arrc.of1 = QVector<double>(NMAX);
    arrc.of2 = QVector<double>(NMAX);
    arrc.of3 = QVector<double>(NMAX);
    arrc.of4 = QVector<double>(NMAX);
    arrc.of5 = QVector<double>(NMAX);
    arrc.of6 = QVector<double>(NMAX);
    ui->setupUi(this);    
    ui->spinBox->setMaximum(NMAX);
    ui->widget->installEventFilter(this);
    tim = new QTimer(this);
    tim->connect(tim, SIGNAL(timeout()), this, SLOT(timerUpdate()));
    tim->setInterval(0);
    ssdtim = new QTimer(this);
    ssdtim->connect(ssdtim, SIGNAL(timeout()), this, SLOT(ssdtimUpdate()));
    ssdtim->setInterval(0);
    estimation=false;   
    extractoptphase=false;
    showprediction=false;
    corrprednumb=0; totalstims=0;
    relations = new int[500];
    for (int i=0; i<500; i++)
        relations[i]=0;
    arrphasediffhr = new double[200];
    arrphasediffar = new double[200];
    arrphasediffzc = new double[200];
    arrphasedifffft = new double[200];
    arrsynchr = new double [200];
    arrsyncar = new double [200];
    arrsynczc = new double [200];
    arrsyncfft = new double [200];
    spectECpre = new double[36];
    spectECpost = new double[36];
    setaddmode(false);
    QTime time = QTime::currentTime();
    qsrand((uint)time.msec());
   // ssdt=0;
}

plotwindow::~plotwindow()
{
    delete ui;
}

bool plotwindow::eventFilter(QObject *target, QEvent *event)
{
    if ((target == ui->widget) && (event->type() == QEvent::MouseMove))
    {
        QMouseEvent* mouseEvent = static_cast<QMouseEvent*>(event);
        QPalette *p = new QPalette;
        p->setColor(QPalette::WindowText, Qt::darkBlue);
        ui->lcdNumber->setPalette(*p);
        ui->lcdNumber_2->setPalette(*p);
        ui->lcdNumber->display((int)ui->widget->xAxis->pixelToCoord(mouseEvent->pos().x()));
        ui->lcdNumber_2->display((int)ui->widget->yAxis->pixelToCoord(mouseEvent->pos().y()));      
    }
    if ((target == ui->widget) && (event->type() == QEvent::MouseButtonPress))
    {
        QMouseEvent* mouseEvent = static_cast<QMouseEvent*>(event);
        QPalette *p = new QPalette;
        p->setColor(QPalette::WindowText, Qt::red);
        ui->lcdNumber->setPalette(*p);
        ui->lcdNumber_2->setPalette(*p);
        ui->lcdNumber->display((int)ui->widget->xAxis->pixelToCoord(mouseEvent->pos().x()));
        ui->lcdNumber_2->display((int)ui->widget->yAxis->pixelToCoord(mouseEvent->pos().y()));
        if ((mouseEvent->button() == Qt::RightButton) && (!cleardata))
        {
            ui->spinBox->setValue((int)ui->widget->xAxis->pixelToCoord(mouseEvent->pos().x()));
            on_pushButton_clicked();
        }
        if (mouseEvent->button() == Qt::LeftButton)
        {           
            ui->widget->setFocus();
        }
        if ((mouseEvent->button() == Qt::MidButton) && (appcn->ready))
        {
            if ( ((phasepr) && (hronar)) || ((phasepr) && (zerocr)) || ((hronar) && (zerocr)) )
            {
                QMessageBox msgBox;
                msgBox.setText("Only one method for stimulation should be chosen!");
                msgBox.exec();
            }
            else
            if ((addmode) && (addmodeon))
                addmodeon=false;
            else
            if ((addmode) && (!addmodeon))
            {                
                stims=0; recparts=0;
                startpos=(int)ui->widget->xAxis->pixelToCoord(mouseEvent->pos().x());                
                addmodeon=true;
                offlinedata=false;
                tim->start();
            }
            else           
            {
                startpos=(int)ui->widget->xAxis->pixelToCoord(mouseEvent->pos().x());                
                hnt->imlength=ui->spinBox_5->value();
                hnt->imstprop=ui->doubleSpinBox_3->value();
                hnt->stlength=hnt->imlength*hnt->imstprop;
                stims=0; counter=0;
                zerophaseinit(lcutoff, hcutoff, butterord, hnt->srfr);
                clearstim(startpos,(2*hnt->imlength+hnt->stlength)*hnt->numst,1);
                cleareeg(startpos,startpos+(2*hnt->imlength+hnt->stlength)*hnt->numst);
               // ui->widget->graph(0)->setData(arrc.xc, arrc.amp1);
                ui->widget->graph(3)->setData(arrc.xc, arrc.of1);
                ui->widget->replot();
                offlinedata=false;
                mw->printdata("Sequence start point: "+QString::number(startpos));
                if (randomstims)
                    setrandomrelations();
                else
                    setfixedrelations();
                tim->start();
            }
        }        
    }
    if ((target == ui->widget) && (event->type() == QEvent::KeyPress))
    {
        QKeyEvent* keyEvent = static_cast<QKeyEvent*>(event);
        if (keyEvent->key() == Qt::Key_Left)
        { double low = ui->widget->xAxis->range().lower;
            double high = ui->widget->xAxis->range().upper;
            ui->widget->xAxis->moveRange(-((high - low) / stepsPerPress));
            ui->widget->replot();
        }
        if (keyEvent->key() == Qt::Key_Right)
        {double low = ui->widget->xAxis->range().lower;
            double high = ui->widget->xAxis->range().upper;
            ui->widget->xAxis->moveRange((high - low) / stepsPerPress);
            ui->widget->replot();
        }
        if (keyEvent->key()==Qt::Key_Control)
        {
            ui->widget->axisRect()->setRangeZoom(Qt::Horizontal);
        }
        if (keyEvent->key()==Qt::Key_Alt)
        {
            ui->widget->axisRect()->setRangeZoom(Qt::Vertical);
        }


    }
    return false;
}

void plotwindow::setrandomrelations()
{
    for (int i=totalstims; i<totalstims+hnt->numst; i++)
        if (i<totalstims+hnt->numst/2)
            relations[i]=1;
        else
            relations[i]=-1;
    random_shuffle(&relations[totalstims], &relations[totalstims+hnt->numst]);
    //for (int i=0; i<hnt->numst; i++)
    //    cout << relations[totalstims+i] << " ";
    //cout << endl;
}

void plotwindow::setfixedrelations()
{
    for (int i=totalstims; i<totalstims+hnt->numst; i++)
        if (ui->radioButton->isChecked())
            relations[i]=1;
        else
            relations[i]=-1;
}

void plotwindow::clean()
{
    for (int i=0; i<NMAX; i++)
    {
        arrc.amp1[i]=0;
        arrc.of1[i]=0;
        arrc.of3[i]=0;
        arrc.of4[i]=0;
        arrc.of5[i]=0;
        arrc.of6[i]=0;
    }
    ui->widget->replot();
}

void plotwindow::printtoresultstring(QString str)
{
    strLst1.push_front(str);
    strLstM1->setStringList(strLst1);
}

void plotwindow::printtoresultbox(QString str)
{
    strLst2.push_front(str);
    strLstM2->setStringList(strLst2);
}

void plotwindow::printtoresultsync(QString str)
{
    strLst3.push_front(str);
    strLstM3->setStringList(strLst3);
}

void plotwindow::doplot()
{           
    ui->widget->setFixedSize(1470,700);
    ui->listView->setGeometry(60,780,860,20);
    ui->listView->setModel(strLstM1);
    ui->listView->show();
    ui->listView->setAutoScroll(true);
    ui->listView_2->setGeometry(1175,753,145,20);
    ui->listView_2->setModel(strLstM2);
    ui->listView_2->show();
    ui->listView_2->setAutoScroll(true);
    ui->listView_3->setGeometry(1330,753,145,20);
    ui->listView_3->setModel(strLstM3);
    ui->listView_3->show();
    ui->listView_3->setAutoScroll(true);
    ui->label->setGeometry(300,710,140,25);
    ui->label_2->setGeometry(520,710,140,25);
    ui->label_3->setGeometry(730,710,140,25);
    ui->label_4->setVisible(false);// 905,710,140,25
    ui->label_5->setGeometry(800,750,140,25);
    ui->label_6->setGeometry(100,750,140,25);
    ui->label_7->setGeometry(665,750,140,25);
    ui->label_8->setGeometry(905,750,140,25);//
    ui->label_8->setVisible(false);
    ui->label_12->setGeometry(510,750,80,25);
    ui->label_13->setGeometry(300,750,130,25);
    ui->spinBox->setGeometry(440,710,60,25);
    ui->spinBox_2->setGeometry(640,710,50,25);
    ui->spinBox_3->setVisible(false);// 815,710,50,25
    ui->spinBox_4->setGeometry(840,750,50,25);
    ui->spinBox_5->setGeometry(215,750,50,25);
    ui->spinBox_6->setGeometry(730,750,50,25);
    ui->spinBox_7->setGeometry(595,750,50,25);
    ui->doubleSpinBox->setGeometry(815,710,50,25);
    ui->doubleSpinBox_2->setGeometry(990,750,50,25);
    ui->doubleSpinBox_3->setGeometry(440,750,50,25);

    ui->checkBox->setGeometry(920,710,100,25);
    ui->checkBox_2->setGeometry(920,732,100,25);
    ui->radioButton->setGeometry(925,755,90,25);
    ui->radioButton_4->setGeometry(925,775,110,25);
    ui->radioButton_2->setGeometry(1010,755,75,25);
    ui->radioButton_5->setGeometry(1090,755,70,25);
    ui->radioButton_3->setGeometry(1050,775,110,25);
    ui->checkBox_3->setGeometry(1040,732,100,25);
    ui->checkBox_4->setGeometry(1040,710,125,25);

    ui->pushButton->setGeometry(1270,780,70,25);
    ui->pushButton_6->setGeometry(1360,780,70,25);
    ui->pushButton_2->setGeometry(20,710,35,30);
    ui->pushButton_3->setGeometry(55,710,35,30);
    ui->pushButton_5->setGeometry(18,745,74,30);
    ui->pushButton_4->setGeometry(1180,780,70,25);
    ui->pushButton_5->setVisible(true);
    ui->widget->setFocus();
    ui->spinBox->setValue(hnt->posstim);
    ui->spinBox_2->setValue(hnt->lfilt);
    ui->spinBox_3->setValue(hnt->phase);
    ui->spinBox_4->setValue(hnt->noise);
    ui->spinBox_5->setValue(hnt->imlength);
    ui->spinBox_6->setValue(hnt->stampl);
    ui->spinBox_7->setValue(hnt->numst);
    ui->doubleSpinBox->setValue(hnt->osfr);    
    ui->doubleSpinBox_2->setValue(hnt->thfr);
    ui->doubleSpinBox_2->setVisible(false);
    ui->doubleSpinBox_3->setValue(hnt->imstprop);
    ui->checkBox_2->setChecked(hnt->twooscil);
    ui->checkBox_3->setChecked(usefiltering);
    ui->lcdNumber->setGeometry(100,35,70,25);
    ui->lcdNumber_2->setGeometry(190,35,70,25);
    QPalette *p = new QPalette(Qt::green);      
    ui->lcdNumber_4->setPalette(*p);    
    ui->spinBox_8->setGeometry(215,710,50,25);
    ui->lcdNumber_4->setGeometry(210,70,50,25);  
    ui->label_11->setGeometry(1175,710,150,40);
    ui->label_14->setGeometry(1330,710,150,40);
    ui->label_9->setGeometry(100,710,140,25);
    ui->label_10->setGeometry(115,70,140,25);
    ui->spinBox_8->setValue(hnt->srfr);
    ui->lcdNumber_4->display(hnt->srfr/hnt->osfr);   
    stepsPerPress = 10;
    startpos=hnt->posstim;
    plot(ui->widget);   
}

void plotwindow::plot(QCustomPlot *customPlot)
{    
    if (start) return;
    sw = new Settings();
    QCPItemText *textLabel1 = new QCPItemText(customPlot);
    customPlot->addItem(textLabel1);
    textLabel1->setPositionAlignment(Qt::AlignTop|Qt::AlignHCenter);
    textLabel1->position->setType(QCPItemPosition::ptAxisRectRatio);
    textLabel1->position->setCoords(0.5, 0); // place position at center/top of axis rect
    textLabel1->setText(" Mouse MiddleButton - online mode with BCI2000 signal to tACS output \n \
MiddleButton in 'Additional mode' - start/stop EEG recording or stimulation \n \
RightButton / RUN Button - offline mode with preloaded data");
    textLabel1->setFont(QFont(font().family(), 10));
    textLabel1->setPen(QPen(Qt::gray));
    customPlot->addGraph();
    customPlot->legend->setVisible(true);
    QFont legendFont = font();
    legendFont.setPointSize(9);
    legendFont.setBold(true);
    customPlot->legend->setFont(legendFont);
    customPlot->legend->setBrush(QBrush(QColor(255,255,255,230)));    
  //  QPen pen1;
  //  pen1.setWidth(1);
  //  pen1.setColor(QColor(Qt::green));
  //  customPlot->graph(0)->setPen(pen1);
    customPlot->graph(0)->setPen(QPen(Qt::green));
    customPlot->graph(0)->setData(arrc.xc, arrc.amp1);
  //  customPlot->graph(0)->setData(datamap_ptr,false);
    customPlot->graph(0)->setName("EEG signal ");
    customPlot->graph(0)->setAdaptiveSampling(true);
    customPlot->addGraph();
  //  QPen pen2;
  //  pen2.setWidth(1);
  //  pen2.setColor(QColor(Qt::darkGreen));
  //  customPlot->graph(1)->setPen(pen2);
    customPlot->graph(1)->setPen(QPen(Qt::darkGreen));
    customPlot->graph(1)->setData(arrc.xc, arrc.amp2);
    customPlot->graph(1)->setName("Filtered data");
    customPlot->addGraph();
  //  QPen pen3;
  //  pen3.setWidth(1);
  //  pen3.setColor(QColor("#ffa500")); // "#ffa500" - orange
  //  customPlot->graph(2)->setPen(pen3);
    customPlot->graph(2)->setData(arrc.xc, arrc.of2);
    customPlot->graph(2)->setPen(QPen("#ffa500"));
    customPlot->graph(2)->setName("Autoregression data (AR)");
    customPlot->addGraph();
  //  QPen pen4;
  //  pen4.setWidth(2.2);
  //  pen4.setColor(QColor(Qt::darkMagenta));
  //  customPlot->graph(3)->setPen(pen4);
    customPlot->graph(3)->setPen(QPen(Qt::darkMagenta));
    customPlot->graph(3)->setData(arrc.xc, arrc.of1);
    customPlot->graph(3)->setName("Stimulation based on Hilbert transform (HT)");
//    customPlot->graph(3)->setPen(Qt::DotLine);
    customPlot->addGraph();
  //  QPen pen5;
  //  pen5.setWidth(2.2);
  //  pen5.setColor(QColor(Qt::darkBlue));
  //  customPlot->graph(4)->setPen(pen5);
    customPlot->graph(4)->setPen(QPen(Qt::darkBlue));
    customPlot->graph(4)->setData(arrc.xc, arrc.of3);
    customPlot->graph(4)->setName("Stimulation based on Autoregression (AR)");
//    customPlot->graph(4)->setPen(Qt::DashLine);
    customPlot->addGraph();
  //  QPen pen6;
  //  pen6.setWidth(2.2);
  //  pen6.setColor(QColor(Qt::red));
  //  customPlot->graph(5)->setPen(pen6);
    customPlot->graph(5)->setPen(QPen(Qt::red));
    customPlot->graph(5)->setData(arrc.xc, arrc.of4);
    customPlot->graph(5)->setName("Stimulation based on Zero Crossing (ZC)");
//    customPlot->graph(5)->setPen(Qt::DashDotLine);
    customPlot->addGraph();
  //  QPen pen7;
  //  pen7.setWidth(2.2);
  //  pen7.setColor(QColor(Qt::black));
  //  customPlot->graph(6)->setPen(pen7);
    customPlot->graph(6)->setPen(QPen(Qt::black));
    customPlot->graph(6)->setData(arrc.xc, arrc.of6);
    customPlot->graph(6)->setName("Stimulation based on FFT values (FFT)");
    customPlot->addGraph();
  //  QPen pen8;
  //  pen8.setWidth(3);
  //  pen8.setColor(QColor("#ffa500"));
  //  customPlot->graph(7)->setPen(pen8);
    customPlot->graph(7)->setPen(QPen("#ffa500"));
    customPlot->graph(7)->setData(arrc.xc, arrc.of5);
    customPlot->graph(7)->setName("Optimal stimulation (offline mode)");
    customPlot->xAxis->setLabel("x");
    customPlot->yAxis->setLabel("y");    
    customPlot->rescaleAxes();
    customPlot->xAxis->setRange(0, hnt->srfr * 6);
    customPlot->yAxis->setRange(-150, 150);
    customPlot->xAxis->moveRange(hnt->posstim - hnt->srfr / 2);
    customPlot->setInteraction(QCP::iRangeDrag, true);
    customPlot->setInteraction(QCP::iRangeZoom, true);
    customPlot->axisRect()->setRangeZoom(Qt::Horizontal);   
    customPlot->replot();
    start=true;
    ssdtim->start();
}

void plotwindow::updatedata()
{
    ui->widget->graph(0)->setData(arrc.xc, arrc.amp1);
    ui->widget->replot();
}

double plotwindow::phasediff(double* stim1, double* stim2, int length)
{
    double tp1[length];
    double tp2[length];
    plv=0;

    hnt->init();
    hnt->npt=length;

    double avg=0;
    for (int i = 0; i<hnt->npt; i++)
        avg+=stim1[i];
    avg/=length;

    for (int i = 0; i<hnt->npt; i++)
        hnt->x[i]=stim1[i]-avg;
    hnt->firhilbert();
   // complex<double>* wavel = morletwavelet(hnt->x,0,hnt->npt);
    for (int i = 0; i<hnt->npt; i++)
        tp1[i]=hnt->phase1[i];
      //  tp1[i] = atan2(wavel[i].imag(),wavel[i].real());



  //  for (int i = 0; i<hnt->npt; i++)
  //      cout<<(int)qRadiansToDegrees(tp1[i])<<" ";
  //  cout<<endl;

    hnt->init();
    hnt->npt=length;

    avg=0;
    for (int i = 0; i<hnt->npt; i++)
        avg+=stim2[i];
    avg/=length;

    for (int i = 0; i<hnt->npt; i++)
        hnt->x[i]=stim2[i]-avg;
    //complex<double>* wavel = morletwavelet(hnt->x,0,hnt->npt);
    hnt->firhilbert();
    for (int i = 0; i<hnt->npt; i++)
        tp2[i]=hnt->phase1[i];
      //  tp2[i] = atan2(wavel[i].imag(),wavel[i].real());


  //  for (int i = 0; i<hnt->npt; i++)
  //      cout<<(int)qRadiansToDegrees(tp2[i])<<" ";
  //  cout<<endl;

    double diff=0;
    complex<double> pv(0.0,0.0);    
    for (int i = hnt->lfilt/2+1; i<hnt->npt-hnt->lfilt/2+1; i++)
    {
        diff+=pow(tp1[i]-tp2[i],2);

        pv+=exp(I*(tp1[i]-tp2[i]));

        hnt->fillphasebins(tp1[i]-tp2[i]);
    } 

    plv=fabs(pv)/(length-hnt->lfilt);
   // qDebug()<<plv<<" ";
    diff=sqrt(diff);
   // qDebug()<<diff;

    return diff;
}

double plotwindow::arphasediff(int posstim, int length)
{
    double st1[length];
    for (int i=0; i<length; i++)
        st1[i]=arrc.amp1[posstim+i];
    double st2[length];
    for (int i=0; i<length; i++)
        st2[i]=arrc.of3[posstim+i];
    hnt->phsdif = phasediff(st1, st2, length);
    return hnt->phsdif;
}

void plotwindow::clearstim(int posstim, int length, int method)
{
    if (hronarrun)
    {
        for (int i = posstim; i < posstim+length; i++)
            arrc.of3[i]=0;
    }
    else
    for (int i = posstim; i < posstim+length; i++)
        if (method==1)
            arrc.of1[i]=0;
        else if (method==2)
            arrc.of3[i]=0;
        else if (method==3)
            arrc.of4[i]=0;
    for (int i = posstim; i < posstim+length; i++)
        arrc.of5[i]=0;
}

void plotwindow::generatesignals(int pos, int length, double fr1, double fr2, double fr3, double fr4, int phase, double noise)
{
    // do stim seq on clear signals, save opt phase shift
    // do stim seq on set of noisy signals, calculate DD to clear opt phase shifts
    for (int i=pos-hnt->imlength; i<pos+length-hnt->imlength; i++)
    {
        arrc.amp1[i]=0;
        arrc.of1[i]=0;
        arrc.of3[i]=0;
        arrc.of4[i]=0;
        arrc.of5[i]=0;
        arrc.of6[i]=0;
    }    

    double amp = ampval;
    double nois, pn, ps, snr;
    pn = 0; ps = 0; snr=0;
    hnt->stampl=10;

    const double mean = 0.0;
    const double stddev = 0.5;
    std::default_random_engine generator;
    std::normal_distribution<double> dist(mean, stddev);

    for (int i = pos-hnt->imlength; i < pos+length-hnt->imlength; i++)
    {
        nois = noise*dist(generator);
     //   arrc.amp1[i] =  hnt->stampl * 1.5 * sin (fr2*2*PI/hnt->srfr*(i+phase)) + nois;
        arrc.amp1[i] =  hnt->stampl / 4 * sin (fr1*2*PI/hnt->srfr*(i+phase)) + hnt->stampl * 1.5 * sin (fr2*2*PI/hnt->srfr*(i+phase)) + hnt->stampl / 2 * sin (fr3*2*PI/hnt->srfr*(i+phase)) + hnt->stampl / 4 * sin (fr4*2*PI/hnt->srfr*(i+phase)) + nois;
        pn += pow(nois,2);
        ps += pow(arrc.amp1[i]-nois,2);
    }

    pn=pow(pn/length,0.5);
    ps=pow(ps/length,0.5);
    if (pn>0)
        snr = 10*log10(ps/pn);
   // qDebug()<<ps<<pn;
    if (snr!=0)
        qDebug()<<"SNR (dB): "<<snr;
    else
        qDebug()<<"= no noise =";
    if (!estimatenoisy)
        estimateclear=true;

    ui->widget->graph(0)->setData(arrc.xc, arrc.amp1);
    ui->widget->replot();
}

void plotwindow::filloptstim(int posstim, int length, int phase, int shift)
{
    double noise;
    for (int i = posstim; i < posstim+length; i++)  // fill stimulation interval
    {
        if (hnt->noise == 0)
            noise = 0;
        else
            noise = (1.5 + hnt->noise ) * (0.5 - ((double)rand()/(double)RAND_MAX));
        noise = 0;

        double amp;
        if (adaptampl)
            amp = arrc.amp1[posstim-1];
        else
            amp = ampval;

        if (extractoptphase)
            arrc.of5[i] = stimshift + shift + amp + amp / tscale * 2 * hnt->stampl * sin (hnt->osfr*2*PI/hnt->srfr*(i+phase)) + noise;
        else if (hnt->twooscil)
            arrc.of5[i] = stimshift + shift + amp + amp / tscale * hnt->stampl / 2 * sin (hnt->osfr*2*PI/hnt->srfr*(i+phase)) + hnt->stampl / 2 * sin (hnt->thfr*2*PI/hnt->srfr*(i+phase)) + noise;
        else
            arrc.of5[i] = stimshift + shift + amp + amp / tscale / 2 * hnt->stampl * sin (hnt->osfr*2*PI/hnt->srfr*(i+phase)) + noise;
    }
}

complex<double>* plotwindow::morletwavelet(float* inp, int poswavelet, int length)
{
  //  timer.restart();
    int nKern=263; //263 for 250, 1549 for 500
    double* time  = new double[nKern]; // best practice is to have time=0 at the center of the wavelet
    for (int i=-nKern/2; i<nKern/2+1; i++)
        time[i+nKern/2]=(double)i/1000;
    double frex  = 10.0; // wavelet frequency

    complex<double>* sine_wave = new complex<double>[nKern];
    for (int i=0; i<nKern; i++)
        sine_wave[i]=exp(I*2.0*M_PI*frex*time[i]);

    //  create Gaussian window
    double s = 7 / (2*M_PI*frex);      // this is the standard deviation of the gaussian
    complex<double>* gaus_win = new complex<double>[nKern];
    for (int i=0; i<nKern; i++)
        gaus_win[i]=exp(-pow(time[i],2)/(2*pow(s,2)));

    // now create Morlet wavelet
    complex<double>* cmw = new complex<double>[nKern];
    for (int i=0; i<nKern; i++)
        cmw[i] = sine_wave[i] * gaus_win[i];

  //  for (int i=0; i<1001; i++)
  //      arrc.amp1[1000+i]=cmw[i].real()*30;
  //  ui->widget->graph(0)->setData(arrc.xc, arrc.amp1);
  //  ui->widget->replot();

    int nData = length;
    int nConv = nData + nKern - 1;

    // zero-padding for Morlet wavelet
    complex<double>* cmwfft = new complex<double>[nConv];
    for (int i=0; i<nConv; i++)
    {
        if (i<nKern)
            cmwfft[i]=cmw[i];
        else
            cmwfft[i]=0;
    }
    // zero-padding for data
    complex<double>* datafft = new complex<double>[nConv];
    for (int i=0; i<nConv; i++)
    {
        if (i<length)
            datafft[i].real(inp[poswavelet+i]);
        else
            datafft[i].real(0);
        //datafft[i].imag(0);
    }

    CArray mfft = CArray(cmwfft,nConv);
    hnt->fft(mfft);

    // normalize fft of wavelet
   /* double maxampl=0;
    for (int i=0; i<nConv; i++)
        if (mfft[i].real()>maxampl)
            maxampl=sqrt(mfft[i].real()*mfft[i].real()+mfft[i].imag()*mfft[i].imag());
    for (int i=0; i<nConv; i++)
        mfft[i]/=maxampl; */

    CArray dfft = CArray(datafft,nConv);
    hnt->fft(dfft);

    // now convolution
    complex<double>* conv = new complex<double>[nConv];
    for (int i=0; i<nConv; i++)
        conv[i]=mfft[i]*dfft[i];

    CArray conv_time = CArray(conv,nConv);
    hnt->ifft(conv_time);
    complex<double>* conv_res = new complex<double>[nConv-nKern+1];
    for (int i=0; i<nConv-nKern+1; i++)
        conv_res[i]=conv_time[i+length];

  //  cout << timer.nsecsElapsed() / 1000000.0 << endl;

    for (int i=0; i<length; i++)
    {
      //  arrc.amp2[poswavelet+i]=conv_res[i].real()/40;
     //   cout<<(int)qRadiansToDegrees(atan2(conv_res[i].imag(),conv_res[i].real()))<<" ";
    }
   // cout<<endl;


   // timer.restart();

  /*  hnt->init();
    hnt->npt=length;
    for (int i = 0; i<length; i++)
        hnt->x[i]=arrc.amp1[poswavelet+i];
    hnt->firhilbert();

  //  cout << timer.nsecsElapsed() / 1000000.0 << endl;

    for (int i = 0; i<length; i++)
        cout<<(int)qRadiansToDegrees(hnt->phase1[i])<<" ";
    cout<<endl; */

   // ui->widget->graph(1)->setData(arrc.xc, arrc.amp2);
  //  ui->widget->replot();

    return conv_res;
}

void plotwindow::fftphasepred(int posstim, int length)
{
    zerophaseinit(7, 11, 1, hnt->srfr);
    //hnt->imlength=hnt->imlength*hnt->extinterval;
    //if ((!zerocr))// || (estimatenoisy))
    zerophasefiltzc(posstim-length,length);

    int lt=512;
    Complex t[lt];
    for (int i=posstim-length; i<posstim; i++)
        t[i-posstim+length].real(arrc.amp2[i]);

    for (int i=0; i<lt; i++)
        t[i].imag(0);

    for (int i=0; i<lt-length; i++)  // zero-padding
        t[length+i].real(0);

    cdata = CArray(t,lt);
    hnt->fft(cdata);
    double* frampl = new double[length];
    double* phs = new double[length];
    for (int i=6; i<15; i++)
    {
        frampl[i]=2*sqrt(cdata[i].real()*cdata[i].real()+cdata[i].imag()*cdata[i].imag());
        frampl[i]/=lt;
        phs[i]=(int)qRadiansToDegrees(atan2(cdata[i].imag(),cdata[i].real()));
    //    cout<<frampl[i]<<" ";
    }
   // cout<<endl;

    int domfreq=0; double maxamp=0;
    for (int i=7; i<14; i++) //
        if (frampl[i]>maxamp)
        {
            maxamp=frampl[i];
            domfreq=i;
        }
  //  cout<<"Dominant freq: "<<domfreq<<" Hz with phase: "<<phs[domfreq]<<endl;
   // fftphase = phs[domfreq-1]; //domfreq+1
    fftphase = phs[(int)round(hnt->osfr-1)];
    int sh=length;
    int tp=0;
    if (sh<hnt->srfr/hnt->osfr)
        tp=sh;
     else
        tp=sh % (int)(hnt->srfr/hnt->osfr);
   // cout<<sh<<" "<<tp<<endl;
    hnt->phase = fftphase / ((360 * hnt->osfr) / hnt->srfr) + tp;
    if (!hnt->inphase)
        hnt->phase = hnt->phase + (int)(hnt->srfr/hnt->osfr)/2;
   // int iaf = hnt->osfr;
  //  hnt->osfr = domfreq;
    fillstimforfft(posstim,hnt->stlength);
  //  hnt->osfr = iaf;
   /* for (int i=0; i<15; i++)
        cout<<frampl[i]<<" ";
    cout<<endl;
    for (int i=0; i<15; i++)
        cout<<phs[i]<<" ";
    cout<<endl; */

  //  ui->widget->graph(0)->setData(arrc.xc, arrc.amp1);
//    ui->widget->graph(3)->setData(arrc.xc, arrc.of1);
//    ui->widget->replot();

 //   hnt->imlength=hnt->imlength/hnt->extinterval;
}

void plotwindow::fillstimforfft(int posstim, int length)
{
    double noise;
    for (int i = posstim; i < posstim+length; i++)  // fill stimulation interval
    {
        if (hnt->noise == 0)
            noise = 0;
        else
            noise = (1.5 + hnt->noise ) * (0.5 - ((double)rand()/(double)RAND_MAX));
        noise=0;

        double amp;
        if (adaptampl)
            amp = arrc.amp1[posstim-1];
        else
            amp = ampval;

         // !!! with external funcs for of3 expression - much much much slower !!!
        if (hnt->twooscil)
            if (offlinedata)
                arrc.of6[i] = stimshift + amp + amp / tscale * hnt->stampl / 2 * sin (hnt->osfr*2*PI/hnt->srfr*(i+hnt->phase)) + hnt->stampl / 2 * sin (hnt->thfr*2*PI/hnt->srfr*(i+hnt->phase)) + noise;
            else
                arrc.of6[i] = stimshift + amp / tscale * hnt->stampl * sin (hnt->osfr*2*PI/hnt->srfr*(i+hnt->phase)) + hnt->stampl * 1.5 * sin (hnt->thfr*2*PI/hnt->srfr*(i+hnt->phase)) + noise;
        else
            if (offlinedata)
                arrc.of6[i] = stimshift + amp + amp / tscale / 2 * hnt->stampl * sin (hnt->osfr*2*PI/hnt->srfr*(i+hnt->phase));// + noise;
            else
                arrc.of6[i] = stimshift + amp / tscale * hnt->stampl * 2 * sin (hnt->osfr*2*PI/hnt->srfr*(i+hnt->phase));// + noise;

    }

    double st1[length];
    for (int i=0; i<length; i++)
        st1[i]=arrc.of6[posstim+i];
    double st2[length];
    for (int i=0; i<length; i++)
        st2[i]=arrc.amp1[posstim+i];
    hnt->phsdif = phasediff(st1, st2, length);
}

void plotwindow::fillrandstim(int posstim, int length)
{
    double noise;
    for (int i = posstim; i < posstim+length; i++)
    {
        noise = (0.5 - ((double)rand()/(double)RAND_MAX));
        double amp = ampval;
        arrc.of1[i] =  88*noise;
       // qDebug()<<arrc.of1[i];

    }
}

void plotwindow::fillstim(int posstim, int length)
{
    double noise;
    for (int i = posstim; i < posstim+length; i++)  // fill stimulation interval
    {
        if (hnt->noise == 0)
            noise = 0;
        else
            noise = (1.5 + hnt->noise ) * (0.5 - ((double)rand()/(double)RAND_MAX));        
        noise=0;

        double amp;
        if (adaptampl)
            amp = arrc.amp1[posstim-1];
        else
            amp = ampval;

         // !!! with external funcs for of3 expression - much much much slower !!!
        if (hnt->twooscil)
            if (offlinedata)
                arrc.of1[i] = stimshift + amp + amp / tscale * hnt->stampl / 2 * sin (hnt->osfr*2*PI/hnt->srfr*(i+hnt->phase)) + hnt->stampl / 2 * sin (hnt->thfr*2*PI/hnt->srfr*(i+hnt->phase)) + noise;
            else
                arrc.of1[i] = stimshift + amp / tscale * hnt->stampl * sin (hnt->osfr*2*PI/hnt->srfr*(i+hnt->phase)) + hnt->stampl * 1.5 * sin (hnt->thfr*2*PI/hnt->srfr*(i+hnt->phase)) + noise;
        else
            if (offlinedata)
                arrc.of1[i] = stimshift + amp + amp / tscale / 2 * hnt->stampl * sin (hnt->osfr*2*PI/hnt->srfr*(i+hnt->phase));// + noise;
            else
                arrc.of1[i] = stimshift + amp / tscale * hnt->stampl * 2 * sin (hnt->osfr*2*PI/hnt->srfr*(i+hnt->phase));// + noise;

    }

    double st1[length];
    for (int i=0; i<length; i++)
        st1[i]=arrc.of1[posstim+i];
    double st2[length];
    for (int i=0; i<length; i++)
        st2[i]=arrc.amp1[posstim+i];
    hnt->phsdif = phasediff(st1, st2, length);
}

void plotwindow::fillstimforstim(int posstim, int length)
{
    double noise;
    for (int i = posstim; i < posstim+length; i++)  // fill stimulation interval
    {
        if (hnt->noise == 0)
            noise = 0;
        else
            noise = (1.5 + hnt->noise ) * (0.5 - ((double)rand()/(double)RAND_MAX));
        noise=0;

        double amp;
        if (adaptampl)
            amp = arrc.amp1[posstim-1];
        else
            amp = ampval;

         // !!! with external funcs for of3 expression - much much much slower !!!
        if (hnt->twooscil)
            if (offlinedata)
                arrc.amp2[i] = stimshift + amp + amp / tscale * hnt->stampl / 2 * sin (hnt->osfr*2*PI/hnt->srfr*(i+hnt->phase)) + hnt->stampl / 2 * sin (hnt->thfr*2*PI/hnt->srfr*(i+hnt->phase)) + noise;
            else
                arrc.amp2[i] = stimshift + amp / tscale * hnt->stampl * sin (hnt->osfr*2*PI/hnt->srfr*(i+hnt->phase)) + hnt->stampl * 1.5 * sin (hnt->thfr*2*PI/hnt->srfr*(i+hnt->phase)) + noise;
        else
            if (offlinedata)
                arrc.amp2[i] = stimshift + amp + amp / tscale / 2 * hnt->stampl * sin (hnt->osfr*2*PI/hnt->srfr*(i+hnt->phase));// + noise;
            else
                arrc.amp2[i] = stimshift + amp / tscale * hnt->stampl * 2 * sin (hnt->osfr*2*PI/hnt->srfr*(i+hnt->phase));// + noise;

    }

    double st1[length];
    for (int i=0; i<length; i++)
        st1[i]=arrc.of1[posstim+i];
    double st2[length];
    for (int i=0; i<length; i++)
        st2[i]=arrc.amp2[posstim+i];
    hnt->phsdif = phasediff(st1, st2, length);
}

void plotwindow::fillstimforzc(int posstim, int length, int phaseshift)
{
    double noise;
    for (int i = posstim; i < posstim+length; i++)  // fill stimulation interval
    {
        if (hnt->noise == 0)
            noise = 0;
        else
            noise = (1.5 + hnt->noise ) * (0.5 - ((double)rand()/(double)RAND_MAX));
        noise=0;

        double amp;
        if (adaptampl)
            amp = arrc.amp1[posstim-1];
        else
            amp = ampval;

         // !!! with external funcs for of3 expression - much much much slower !!!
        if (hnt->twooscil)
            if (offlinedata)
                arrc.of4[i] = stimshift + amp + amp / tscale * hnt->stampl / 2 * sin (hnt->osfr*2*PI/hnt->srfr*(i+phaseshift)) + hnt->stampl / 2 * sin (hnt->thfr*2*PI/hnt->srfr*(i+phaseshift)) + noise;
            else
                arrc.of4[i] = stimshift + amp / tscale * hnt->stampl * sin (hnt->osfr*2*PI/hnt->srfr*(i+phaseshift)) + hnt->stampl * 1.5 * sin (hnt->thfr*2*PI/hnt->srfr*(i+phaseshift)) + noise;
        else
            if (offlinedata)
                arrc.of4[i] = stimshift + amp + amp / tscale / 2 * hnt->stampl * sin (hnt->osfr*2*PI/hnt->srfr*(i+phaseshift));// + noise;
            else
                arrc.of4[i] = stimshift + amp / tscale * hnt->stampl * 2 * sin (hnt->osfr*2*PI/hnt->srfr*(i+phaseshift));// + noise;

    }

    double st1[length];
    for (int i=0; i<length; i++)
        st1[i]=arrc.of4[posstim+i];
    double st2[length];
    for (int i=0; i<length; i++)
        st2[i]=arrc.amp1[posstim+i];
    hnt->phsdif = phasediff(st1, st2, length);
}

void plotwindow::fillstimforregr(int posstim, int length)
{
    double noise;    
    for (int i = posstim; i < posstim+length; i++)  // fill stimulation interval
    {
        if (hnt->noise == 0)
            noise = 0;
        else
            noise = (1.5 + hnt->noise ) * (0.5 - ((double)rand()/(double)RAND_MAX));
        noise=0;

        double amp;
        if (adaptampl)
            amp = arrc.of2[posstim+1] ;
        else
            amp = ampval ;

        if (hnt->twooscil)
            if (offlinedata)
                arrc.of3[i] = stimshift + amp + amp / tscale * hnt->stampl / 2 * sin (hnt->osfr*2*PI/hnt->srfr*(i+hnt->phase)) + hnt->stampl / 2 * sin (hnt->thfr*2*PI/hnt->srfr*(i+hnt->phase)) + noise;
            else
                arrc.of3[i] = stimshift + amp / tscale * hnt->stampl * sin (hnt->osfr*2*PI/hnt->srfr*(i+hnt->phase)) + hnt->stampl * 1.5 * sin (hnt->thfr*2*PI/hnt->srfr*(i+hnt->phase)) + noise;
        else
            if (offlinedata)
                arrc.of3[i] = stimshift + amp + amp / tscale / 2 * hnt->stampl * sin (hnt->osfr*2*PI/hnt->srfr*(i+hnt->phase)) + noise;
            else
                arrc.of3[i] = stimshift + amp / tscale * hnt->stampl * 2 * sin (hnt->osfr*2*PI/hnt->srfr*(i+hnt->phase)) + noise;      

    }

    double st1[length];
    for (int i=0; i<length; i++)
        st1[i]=arrc.of3[posstim+i];
    double st2[length];
    if (hronarrun)
    {
        for (int i=0; i<length; i++)
            st2[i]=arrc.of2[posstim+i];
    } else if (offlinedata)
    {
        for (int i=0; i<length; i++)
            st2[i]=arrc.amp1[posstim+i];
    }
    hnt->phsdif = phasediff(st1, st2, length);
}

void plotwindow::minphasediff(int posstim, int length)
{   
    double t;    
    hnt->diff=100000;
    int optd=0;
    // qDebug()<<hnt->mind<<" "<<hnt->maxd;
    t=hnt->mind;
    while (t<=hnt->maxd)
    {
        hnt->phase=t;
        if (extractoptphase)
            fillstimforstim(posstim,length);
        else if (!hronarrun)
            fillstim(posstim,length);
        else
            fillstimforregr(posstim,length);
        if (hnt->phsdif<hnt->diff)
        {
            hnt->diff=hnt->phsdif;
            optd=t;
        }
        t+=1;
    }
    hnt->phase=optd;
    if (hnt->twooscil)
    {
        hnt->diff=100000;
        t=hnt->mintheta;
        double opttheta=0.0;
        while (t<=hnt->maxtheta)
        {
            hnt->thfr=t;
            if (extractoptphase)
                fillstimforstim(posstim,length);
            else if (!hronarrun)
                fillstim(posstim,length);
            else
                fillstimforregr(posstim,length);
            if (hnt->phsdif<hnt->diff)
            {
                hnt->diff=hnt->phsdif;
                opttheta=t;
            }
            t+=0.1;
        }
        hnt->thfr=opttheta;
    }
}

void plotwindow::maxphasediff(int posstim, int length)
{
    double t;
    hnt->diff=0.0001;
    int optd=0;
    t=hnt->mind;
    while (t<=hnt->maxd)
    {
        hnt->phase=t;
        if (extractoptphase)
            fillstimforstim(posstim,length);
        else if (!hronarrun)
            fillstim(posstim,length);
        else
            fillstimforregr(posstim,length);
        if (hnt->phsdif>hnt->diff)
        {
            hnt->diff=hnt->phsdif;
            optd=t;
        }
        t+=1;
    }
    hnt->phase=optd;
    if (hnt->twooscil)
    {
        hnt->diff=0.0001;
        t=hnt->mintheta;
        double opttheta=0.0;
        while (t<=hnt->maxtheta)
        {
            hnt->thfr=t;
            if (extractoptphase)
                fillstimforstim(posstim,length);
            else if (!hronarrun)
                fillstim(posstim,length);
            else
                fillstimforregr(posstim,length);
            if (hnt->phsdif>hnt->diff)
            {
                hnt->diff=hnt->phsdif;
                opttheta=t;
            }
            t+=0.1;
        }
        hnt->thfr=opttheta;
    }
}

void plotwindow::fordemoshow()
{
    Sleep(500);
    ui->widget->graph(3)->setData(arrc.xc, arrc.of1);
    ui->widget->replot();
    ui->widget->xAxis->moveRange(3*hnt->imlength/2);
}

double plotwindow::offlinevaluation(int pos, int i, int method)
// find optimal phase and phasediff for stimulation interval
{
    double k=0; double t=0; double diff;
    // qDebug()<<"phsdif prd: "<<hnt->phsdif<<" phase prd: "<<hnt->phase;
    diff = hnt->phsdif; k = hnt->phase; t = hnt->thfr;
    if (hnt->inphase)
    {
         if (estimatenoisy)
         {
              hnt->phase=clearoptphase[i];
              fillstim(pos+(hnt->imlength+hnt->stlength)*i,hnt->stlength);
              hnt->diff=hnt->phsdif;
         } else
              minphasediff(pos+(hnt->imlength+hnt->stlength)*i,hnt->stlength);
        hnt->optphasediff=hnt->diff;
      //  qDebug()<<hnt->optphasediff<<" "<<diff;
       //   qDebug()<<"phsdif opt: "<<hnt->diff<<" phase opt: "<<hnt->phase;;
        diff=hnt->diff/diff;
    }
    else
    {
        if (estimatenoisy)
        {
            hnt->phase=clearoptphase[i];
            fillstim(pos+(hnt->imlength+hnt->stlength)*i,hnt->stlength);
            hnt->diff=hnt->phsdif;
        } else
            maxphasediff(pos+(hnt->imlength+hnt->stlength)*i,hnt->stlength);
        hnt->optphasediff=hnt->diff;
       // qDebug()<<hnt->optphasediff<<" "<<diff;
        //   qDebug()<<"phsdif opt: "<<hnt->diff<<" phase opt: "<<hnt->phase;
        diff/=hnt->diff;
    } 
    clearstim(pos+(hnt->imlength+hnt->stlength)*i,hnt->stlength,1);
    hnt->optphase=hnt->phase;
    hnt->phsdiffdeg = (360 / (hnt->srfr/hnt->osfr))*abs(k-hnt->optphase);
    if (hnt->phsdiffdeg>180)
        hnt->phsdiffdeg=abs(360-hnt->phsdiffdeg);

  //  qDebug()<<"predic acc% "<<diff<<"phase diff deg"<<hnt->phsdiffdeg;
    hnt->phase = k; hnt->thfr = t;// int fl = hnt->lfilt;
    if (method==1)
    {
      //  hnt->lfilt=fl*2; // for more correct plv and synchrony estimation
        fillstim(pos+(hnt->imlength+hnt->stlength)*i,hnt->stlength);
      //  hnt->lfilt=fl;
       // for (int i = 0; i<hnt->binnums; i++)
       //     cout<<hnt->phasebins[i]<<" ";
       // cout<<endl;
    }
    return diff;
}

double plotwindow::zcdifference(int pos)
{
    int length=hnt->stlength;
    //double k=hnt->phsdif;
    double st1[length];
    double st2[length];
    for (int i=0; i<length; i++)
        st1[i]=arrc.amp1[pos+i];
    for (int i=0; i<length; i++)
        st2[i]=arrc.of4[pos+i];
    return phasediff(st1, st2, length);
  //  qDebug()<<"Phase diff HR: "<<k/length<<" Phase diff AR: "<<phsdAR/length;
  // return 0;
}

double plotwindow::gethindex()
{
    double hindex=0;
    double pbin[hnt->binnums];
    for (int i = 0; i<hnt->binnums; i++)
        pbin[i]=(double)hnt->phasebins[i]/(hnt->stlength-hnt->lfilt);
    for (int i = 0; i<hnt->binnums; i++)
        if (pbin[i]>0)
            hindex-=pbin[i]*log(pbin[i]);
    double synch=(log(hnt->binnums)-hindex)/log(hnt->binnums);
    return synch;
}

void plotwindow::strategy(int pos) // MAIN PROCEDURE
{
    if (offlinedata)
        mw->printdata("Sequence start point: "+QString::number(hnt->posstim));
    double procpred=0; double procregr=0; double proczc = 0; double procfft=0;
    double diffHR=0; double diffZC=0; double diffHRonAR=0; double diffFFT=0;
    plvhr=0; plvzc=0; plvHRonAR=0; plvfft=0;
    double currtime = 0;
    if (offlinedata)
    {
        hnt->totalprocpred=0;
        hnt->totalprocregr=0;
        hnt->totalproczc=0;
        hnt->totalprocfft=0;
        hnt->totaldiffHR=0;
        hnt->totaldiffZC=0;
        hnt->totaldiffHRonAR=0;
        hnt->totaldiffFFT=0;
        hnt->averagedegdiffHR=0;
        hnt->averagedegdiffAR=0;
        hnt->averagedegdiffZC=0;
        hnt->averagedegdiffFFT=0;
        hnt->optphasehr=0;
        hnt->optphasear=0;
        hnt->optphasezc=0;
        hnt->optphasefft=0;
        hnt->optphase=0;
        hnt->degdevHR=0; hnt->degdevAR=0; hnt->degdevZC=0; hnt->degdevFFT=0;
        hnt->totalsynchr=0; hnt->totalsyncar=0; hnt->totalsynczc=0; hnt->totalsyncfft=0;
        avgdevhr=0;
    }
    int nms=1;
    if (offlinedata)
    {
        averopttime=0;
        nms=hnt->numst;    
        for (int i=0; i<hnt->numst; i++)
        {
            arrphasediffhr[i]=0;
            arrphasediffar[i]=0;
            arrphasediffzc[i]=0;
            arrphasedifffft[i]=0;
            arrsynchr[i]=0;
            arrsyncar[i]=0;
            arrsynczc[i]=0;
            arrsyncfft[i]=0;
        }
    }
    for (int i=0; i<nms; i++)
    {   // cycle works for offline mode, in online mode - function is called every time for new eeg segment
        // find optimal parameters on imaging interval, fill stimulation interval

        if ((offlinedata) && (randomstims))
        {
            if (relations[i]==1)
                hnt->inphase=true;
            else
                hnt->inphase=false;
        }

        hnt->optphase=0;
        timer.restart();

        if (usefiltering) {
            if ((filterstims) && (offlinedata))
                zerophasefilt(pos+(hnt->imlength+hnt->stlength)*i,hnt->stlength);
               //filterdata(pos+(hnt->imlength+hnt->stlength)*i,hnt->stlength); // for plv test
            if (filtfir)
                filterdata(pos+(hnt->imlength+hnt->stlength)*i-hnt->imlength,hnt->imlength);
            else if (recurbutter)
                recurbutterfilter(pos+(hnt->imlength+hnt->stlength)*i-hnt->imlength,hnt->imlength);
            else if (zerobutter)
                zerophasefilt(pos+(hnt->imlength+hnt->stlength)*i-hnt->imlength,hnt->imlength);
            }
        if (zerocr) // 3rd method - phase prediction based on zero crossing
        {
            zerophasefiltzc(pos+(hnt->imlength+hnt->stlength)*i-hnt->imlength,hnt->imlength);
            zerocrossingstim(pos+(hnt->imlength+hnt->stlength)*i,hnt->imlength);
            if (offlinedata)
                arrsynczc[i]=gethindex();
            hnt->optphasezc=hnt->phase;
            int tp = hnt->srfr/hnt->osfr;
            if (hnt->optphasezc > tp)
                hnt->optphasezc = hnt->optphasezc % tp;
        }
        if ((phasepr) && (!randnoisestim) && (!randphasestim)) // 1st method - phase prediction based on imaging interval
        {                        
            if (ui->checkBox->isChecked()) // alpha frequency estimation - does not work nice !
            {
                hnt->averagefreq();
                ui->doubleSpinBox->setValue(hnt->avg1);
                hnt->osfr=hnt->avg1;
                hnt->mind=(int)(-hnt->srfr/hnt->osfr/2);
                hnt->maxd=(int)(hnt->srfr/hnt->osfr/2);
            }
            if (hnt->inphase) // find optimal phase value by iterative search
                minphasediff(pos+(hnt->imlength+hnt->stlength)*i-hnt->imlength*hnt->extinterval,hnt->imlength*hnt->extinterval);
            else
                maxphasediff(pos+(hnt->imlength+hnt->stlength)*i-hnt->imlength*hnt->extinterval,hnt->imlength*hnt->extinterval);

            hnt->optphasehr = hnt->phase;

          //  qDebug()<<hnt->diff;
            clearstim(pos+(hnt->imlength+hnt->stlength)*i-hnt->imlength*hnt->extinterval,hnt->imlength*hnt->extinterval,1); // clear stim after last iteration
            if (showprediction)
                fillstim(pos+(hnt->imlength+hnt->stlength)*i-hnt->imlength*hnt->extinterval,hnt->imlength*hnt->extinterval); // plot optimal stim on imaging interval
            if ((!offlinedata) && (correctdelay))
                hnt->phase+=(timer.nsecsElapsed() / 1000000.0) / 2 + transdelay + ssdt; // optimization delay for online mode
           // if ((offlinedata) && (compavgdeviation))
           //     hnt->phase+=avgdevhr;

            fillstim(pos+(hnt->imlength+hnt->stlength)*i,hnt->stlength);

            if ((offlinedata) && (estimpredstim))
            {
                arrsynchr[i]=gethindex();
                procpred=offlinevaluation(pos,i,1); // gives hnt->optphasediff value for % of relative accuracy
                hnt->averagedegdiffHR+=hnt->phsdiffdeg;
                arrphasediffhr[i]=hnt->phsdiffdeg;
                plvhr+=plv;
                diffHR=hnt->phsdif;             
              //  procpred=(double)(180-hnt->phsdiffdeg)/180;
            }
        }
        if (hronar) // 2nd method - phase prediction based on autoregression
        {
            runautoreg(pos+(hnt->imlength+hnt->stlength)*i);
            if ((!usefiltering) && (filterardata)) {
                if (filtfir)
                    filterar(pos+(hnt->imlength+hnt->stlength)*i,hnt->stlength);
                else if (recurbutter)
                    recurbutterfilterar(pos+(hnt->imlength+hnt->stlength)*i,hnt->stlength);
                else if (zerobutter)
                    zerophasefiltar(pos+(hnt->imlength+hnt->stlength)*i,hnt->stlength);
            }
            hilbonar(pos+(hnt->imlength+hnt->stlength)*i);
            hnt->optphasear=hnt->phase;
            if ((offlinedata) && (estimpredstim))
            {
                arrsyncar[i]=gethindex();
                diffHRonAR=arphasediff(pos+(hnt->imlength+hnt->stlength)*i,hnt->stlength);
                if (!phasepr)
                {
                    offlinevaluation(pos,i,2);
                    hnt->averagedegdiffAR+=hnt->phsdiffdeg;
                } else
                {                    
                    hnt->phsdiffdeg = (360 / (hnt->srfr/hnt->osfr))*abs(hnt->optphasear-hnt->optphase);
                    if (hnt->phsdiffdeg>180)
                        hnt->phsdiffdeg=abs(360-hnt->phsdiffdeg);
                    hnt->averagedegdiffAR+=hnt->phsdiffdeg;
                }
               // qDebug()<<hnt->optphasediff<<" "<<diffHRonAR;
                arrphasediffar[i]=hnt->phsdiffdeg;
                plvHRonAR+=plv;
                if (hnt->inphase)
                    procregr=hnt->optphasediff/diffHRonAR;
                else if (!hnt->inphase)
                    procregr=diffHRonAR/hnt->optphasediff;
             //   procregr = (double)(180-hnt->phsdiffdeg)/180;
            }
        }
        if ((offlinedata) && (estimpredstim) && (zerocr)) // 3nd method - end
        {            
            diffZC=zcdifference(pos+(hnt->imlength+hnt->stlength)*i);
            plvzc+=plv;
            if ((!phasepr) && (!hronar))
            {
                offlinevaluation(pos,i,3);
                hnt->averagedegdiffZC+=hnt->phsdiffdeg;
            } else
            {
                hnt->phsdiffdeg = (360 / (hnt->srfr/hnt->osfr))*abs(hnt->optphase-hnt->optphasezc);
                if (hnt->phsdiffdeg>180)
                    hnt->phsdiffdeg=abs(360-hnt->phsdiffdeg);
                hnt->averagedegdiffZC+=hnt->phsdiffdeg;
                arrphasediffzc[i]=hnt->phsdiffdeg;
            }
            if (hnt->inphase)
                proczc=hnt->optphasediff/diffZC;
            else
                proczc=diffZC/hnt->optphasediff;
           //  proczc = (double)(180-hnt->phsdiffdeg)/180;
        }
        if (fftpredict) // 4th method - prediction by phase of dominant frequency from FFT
        {
           fftphasepred(pos+(hnt->imlength+hnt->stlength)*i,hnt->imlength*hnt->extinterval);
           hnt->optphasefft=hnt->phase;
           int tp = hnt->srfr/hnt->osfr;
           if (hnt->optphasefft > tp)
               hnt->optphasefft = hnt->optphasefft % tp;
           if (hnt->optphasefft < 0)
               hnt->optphasefft = tp - abs(hnt->optphasefft);
           if ((offlinedata) && (estimpredstim))
           {
               arrsyncfft[i]=gethindex();
               plvfft+=plv;
               diffFFT=hnt->phsdif;
               /* if ((!phasepr) && (!hronar) && (!zerocr))
               {
                   offlinevaluation(pos,i,4);
                   hnt->averagedegdiffFFT+=hnt->phsdiffdeg;
               } else */
               {
                   hnt->phsdiffdeg = (360 / (hnt->srfr/hnt->osfr))*abs(hnt->optphase-hnt->optphasefft);
                   if (hnt->phsdiffdeg>180)
                       hnt->phsdiffdeg=abs(360-hnt->phsdiffdeg);
                   hnt->averagedegdiffFFT+=hnt->phsdiffdeg;
                   arrphasedifffft[i]=hnt->phsdiffdeg;
               }
               if (hnt->inphase)
                   procfft=hnt->optphasediff/diffFFT;
               else
                   procfft=diffFFT/hnt->optphasediff;
           }
        }
        if (randnoisestim)
        {
            hnt->phase = 0;
            hnt->optphasehr = hnt->phase;
            fillrandstim(pos+(hnt->imlength+hnt->stlength)*i,hnt->stlength);
        }
        if (randphasestim)
        {
            hnt->phase = qrand()%((int)(hnt->srfr/hnt->osfr));
            hnt->optphasehr = hnt->phase;
            fillstim(pos+(hnt->imlength+hnt->stlength)*i,hnt->stlength);
        }
   /**/ currtime = ssdt + timer.nsecsElapsed() / 1000000.0;
        averopttime+=currtime;
        QString rel;
        if (offlinedata)
        {
            if (relations[i]==1)
                rel="  In   |";
            else
                rel="Anti |";
        } else
        {
            if (relations[totalstims+stims]==1)
                rel="  In   |";
            else
                rel="Anti |";
        }
        if (ui->radioButton_3->isChecked())
            rel="Noise |";
        if (ui->radioButton_4->isChecked())
            rel="Rand |";
        if (offlinedata)
            mw->printdata("TR_"+QString::number(i+1)+" "+rel+" Optimization time (msec): " + QString::number(currtime,4,2)+"     Opt Phase shift: "+QString::number(hnt->optphase)+"     Pred Phase shifts: "+QString::number(hnt->optphasehr)+" "+QString::number(hnt->optphasear)+" "+QString::number(hnt->optphasezc)+" "+QString::number(hnt->optphasefft));
        else
            mw->printdata("TR_"+QString::number(stims+1)+" "+rel+" Optimization time (msec): " + QString::number(currtime,4,2)+"     Opt Phase shift: "+QString::number(hnt->optphase)+"     Pred Phase shifts: "+QString::number(hnt->optphasehr)+" "+QString::number(hnt->optphasear)+" "+QString::number(hnt->optphasezc)+" "+QString::number(hnt->optphasefft));
        if (!offlinedata)
        {
            // transdelay ~ 37 points, should be estimated from differebce between data here and raw data
            receivedelay=counter-hnt->imlength + transdelay;
            readypredict=true;
          //  qDebug()<<counter-hnt->imlength;
            writedaq(pos);
        }
        if ((offlinedata) && (estimateclear))
            clearoptphase[i]=hnt->optphase;
        hnt->totalprocregr+=procregr;
        hnt->totalprocpred+=procpred;
        hnt->totalproczc+=proczc;
        hnt->totalprocfft+=procfft;
        hnt->totaldiffHR+=diffHR;
        hnt->totaldiffZC+=diffZC;
        hnt->totaldiffHRonAR+=diffHRonAR;
        hnt->totaldiffFFT+=diffFFT;
        hnt->totalsynchr+=arrsynchr[i];
        hnt->totalsyncar+=arrsyncar[i];
        hnt->totalsynczc+=arrsynczc[i];
        hnt->totalsyncfft+=arrsyncfft[i];
        if (estimation)
            clearstim(pos+(hnt->imlength+hnt->stlength)*i,hnt->stlength,1);
        if ((offlinedata) && (optstim))
            filloptstim(pos+(hnt->imlength+hnt->stlength)*i,hnt->stlength,hnt->optphase, 6);
    }
    // end of stim sequence


    if ((offlinedata) && (estimateclear))
    {
        estimateclear=false;
        estimatenoisy=true;
    }
    if (filtereddata)
          ui->widget->graph(1)->setData(arrc.xc, arrc.amp2);
    if (offlinedata)
    {
    hnt->totalprocpred/=hnt->numst;   
    hnt->totalprocregr/=hnt->numst;
    hnt->totalproczc/=hnt->numst;
    hnt->totalprocfft/=hnt->numst;
    plvzc/=hnt->numst; plvhr/=hnt->numst; plvHRonAR/=hnt->numst; plvfft/=hnt->numst;
    hnt->totaldiffHR/=hnt->stlength;
    hnt->totaldiffZC/=hnt->stlength;
    hnt->totaldiffHRonAR/=hnt->stlength;
    hnt->totaldiffFFT/=hnt->stlength;
    hnt->averagedegdiffHR/=hnt->numst;
    hnt->averagedegdiffAR/=hnt->numst;
    hnt->averagedegdiffZC/=hnt->numst;
    hnt->averagedegdiffFFT/=hnt->numst;
    hnt->totalsynchr/=hnt->numst;
    hnt->totalsyncar/=hnt->numst;
    hnt->totalsynczc/=hnt->numst;
    hnt->totalsyncfft/=hnt->numst;
    for (int i=0; i<nms; i++)
    {        
        arrphasediffhr[i]=abs(arrphasediffhr[i]-hnt->averagedegdiffHR);
        arrphasediffar[i]=abs(arrphasediffar[i]-hnt->averagedegdiffAR);
        arrphasediffzc[i]=abs(arrphasediffzc[i]-hnt->averagedegdiffZC);
        arrphasedifffft[i]=abs(arrphasedifffft[i]-hnt->averagedegdiffFFT);
        hnt->degdevHR+=arrphasediffhr[i];
        hnt->degdevAR+=arrphasediffar[i];
        hnt->degdevZC+=arrphasediffzc[i];
        hnt->degdevFFT+=arrphasedifffft[i];
    }

 /*   for (int i=0; i<nms; i++)
        hnt->fillphasebins(qDegreesToRadians((double)arrphasediffhr[i]));
    cout<<endl;
    for (int i = hnt->binnums/2; i < hnt->binnums; i++)
        cout<<hnt->phasebins[i]<<" ";
    cout<<endl;*/

    hnt->degdevHR/=(nms); hnt->degdevAR/=(nms); hnt->degdevZC/=(nms); hnt->degdevFFT/=(nms);
    }
    //hnt->degdevHR=sqrt(hnt->degdevHR);
    //hnt->degdevAR=sqrt(hnt->degdevAR);
    //hnt->degdevZC=sqrt(hnt->degdevZC);
   // qDebug()<<endl<<hnt->degdevHR<<" "<<hnt->degdevAR<<" "<<hnt->degdevZC;
    if (offlinedata)
    {        
        mw->printdata("Average opt. time (msec): "+ QString::number(averopttime/hnt->numst));
        mw->printdata("==== End of " + QString::number(hnt->numst) + " stims sequence ====");
        mw->printdata("");
        averopttime=0;
    }
    else if (stims==hnt->numst-1)
    {
        mw->printdata("");
        mw->printdata("Average opt. time (msec): "+ QString::number(averopttime/hnt->numst));        
        averopttime=0;
    }    
    if (offlinedata)
    {
        ui->widget->graph(0)->setData(arrc.xc, arrc.amp1);
        ui->widget->graph(3)->setData(arrc.xc, arrc.of1); // HT
        if (!hideardata)
            ui->widget->graph(2)->setData(arrc.xc, arrc.of2);
        ui->widget->graph(4)->setData(arrc.xc, arrc.of3); // AR
        ui->widget->graph(5)->setData(arrc.xc, arrc.of4); // ZC
        ui->widget->graph(6)->setData(arrc.xc, arrc.of6); // FFT
        ui->widget->graph(7)->setData(arrc.xc, arrc.of5); // opt
        printtoresultsync(QString::number(hnt->totalsynchr,'f',2)+"  "+QString::number(hnt->totalsyncar,'f',2)+"  "+QString::number(hnt->totalsynczc,'f',2)+"  "+QString::number(hnt->totalsyncfft,'f',2));
        printtoresultbox(QString::number(100*hnt->totalprocpred,'f',1)+"  "+QString::number(100*hnt->totalprocregr,'f',1)+"  "+QString::number(100*hnt->totalproczc,'f',1)+"  "+QString::number(100*hnt->totalprocfft,'f',1));
        printtoresultstring("ADD   HT: "+QString::number(hnt->averagedegdiffHR,'f',0)+"+/-"+QString::number((int)(hnt->degdevHR))+"    AR: "+QString::number(hnt->averagedegdiffAR,'f',0)+"+/-"+QString::number((int)(hnt->degdevAR))+"    ZC: "+QString::number(hnt->averagedegdiffZC,'f',0)+"+/-"+QString::number((int)(hnt->degdevZC))+"    FFT: "+QString::number(hnt->averagedegdiffFFT,'f',0)+"+/-"+QString::number((int)(hnt->degdevFFT))+"          PLV   HT: "+QString::number(plvhr,'f',3)+"    AR: "+QString::number(plvHRonAR,'f',3)+"    ZC: "+QString::number(plvzc,'f',3)+"    FFT: "+QString::number(plvfft,'f',3));
        saverestofile();
    }
}

void plotwindow::parseoptparams()
{
    QFile inputFile("C:/optparams.txt");
    inputFile.open(QIODevice::ReadOnly);
    QTextStream fin(&inputFile);
    QStringList sl;

    int subj=0;
    while (subj<20)
    {
        QString l1 = fin.readLine();
        QString l2 = fin.readLine();
        QString l3 = fin.readLine();
        QString l4 = fin.readLine();

        sl = l1.split(QRegExp("\\s+"), QString::SkipEmptyParts);
        for (int k=0; k < 4; k++)
            arrres[subj][k]=sl[k].toInt();
        sl = l2.split(QRegExp("\\s+"), QString::SkipEmptyParts);
        for (int k=0; k < 4; k++)
            arrres[subj][k+4]=sl[k].toInt();
        sl = l3.split(QRegExp("\\s+"), QString::SkipEmptyParts);
        for (int k=0; k < 4; k++)
            arrres[subj][k+8]=sl[k].toInt();
        sl = l4.split(QRegExp("\\s+"), QString::SkipEmptyParts);
        for (int k=0; k < 4; k++)
            arrres[subj][k+12]=sl[k].toInt();

        subj++;
    }

    inputFile.close();

    QFile outputFile("C:/optparamslines.csv");
    outputFile.open(QIODevice::Append);
    QTextStream fout(&outputFile);
    for (int i=0; i<subj; i++)
    {
        for (int j=0; j<16; j++)
        {
            fout << arrres[i][j];
            if (j<15) fout << ',';
        }
        fout << endl;
    }

    outputFile.close();
}

void plotwindow::parseresfromfile()
{
    QFile inputFile("C:/closedloopres.dat");
    inputFile.open(QIODevice::ReadOnly);
    QTextStream fin(&inputFile);
    QStringList sl;

    int subj=0;
    while (subj<17)
    {
        {
          //  for (int i=0; i<44; i++)
          //      fin.readLine();

            QString l1 = fin.readLine();
            QString l2 = fin.readLine();
            QString l3 = fin.readLine();
            QString l4 = fin.readLine();

            sl = l1.split(QRegExp("\\s+"), QString::SkipEmptyParts);
            for (int k=0; k < 4; k++)
                arrres[subj][k]=sl[k].toDouble();
            sl = l2.split(QRegExp("\\s+"), QString::SkipEmptyParts);
            for (int k=0; k < 4; k++)
                arrres[subj][k+4]=sl[k].toDouble();
            sl = l3.split(QRegExp("\\s+"), QString::SkipEmptyParts);
            for (int k=0; k < 4; k++)
                arrres[subj][k+8]=sl[k].toDouble();
            sl = l4.split(QRegExp("\\s+"), QString::SkipEmptyParts);
            for (int k=0; k < 4; k++)
                arrres[subj][k+12]=sl[k].toDouble();

         //   for (int i=0; i<12; i++)
         //       fin.readLine();

        }
        subj++;
    }

    inputFile.close();

    QFile outputFile("C:/arteegres.txt");
    outputFile.open(QIODevice::Append);
    QTextStream fout(&outputFile);
    for (int i=0; i<subj; i++)
    {
        for (int j=0; j<16; j++)
        {
            fout << arrres[i][j];
            if (j<15) fout << ',';
        }
        fout << endl;
    }

    outputFile.close();

}


void plotwindow::saverestofile()
{
    QFile outputFile("C:/closedloopres.dat");
    outputFile.open(QIODevice::Append);
    QTextStream fout(&outputFile);
  //  fout << endl;       
    fout << QString::number(hnt->averagedegdiffHR,'f',1) << " " << QString::number(hnt->averagedegdiffAR,'f',1) << " " << QString::number(hnt->averagedegdiffZC,'f',1) << " " << QString::number(hnt->averagedegdiffFFT,'f',1) << endl;
    fout << QString::number(100*hnt->totalprocpred,'f',2) << " " << QString::number(100*hnt->totalprocregr,'f',2) << " " << QString::number(100*hnt->totalproczc,'f',2) << " " << QString::number(100*hnt->totalprocfft,'f',2) << endl;
    fout << QString::number(plvhr,'f',3) << " " << QString::number(plvHRonAR,'f',3) <<" "<< QString::number(plvzc,'f',3) <<" "<< QString::number(plvfft,'f',3) << endl;
    fout << QString::number(hnt->totalsynchr,'f',3) << " " << QString::number(hnt->totalsyncar,'f',3) << " " << QString::number(hnt->totalsynczc,'f',3) << " " << QString::number(hnt->totalsyncfft,'f',3) << endl;
    outputFile.close();
}

void plotwindow::sequencestim()
{    
    hnt->init();
    if ((autoregdegree>=hnt->imlength) && (hronar))
    {
        QMessageBox msgBox;
        msgBox.setText("Autoregression order must be < imaging interval!");
        msgBox.exec();
    } else
    if (hnt->imlength*hnt->extinterval<=hnt->lfilt)
    {
        QMessageBox msgBox;
        msgBox.setText("Filter length must be < extraction interval!");
        msgBox.exec();
    }
    else
    {
        corrprednumb=0;
        hnt->posstim=ui->spinBox->value();
        hnt->lfilt=ui->spinBox_2->value();
        hnt->imlength=ui->spinBox_5->value();
        hnt->imstprop=ui->doubleSpinBox_3->value();
        hnt->stlength=hnt->imlength*hnt->imstprop;
        hnt->osfr=ui->doubleSpinBox->value();
        hnt->phase=ui->spinBox_3->value();
       // hnt->initdelay=ui->spinBox_3->value();
        hnt->noise=ui->spinBox_4->value();
        hnt->stampl=ui->spinBox_6->value();
        hnt->thfr=ui->doubleSpinBox_2->value();        
        hnt->srfr=ui->spinBox_8->value();
        hnt->mind=0;//(int)(-hnt->srfr/hnt->osfr/2);
        hnt->maxd=(int)(hnt->srfr/hnt->osfr);

        if (offlinedata)
        {
            totalstims=0;
            if (randomstims)
                setrandomrelations();
            else
            if (hnt->inphase)
                for (int i=0; i<hnt->numst; i++)
                    relations[i]=1;
            else
                for (int i=0; i<hnt->numst; i++)
                    relations[i]=-1;
        }        

        strategy(hnt->posstim);

        ui->doubleSpinBox_2->setValue(hnt->thfr);
        ui->lcdNumber_4->display(hnt->srfr/hnt->osfr);
        if (offlinedata)                 
            ui->widget->replot();        
        if (offlinedata)
            offlinedata=false;
    }
}
void plotwindow::clearstim()
{
    for (int i=0; i<NMAX; i++) {
        arrc.of1[i]=0;
        arrc.of2[i]=0;
        arrc.of3[i]=0;
        arrc.of4[i]=0;
        arrc.of5[i]=0;
        arrc.of6[i]=0;
    }
    ui->widget->graph(3)->setData(arrc.xc, arrc.of1);
    ui->widget->graph(2)->setData(arrc.xc, arrc.of2);
    ui->widget->graph(4)->setData(arrc.xc, arrc.of3);
    ui->widget->graph(5)->setData(arrc.xc, arrc.of4);
    ui->widget->graph(6)->setData(arrc.xc, arrc.of6);
    ui->widget->graph(7)->setData(arrc.xc, arrc.of5);
    ui->widget->replot();
}

void plotwindow::cleareeg(int start, int end)
{
    for (int i=start; i<end; i++)
        arrc.amp1[i]=0;
    ui->widget->graph(0)->setData(arrc.xc, arrc.amp1);
    ui->widget->replot();
}

void plotwindow::clearfiltered()
{
    for (int i=0; i<NMAX; i++)
        arrc.amp2[i]=0;
    ui->widget->graph(1)->setData(arrc.xc, arrc.amp2);
    ui->widget->replot();
}

void plotwindow::draw(int start)
{
    int k=hnt->imlength;  
    for (int i=0; i<k; i++)
        arrc.amp1[start+i]=arrc.amp0[i];
    if (usefiltering) {
        if (filtfir)
            filterdata(start,hnt->imlength);
        else if (recurbutter)
            recurbutterfilter(start,hnt->imlength);
        else if (zerobutter)
            zerophasefilt(start,hnt->imlength);
        }
}

void plotwindow::delayt(int temp)
{
    QTime dieTime = QTime::currentTime().addMSecs(temp);
    while (QTime::currentTime() < dieTime)
        QCoreApplication::processEvents(QEventLoop::AllEvents, 100);
}

void plotwindow::bcidata(double d1)
{
   // if ((artifactremove) && (abs(d1)>1000))
   //     arrc.amp0[counter]=0;
   // else
    arrc.amp0[counter]=d1;
    counter++;
    // restore data before exact start of stimulation on electrodes, add and shift pre-stim to left
    if ((!(imagingafterstim)) && (counter==hnt->imlength+receivedelay) && (readypredict)) // !! try transdelay instead of receive here !!
    {
        for (int i=0; i<hnt->imlength; i++)
            arrc.amp1[startpos+stims*(2*hnt->imlength+hnt->stlength)+i]=arrc.amp0[receivedelay+i];
        if (usefiltering)
            zerophasefilt(startpos+stims*(2*hnt->imlength+hnt->stlength),hnt->imlength);
        readypredict=false;
    }
    // skip delayed data after stim is finished
    if ((imagingafterstim) && (counter==receivedelay) && (skipdelayedata))
    {
        counter=0;
        skipdelayedata=false;
    }
}

void plotwindow::makestimoneeg(int start)
{          
    if ((stims<hnt->numst) && (!addmode))
    {
       /* if (applyssd)
        {
            timer.restart();
            getssd(ssdcomp,hnt->imlength,"",-1);
            ssdt = timer.nsecsElapsed() / 1000000.0;
            //mw->printdata("SSD projection time (msec): " + QString::number(ssdt));
        } */
        if (imagingafterstim)
        {
            draw(start+stims*(2*hnt->imlength+hnt->stlength)-hnt->imlength);
            if (stims<hnt->numst)
            {
                int t = qrand() % iti1+iti2;
                mw->printdata("Inter Trial Interval (msec): "+QString::number(t));
                delayt(t);
            }
        }
        if (!imagingafterstim)
        {
            draw(start+stims*(2*hnt->imlength+hnt->stlength));
            hnt->posstim=start+(stims+1)*(2*hnt->imlength+hnt->stlength)-hnt->stlength-hnt->imlength;
            ui->spinBox->setValue(hnt->posstim);
           // if (relations[totalstims+stims]==1)
           // {
           //     ui->radioButton->setChecked(true);
              //  mw->printdata("--------- In Phase ---------");
          //  }
          //  else if (relations[totalstims+stims]==-1)
          //  {
          //      ui->radioButton_2->setChecked(true);
              //  mw->printdata("-------- In Opposite --------");
          //  }
            sequencestim();
            stims++;
        }
    } else
    if (stims==hnt->numst)
    {
        draw(start+stims*(2*hnt->imlength+hnt->stlength)-hnt->imlength);

        tim->stop();

        ui->widget->graph(0)->setData(arrc.xc, arrc.amp1);
        ui->widget->graph(3)->setData(arrc.xc, arrc.of1);
        if (!hideardata)
            ui->widget->graph(2)->setData(arrc.xc, arrc.of2);
        ui->widget->graph(4)->setData(arrc.xc, arrc.of3);
        ui->widget->graph(5)->setData(arrc.xc, arrc.of4);
        ui->widget->replot();
        mw->printdata("==== End of " + QString::number(stims) + " stims sequence ====");
        mw->printdata("");
        totalstims+=hnt->numst;
    };
    if ((addmodeon) && (addmoderec))
    {        
        draw(start+recparts*hnt->imlength);
        recparts++;
        ui->widget->graph(0)->setData(arrc.xc, arrc.amp1);
        ui->widget->replot();
    };
}

void plotwindow::timerUpdate()
{
     // wait till new fragment
     if (counter>=hnt->imlength)
     {        
         makestimoneeg(startpos);
         imagingafterstim=!imagingafterstim;
         counter=0; skipdelayedata=true;
     }
     if ((addmodeon) && (!addmoderec))
     {
         stimulationonly(startpos+stims*hnt->stlength);
         stims++;
         if (stims==hnt->numst)
         {
             tim->stop();
             addmodeon=false;
         }
     }
}

void plotwindow::savesingledatatofile(QString fname, double iaf)
{
    QFile outputFile(fname);
    outputFile.open(QIODevice::WriteOnly);
    QTextStream fout(&outputFile);
    fout << iaf << " " << trialnums << endl;

    for (int i=0; i<310000; i++)
        fout << arrc.amp1[i] << " ";
    fout << endl;
}

void plotwindow::savecspdatatofile(QString fname, double iaf)
{
    QFile outputFile(fname);
    outputFile.open(QIODevice::WriteOnly);
    QTextStream fout(&outputFile);
    fout << iaf << " " << trialnums << endl;

    for (int i=0; i<310000; i++)
        fout << arrc.amp2[i] << " ";
    fout << endl;
}

void plotwindow::savedatatofile(QString fname)
{
    QFile outputFile(fname);
    outputFile.open(QIODevice::WriteOnly);
    QTextStream fout(&outputFile);  
    fout << endl;

    fout << QString::fromStdString(hnt->channel)<<" ";
    for (int i=0; i<NMAX; i++)
        fout << arrc.amp1[i] << " ";
    fout << endl;


    fout << "stim1: ";
    for (int i=0; i<NMAX; i++)
        fout << arrc.of1[i] << " ";
    fout << endl;

    fout << "stim2: ";
    for (int i=0; i<NMAX; i++)
        fout << arrc.of3[i] << " ";
    fout << endl;

    fout << "stim3: ";
    for (int i=0; i<NMAX; i++)
        fout << arrc.of4[i] << " ";
    fout << endl;

    fout << "relations: ";
    for (int i=0; i<totalstims; i++)
    {
        fout << relations[i] << " ";
      //  cout << relations[i] << " ";
    }
   // cout << endl;
    fout << endl;

    if (estimatenoisy)
    {
        fout << endl;
        for (int i = 0; i < strLst1.size(); ++i)
            fout << strLst1.at(i) << '\n';
        fout << '\n';
        for (int i = 0; i < strLst2.size(); ++i)
            fout << strLst2.at(i) << '\n';
        fout << '\n';
        for (int i = 0; i < strLst3.size(); ++i)
            fout << strLst3.at(i) << '\n';
        fout << '\n';
    }

    outputFile.close();
}

void plotwindow::writerelationstofile(QString fname)
{
    QFile outputFile(fname);
    outputFile.open(QIODevice::WriteOnly);
    QTextStream fout(&outputFile);
    for (int i=0; i<500; i++)
    {
        fout << relations[i] << " ";
        if ((i+1) % 50 == 0)
            fout << endl;
    }
    outputFile.close();
}

void plotwindow::loadstimfromfile(QString fname)
{
    QFile inputFile(fname);
    inputFile.open(QIODevice::ReadOnly);
    QTextStream fin(&inputFile);        
    QStringList sl;

    QString line = fin.readLine();
    line = fin.readLine();    
    int leng=0;

    if (line.length()>0)
    {
        sl = line.split(QRegExp("\\s+"), QString::SkipEmptyParts);
        leng=sl.size();
        if (leng>NMAX)
            leng=NMAX;
        for (int k=0; k < leng; k++)
            arrc.amp1[k]=sl[k].toDouble();
      /*  QMessageBox msgBox;
        msgBox.setText("1");
        msgBox.exec(); */
    }        

    line = fin.readLine();    
    if (line.length()>0)
    {
        sl = line.split(QRegExp("\\s+"), QString::SkipEmptyParts);
        leng=sl.size();
        if (leng>NMAX)
            leng=NMAX;
        for (int k=1; k < leng; k++)
            arrc.of1[k-1]=sl[k].toDouble();
    }

    line = fin.readLine();
    if (line.length()>0)
    {
        sl = line.split(QRegExp("\\s+"), QString::SkipEmptyParts);
        leng=sl.size();
        if (leng>NMAX)
            leng=NMAX;
        for (int k=1; k < leng; k++)
            arrc.of3[k-1]=sl[k].toDouble();
    }

    line = fin.readLine();
    if (line.length()>0)
    {        
        sl = line.split(QRegExp("\\s+"), QString::SkipEmptyParts);
        leng=sl.size();
        if (leng>NMAX)
            leng=NMAX;
        for (int k=1; k < leng; k++)
            arrc.of4[k-1]=sl[k].toDouble();
    }   

    for (int k=1; k < NMAX; k++)
        arrc.of5[k-1]=0;

    line = fin.readLine();
    if (line.length()>0)
    {
        sl = line.split(QRegExp("\\s+"), QString::SkipEmptyParts);
        leng=sl.size();
        totalstims = leng-1;
        for (int k=1; k < leng; k++)
            relations[k-1]=sl[k].toInt();
    }

//    for (int i=0; i<leng-1; i++)
//        cout<<relations[i]<<" ";
    inputFile.close();
}

void plotwindow::savescreen()
{
    QString fname=QFileDialog::getSaveFileName(this,tr("Save File as"),"C://","Data file (*.dat);;All files (*.*)");
    if (fname!="")
    {
       ui->widget->savePng(fname);
    }
}

void plotwindow::refresh()
{
    ui->doubleSpinBox->setValue(hnt->osfr);
    ui->doubleSpinBox_3->setValue(hnt->imstprop);
    ui->spinBox_5->setValue(hnt->imlength);
    ui->spinBox_7->setValue(hnt->numst);
}

void plotwindow::slSettings()
{          
   sw->pwd=this;
   sw->setFixedSize(666,694);
   sw->init();
   sw->show();
}

int plotwindow::estimateoptlength(int n, int l1, int l2, int pos)
{
    estimation=true;
    int k = hnt->numst;
    int p = hnt->posstim;
    int iml = hnt->imlength;
    hnt->numst=n;
    hnt->posstim=pos;
    double err;
    if (hnt->inphase)
        err=1000000;
    else
        err=0;
    int optl=0;
    for (int i=l1; i<=l2; i+=5)
    {        
        hnt->imlength=i;
        hnt->stlength=hnt->imlength*hnt->imstprop;
        offlinedata=true;
        strategy(hnt->posstim);
        if ((hnt->inphase) && (hnt->totalprocpred<err))
        {
            err=hnt->totalprocpred;
            optl=i;
        }  else
        if (!(hnt->inphase) && (hnt->totalprocpred>err))
        {
            err=hnt->totalprocpred;
            optl=i;
        }
    }
    hnt->numst=k;
    hnt->posstim=p;
    hnt->imlength=iml;
    estimation=false;
    return optl;
}

double plotwindow::estimateoptprop(int n, double p1, double p2, int pos)
{
    estimation=true;
    int k = hnt->numst;
    int p = hnt->posstim;
    double prp=hnt->imstprop;
    hnt->numst=n;
    hnt->posstim=pos;
    double err;
    if (hnt->inphase)
        err=1000000;
    else
        err=0;
    double optprop=0;
    double pt=p1;
    while (pt<=p2)
    {
        hnt->imstprop=pt;
        hnt->stlength=hnt->imlength*hnt->imstprop;
        offlinedata=true;
        strategy(hnt->posstim);
        if ((hnt->inphase) && (hnt->totalprocpred<err))
        {
            err=hnt->totalprocpred;
            optprop=pt;
        } else
        if (!(hnt->inphase) && (hnt->totalprocpred>err))
        {
            err=hnt->totalprocpred;
            optprop=pt;
        }
        pt+=0.05;
    }
    hnt->numst=k;
    hnt->posstim=p;
    hnt->imstprop=prp;
    estimation=false;
    return optprop;
}

void plotwindow::setaddmode(bool f)
{
    addmoderec=f;
    if (f)
        ui->checkBox_4->setText("Record-only mode");
    else
        ui->checkBox_4->setText("Stimulation-only mode");
}

int plotwindow::maxabsstim(int pos)
{
    int maxel=0;
    for (int j=pos; j<pos+hnt->stlength; j++)
        if (abs(arrc.of1[j])>maxel)
            maxel=abs(arrc.of1[j])+1;
    return maxel;
}

void plotwindow::stimulationonly(int pos)
{
    hnt->stlength=ui->spinBox_5->value()*ui->doubleSpinBox_3->value();
    fillstim(pos,hnt->stlength);
    writedaq(pos);
    ui->widget->graph(3)->setData(arrc.xc, arrc.of1);
    ui->widget->replot();
}

void plotwindow::delay(int temp)
{
    QTime dieTime = QTime::currentTime().addMSecs(temp);
    while (QTime::currentTime() < dieTime)
        QCoreApplication::processEvents(QEventLoop::AllEvents, 1);
}

void plotwindow::writedaq(int pos)
{
    timer.restart();
    int length = hnt->stlength;
    if (!offlinedata)
    {
        int32* k;

    //    DAQmxStopTask(taskHandle);
        float64* arr = new float64[length];
        for (int i=pos; i<pos+length-1; i++)
            arr[i-pos]=offsetcomp+arrc.of1[i]/daqscalevalue;
        arr[length-1]=0;
      //  DAQmxWriteAnalogF64(taskHandle,length,1,0,DAQmx_Val_GroupByChannel,arr,k,NULL);
        int t = timer.nsecsElapsed() / 1000000.0;
        Sleep(3000-t);
      //  qDebug()<<t;
    }
}

void plotwindow::maxentropyf(double *input, int length, int degree, double **ar, double *per, double *pef, double *h, double *g)
{
    int j,n,nn,jj;
    double sn,sd;
    double t1,t2;

    for (j=1; j<=length; j++) {
        pef[j] = 0;
        per[j] = 0;
    }

    for (nn=2; nn<=degree+1; nn++) {
        n  = nn - 2;
        sn = 0.0;
        sd = 0.0;
        jj = length - n - 1;
        for (j=1; j<=jj; j++) {
            t1 = input[j+n] + pef[j];
            t2 = input[j-1] + per[j];
            sn -= 2.0 * t1 * t2;
            sd += (t1 * t1) + (t2 * t2);
        }
        g[nn] = sn / sd;
        t1 = g[nn];
        if (n != 0) {
            for (j=2; j<nn; j++)
                h[j] = g[j] + (t1 * g[n - j + 3]);
            for (j=2; j<nn; j++)
                g[j] = h[j];
            jj--;
        }
        for (j=1; j<=jj; j++) {
            per[j] += (t1 * pef[j]) + (t1 * input[j+nn-2]);
            pef[j]  = pef[j+1] + (t1 * per[j+1]) + (t1 * input[j]);
        }

        for (j=2; j<=nn; j++)
            ar[nn-1][j-1] = g[j];
    }
}

void plotwindow::autoregress(double *input, int length, int degree, double *coeffs)
{
    double mean;
    int i, t;
    double *w=NULL;      // Input series - mean
    double *h=NULL;
    double *g=NULL;
    double *per=NULL;
    double *pef=NULL;
    double **ar=NULL;

    w = (double *)malloc(length*sizeof(double));
    h = (double *)malloc((degree+1)*sizeof(double));
    g = (double *)malloc((degree+2)*sizeof(double));
    per = (double *)malloc((length+1)*sizeof(double));
    pef = (double *)malloc((length+1)*sizeof(double));
    ar = (double **)malloc((degree+1)*sizeof(double*));

    for (i=0; i<degree+1; i++)
        ar[i] = (double *)malloc((degree+1)*sizeof(double));
    mean = 0.0;
    for (t=0; t<length; t++)
       mean += input[t];
    mean /= (double)length;
    for (t=0; t<length; t++)
       w[t] = input[t] - mean;

    maxentropyf(input,length,degree,ar,per,pef,h,g);
    for (i=1; i<=degree; i++)
        coeffs[i-1] = -ar[degree][i];
}

void plotwindow::runautoreg(int pos)
{
    maxentropy=true; leastsqr=false;
   // qDebug()<<autoregdegree;
  //  autoregdegree=hnt->imlength/4;
    double * inp = new double [hnt->imlength];
    for (int i=0; i<hnt->imlength; i++)    
         inp[i]=arrc.amp1[pos-hnt->imlength+i]; 
    autoregress(inp,hnt->imlength,autoregdegree,autoregcoeffs);
 //   for (int i=0; i<autoregdegree; i++)
 //       qDebug()<<autoregcoeffs[i];
    double nextpoint = 0;
    ar_error = 0;
    int np = hnt->stlength*hnt->arinterval;
    for (int i=0; i<autoregdegree; i++)
        arrc.of2[pos-i]=arrc.amp1[pos-i];
    for (int i=pos; i<pos+np; i++)
    {
        for (int j=0; j<autoregdegree; j++)
            nextpoint+=autoregcoeffs[j]*arrc.of2[i-j];
      //  qDebug()<<nextpoint;
        arrc.of2[i+1]=nextpoint;
        ar_error += pow(arrc.of2[i+1] - arrc.amp1[i+1],2);
        nextpoint=0;
    }
    ar_error = pow(ar_error,0.5)/np;
  //  qDebug()<<ar_error;

    if ((zerobutter) && (usefiltering))
    for (int i=pos; i<pos+hnt->stlength*hnt->arinterval+1; i++)
        arrc.of2[i]/=10000;
    for (int i=0; i<hnt->imlength; i++)
        arrc.of2[pos-i]=0;
  //  ui->widget->graph(2)->setData(arrc.xc, arrc.of2);

}

void plotwindow::hilbonar(int pos)
{
    hronarrun=true;
    if (hnt->inphase)
        minphasediff(pos,hnt->stlength*hnt->arinterval);
    else
        maxphasediff(pos,hnt->stlength*hnt->arinterval);
    clearstim(pos,hnt->stlength,0);
    hronarrun=false;
    fillstimforregr(pos,hnt->stlength);    
}

void plotwindow::ssdtimUpdate()
{
    if ((extractingssd) && (indexes[chnums-1]==hnt->imlength))
    {
        timer.restart();
        ssd("",0);
        mw->printdata("SSD extraction time (msec): " + QString::number(timer.nsecsElapsed() / 1000000.0));
        mw->printdata("SSD extracted.");
        hnt->imlength=sw->siml;
        extractingssd=false;
        sw->ssdready();
    }
}

void plotwindow::getrawdata(int chn, double val)
{    
    if ((chn==chnums-1) && (indexes[chn]==hnt->imlength))
    {                                       // current index of data in channel "chn"
      /*  cout<<endl;                       // if last channel obtained "imlength" number of data - resetting indexes
        for (int j=0; j<chnums; j++)
        {
            for (int k=0; k<10; k++)
                cout<<rawdata[j][k]<<" ";
            cout<<endl;
        }
        cout<<endl;
        QApplication::quit(); */      
        for (int i=0; i<chnums; i++)
            indexes[i]=0;
    }
    rawdata[chn][indexes[chn]]=val;
    indexes[chn]++;
}

void plotwindow::filterdata(int posstim, int length)
{
    Filter *myfilt;
    myfilt = new Filter(LPF, flength, fsampr, flowpass);
    double temp[flength/2];
    for (int i = 0; i < flength/2; i++)
        temp[i]=arrc.amp1[posstim+i-flength/2];
    for (int i = 0; i < length; i++)
        arrc.amp1[posstim+i-flength/2] = myfilt->do_sample(arrc.amp1[posstim+i]);
    for (int i = 0; i < flength/2; i++)
        arrc.amp1[posstim+i-flength/2]=temp[i];
}

void plotwindow::filterar(int posstim, int length)
{
    Filter *myfilt;
    myfilt = new Filter(LPF, flength, fsampr, flowpass);
    double temp[flength/2];
    for (int i = 0; i < flength/2; i++)
        temp[i]=arrc.of2[posstim+i-flength/2];
    for (int i = 0; i < length; i++)
        arrc.of2[posstim+i-flength/2] = myfilt->do_sample(arrc.of2[posstim+i]);
    for (int i = 0; i < flength/2; i++)
        arrc.of2[posstim+i-flength/2]=temp[i];   
}

void plotwindow::zerocrossingstim(int posstim, int length)
{
    int prevx0=0;
    int lastx0=0;
    int phase=0; int t=0;
    int t_amp=3;

    for (int i=posstim-length; i<posstim; i++)    
    if ((arrc.amp2[i]>0) && (arrc.amp2[i+1]<0))
    {
        prevx0=lastx0;
        lastx0=i;
        t=(hnt->srfr/hnt->osfr)/2;
    } else
    if ((arrc.amp2[i]<0) && (arrc.amp2[i+1]>0))
    {
        prevx0=lastx0;
        lastx0=i;
        t=0;
    }

  //  qDebug()<<prevx0<<lastx0;

    if (prevx0 > 0)
    {
        bool fl = true;
        for (int i=lastx0; i<posstim; i++)
            if (abs(arrc.amp2[i])>t_amp)
                fl = false;

        if (fl)
        {
            lastx0 = prevx0;
            t = (t>0) ? 0 : (hnt->srfr/hnt->osfr)/2;
        }
    }
  //  qDebug()<<lastx0;

    phase = posstim-lastx0+t;
    if (!hnt->inphase)
    {
        if (t>0)
            phase=posstim-lastx0;
        else if (t==0)
            phase=(hnt->srfr/hnt->osfr)/2+(posstim-lastx0);
    }
    hnt->phase=phase;
    fillstimforzc(posstim,hnt->stlength,phase);
}

void plotwindow::recurbutterf(int order, double sfr, double hpf, double* x, int length)
{
    // http://www.exstrom.com/journal/sigproc/ - recursive implementation of Butterworth filter
    int i;
    int n = order/2;
    double s = sfr;
    double f = hpf;
    double a = tan(M_PI*f/s);
    double a2 = a*a;
    double r;
    double *A = (double *)malloc(n*sizeof(double));
    double *d1 = (double *)malloc(n*sizeof(double));
    double *d2 = (double *)malloc(n*sizeof(double));
    double *w0 = (double *)calloc(n, sizeof(double));
    double *w1 = (double *)calloc(n, sizeof(double));
    double *w2 = (double *)calloc(n, sizeof(double));

   // int length=50;
   // double x[length];
   // for (int i=0; i<length; i++)
   // {
   //     x[i]=sin(9*2*PI/500*i);
   //     cout<<x[i]<<" ";
   // }
   // cout<<endl<<endl;

    for(i=0; i<n; i++)
    {
        r = sin(M_PI*(2.0*i+1.0)/(4.0*n));
        s = a2 + 2.0*a*r + 1.0;
        A[i] = a2/s;
        d1[i] = 2.0*(1-a2)/s;
        d2[i] = -(a2 - 2.0*a*r + 1.0)/s;   
    }

    for (int k=0; k<length; k++)
    {
        for(i=0; i<n; i++)
        {
            w0[i] = d1[i]*w1[i] + d2[i]*w2[i] + x[k];
            x[k] = A[i]*(w0[i] + 2.0*w1[i] + w2[i]);
            w2[i] = w1[i];
            w1[i] = w0[i];
        }     
    }  
}

void plotwindow::recurbutterfilter(int posstim, int length)
{
    double t1[length];
    for (int i=0; i<length; i++)
         t1[i]=arrc.amp1[posstim+i];
    recurbutterf(butterord,hnt->srfr,flowpass,t1,length);
    for (int i=20; i<length; i++)
         arrc.amp1[posstim+i]=t1[i];
}

void plotwindow::recurbutterfilterar(int posstim, int length)
{
    double t1[length];
    for (int i=0; i<length; i++)
         t1[i]=arrc.of2[posstim+i];
    recurbutterf(butterord,hnt->srfr,flowpass,t1,length);
    for (int i=20; i<length; i++)
         arrc.of2[posstim+i]=t1[i];

}

// ==== butterfiltcoeffs implementation ==== //

double* plotwindow::ComputeLP( int FilterOrder )
{
    double *NumCoeffs;
    int m;
    int i;

    NumCoeffs = (double *)calloc( FilterOrder+1, sizeof(double) );
    if( NumCoeffs == NULL ) return( NULL );

    NumCoeffs[0] = 1;
    NumCoeffs[1] = FilterOrder;
    m = FilterOrder/2;
    for( i=2; i <= m; ++i)
    {
        NumCoeffs[i] =(double) (FilterOrder-i+1)*NumCoeffs[i-1]/i;
        NumCoeffs[FilterOrder-i]= NumCoeffs[i];
    }
    NumCoeffs[FilterOrder-1] = FilterOrder;
    NumCoeffs[FilterOrder] = 1;

    return NumCoeffs;
}

double* plotwindow::ComputeHP( int FilterOrder )
{
    double *NumCoeffs;
    int i;

    NumCoeffs = ComputeLP(FilterOrder);
    if(NumCoeffs == NULL ) return( NULL );

    for( i = 0; i <= FilterOrder; ++i)
        if( i % 2 ) NumCoeffs[i] = -NumCoeffs[i];

    return NumCoeffs;
}

double* plotwindow::TrinomialMultiply( int FilterOrder, double *b, double *c )
{
    int i, j;
    double *RetVal;

    RetVal = (double *)calloc( 4 * FilterOrder, sizeof(double) );
    if( RetVal == NULL ) return( NULL );

    RetVal[2] = c[0];
    RetVal[3] = c[1];
    RetVal[0] = b[0];
    RetVal[1] = b[1];

    for( i = 1; i < FilterOrder; ++i )
    {
        RetVal[2*(2*i+1)]   += c[2*i] * RetVal[2*(2*i-1)]   - c[2*i+1] * RetVal[2*(2*i-1)+1];
        RetVal[2*(2*i+1)+1] += c[2*i] * RetVal[2*(2*i-1)+1] + c[2*i+1] * RetVal[2*(2*i-1)];

        for( j = 2*i; j > 1; --j )
        {
            RetVal[2*j]   += b[2*i] * RetVal[2*(j-1)]   - b[2*i+1] * RetVal[2*(j-1)+1] +
                c[2*i] * RetVal[2*(j-2)]   - c[2*i+1] * RetVal[2*(j-2)+1];
            RetVal[2*j+1] += b[2*i] * RetVal[2*(j-1)+1] + b[2*i+1] * RetVal[2*(j-1)] +
                c[2*i] * RetVal[2*(j-2)+1] + c[2*i+1] * RetVal[2*(j-2)];
        }

        RetVal[2] += b[2*i] * RetVal[0] - b[2*i+1] * RetVal[1] + c[2*i];
        RetVal[3] += b[2*i] * RetVal[1] + b[2*i+1] * RetVal[0] + c[2*i+1];
        RetVal[0] += b[2*i];
        RetVal[1] += b[2*i+1];
    }

    return RetVal;
}

double* plotwindow::ComputeNumCoeffs(int FilterOrder)
{
    double *TCoeffs;
    double *NumCoeffs;
    int i;

    NumCoeffs = (double *)calloc( 2*FilterOrder+1, sizeof(double) );
    if( NumCoeffs == NULL ) return( NULL );

    TCoeffs = ComputeHP(FilterOrder);
    if( TCoeffs == NULL ) return( NULL );

    for( i = 0; i < FilterOrder; ++i)
    {
        NumCoeffs[2*i] = TCoeffs[i];
        NumCoeffs[2*i+1] = 0.0;
    }
    NumCoeffs[2*FilterOrder] = TCoeffs[FilterOrder];

    free(TCoeffs);

    return NumCoeffs;
}

double* plotwindow::ComputeNumCoeffs(int FilterOrder, double Lcutoff, double Ucutoff, double *DenC)
{
    double *TCoeffs;
    double *NumCoeffs;
    std::complex<double> *NormalizedKernel;
    double Numbers[17]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16};
    int i;

    NumCoeffs = (double *)calloc( 2*FilterOrder+1, sizeof(double) );
    if( NumCoeffs == NULL ) return( NULL );

    TCoeffs = ComputeHP(FilterOrder);
    if( TCoeffs == NULL ) return( NULL );

    NormalizedKernel = (std::complex<double> *)calloc( 2*FilterOrder+1, sizeof(std::complex<double>) );
    if( NormalizedKernel == NULL ) return( NULL );

    for( i = 0; i < FilterOrder; ++i)
    {
        NumCoeffs[2*i] = TCoeffs[i];
        NumCoeffs[2*i+1] = 0.0;
    }
    NumCoeffs[2*FilterOrder] = TCoeffs[FilterOrder];
    double cp[2];
    double Wn;
    cp[0] = 2*2.0*tan(PI * Lcutoff / 2.0);
    cp[1] = 2*2.0*tan(PI * Ucutoff / 2.0);
    //double Bw;
   // Bw = cp[1] - cp[0];
    //center frequency
    Wn = sqrt(cp[0]*cp[1]);
    Wn = 2*atan2(Wn,4);
    const std::complex<double> result = std::complex<double>(-1,0);

    for(int k = 0; k<FilterOrder*2+1; k++)
    {
        NormalizedKernel[k] = std::exp(-sqrt(result)*Wn*Numbers[k]);
    }

    double b=0;
    double den=0;
    for(int d = 0; d<FilterOrder*2+1; d++)
    {
        b+=real(NormalizedKernel[d]*NumCoeffs[d]);
        den+=real(NormalizedKernel[d]*DenC[d]);
    }
    for(int c = 0; c<FilterOrder*2+1; c++)
    {
        NumCoeffs[c]=(NumCoeffs[c]*den)/b;
    }

    free(TCoeffs);
    return NumCoeffs;
}

double* plotwindow::ComputeDenCoeffs( int FilterOrder, double Lcutoff, double Ucutoff )
{
    int k;            // loop variables
    double theta;     // PI * (Ucutoff - Lcutoff) / 2.0
    double cp;        // cosine of phi
    double st;        // sine of theta
    double ct;        // cosine of theta
    double s2t;       // sine of 2*theta
    double c2t;       // cosine 0f 2*theta
    double *RCoeffs;     // z^-2 coefficients
    double *TCoeffs;     // z^-1 coefficients
    double *DenomCoeffs;     // dk coefficients
    double PoleAngle;      // pole angle
    double SinPoleAngle;     // sine of pole angle
    double CosPoleAngle;     // cosine of pole angle
    double a;         // workspace variables

    cp = cos(PI * (Ucutoff + Lcutoff) / 2.0);
    theta = PI * (Ucutoff - Lcutoff) / 2.0;
    st = sin(theta);
    ct = cos(theta);
    s2t = 2.0*st*ct;        // sine of 2*theta
    c2t = 2.0*ct*ct - 1.0;  // cosine of 2*theta

    RCoeffs = (double *)calloc( 2 * FilterOrder, sizeof(double) );
    TCoeffs = (double *)calloc( 2 * FilterOrder, sizeof(double) );

    for( k = 0; k < FilterOrder; ++k )
    {
        PoleAngle = PI * (double)(2*k+1)/(double)(2*FilterOrder);
        SinPoleAngle = sin(PoleAngle);
        CosPoleAngle = cos(PoleAngle);
        a = 1.0 + s2t*SinPoleAngle;
        RCoeffs[2*k] = c2t/a;
        RCoeffs[2*k+1] = s2t*CosPoleAngle/a;
        TCoeffs[2*k] = -2.0*cp*(ct+st*SinPoleAngle)/a;
        TCoeffs[2*k+1] = -2.0*cp*st*CosPoleAngle/a;
    }

    DenomCoeffs = TrinomialMultiply(FilterOrder, TCoeffs, RCoeffs );
    free(TCoeffs);
    free(RCoeffs);

    DenomCoeffs[1] = DenomCoeffs[0];
    DenomCoeffs[0] = 1.0;
    for( k = 3; k <= 2*FilterOrder; ++k )
        DenomCoeffs[k] = DenomCoeffs[2*k-2];

    return DenomCoeffs;
}

void plotwindow::filter(int ord, double *a, double *b, int np, double *x, double *y)
{
    int i,j;
    y[0]=b[0] * x[0];
    for (i=1;i<ord+1;i++)
    {
        y[i]=0.0;
        for (j=0;j<i+1;j++)
            y[i]=y[i]+b[j]*x[i-j];
        for (j=0;j<i;j++)
            y[i]=y[i]-a[j+1]*y[i-j-1];
    }
    for (i=ord+1;i<np+1;i++)
    {
        y[i]=0.0;
        for (j=0;j<ord+1;j++)
            y[i]=y[i]+b[j]*x[i-j];
        for (j=0;j<ord;j++)
            y[i]=y[i]-a[j+1]*y[i-j-1];
    }
}

void plotwindow::butterfiltcoefs(int lcut, int hcut, int order, int sampr)
{
    // https://stackoverflow.com/questions/10373184/bandpass-butterworth-filter-implementation-in-c
    //Frequency bands is a vector of values - Lower Frequency Band and Higher Frequency Band
    //First value is lower cutoff and second value is higher cutoff
    // 10 Hz = 10/(500/2), 50 Hz = 50/(500/2)

    double lc=(double)lcut/(sampr/2);
    double hc=(double)hcut/(sampr/2);
    double FrequencyBands[2] = {lc,hc};   // these values are as a ratio of f/fs, where fs is sampling rate, and f is cutoff frequency
                                             // and therefore should lie in the range [0 1]

    DenC = ComputeDenCoeffs(order, FrequencyBands[0], FrequencyBands[1]); // is A in matlab function
  /*  cout<<"DenC is: ";
    for (int k = 0; k<butterord*2+1; k++)
        cout<<DenC[k]<<" ";    
    cout<<endl; */

    NumC = ComputeNumCoeffs(order,lc,hc,DenC); // is B in matlab
  /*  cout<<"NumC is: ";
    for (int k = 0; k<butterord*2+1; k++)
        cout<<NumC[k]<<" ";    
    cout<<endl; */

    double y[order];
    double x[16]={1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16};
    filter(order, DenC, NumC, order, x, y);
}

// = end of butterfiltcoeffs implementation = //

// ==== filtfilt implementation ==== //
// https://stackoverflow.com/questions/17675053/matlabs-filtfilt-algorithm

void plotwindow::add_index_range(vectori &indices, int beg, int end)
{
    int inc=1;
    for (int i = beg; i <= end; i += inc)
       indices.push_back(i);
}

void plotwindow::add_index_const(vectori &indices, int value, size_t numel)
{
    while (numel--)
        indices.push_back(value);
}

void plotwindow::append_vector(vectord &vec, const vectord &tail)
{
    vec.insert(vec.end(), tail.begin(), tail.end());
}

vectord plotwindow::subvector_reverse(const vectord &vec, int idx_end, int idx_start)
{
    vectord result(&vec[idx_start], &vec[idx_end+1]);
    std::reverse(result.begin(), result.end());
    return result;
}

inline int max_val(const vectori& vec)
{
    return std::max_element(vec.begin(), vec.end())[0];
}

void plotwindow::filter(vectord B, vectord A, const vectord &X, vectord &Y, vectord &Zi)
{
    if (A.empty())
        throw std::domain_error("The feedback filter coefficients are empty.");
    if (std::all_of(A.begin(), A.end(), [](double coef){ return coef == 0; }))
        throw std::domain_error("At least one of the feedback filter coefficients has to be non-zero.");
    if (A[0] == 0)
        throw std::domain_error("First feedback coefficient has to be non-zero.");

    // Normalize feedback coefficients if a[0] != 1;
    auto a0 = A[0];
    if (a0 != 1.0)
    {
        std::transform(A.begin(), A.end(), A.begin(), [a0](double v) { return v / a0; });
        std::transform(B.begin(), B.end(), B.begin(), [a0](double v) { return v / a0; });
    }

    size_t input_size = X.size();
    size_t filter_order = std::max(A.size(), B.size());
    B.resize(filter_order, 0);
    A.resize(filter_order, 0);
    Zi.resize(filter_order, 0);
    Y.resize(input_size);

    const double *x = &X[0];
    const double *b = &B[0];
    const double *a = &A[0];
    double *z = &Zi[0];
    double *y = &Y[0];

    for (size_t i = 0; i < input_size; ++i)
    {
        size_t order = filter_order - 1;
        while (order)
        {
            if (i >= order)
                z[order - 1] = b[order] * x[i - order] - a[order] * y[i - order] + z[order];
            --order;
        }
        y[i] = b[0] * x[i] + z[0];
    }
    Zi.resize(filter_order - 1);
}

void plotwindow::filtfilt(vectord B, vectord A, const vectord &X, vectord &Y)
{

    int len = X.size();     // length of input
    int na = A.size();
    int nb = B.size();
    int nfilt = (nb > na) ? nb : na;
    int nfact = 3 * (nfilt - 1); // length of edge transients

    if (len <= nfact)
        throw std::domain_error("Input data too short! Data must have length more than 3 times filter order.");

    // set up filter's initial conditions to remove DC offset problems at the
    // beginning and end of the sequence
    B.resize(nfilt, 0);
    A.resize(nfilt, 0);

    vectori rows, cols;
    //rows = [1:nfilt-1           2:nfilt-1             1:nfilt-2];
    add_index_range(rows, 0, nfilt - 2);
    if (nfilt > 2)
    {
        add_index_range(rows, 1, nfilt - 2);
        add_index_range(rows, 0, nfilt - 3);
    }
    //cols = [ones(1,nfilt-1)         2:nfilt-1          2:nfilt-1];
    add_index_const(cols, 0, nfilt - 1);
    if (nfilt > 2)
    {
        add_index_range(cols, 1, nfilt - 2);
        add_index_range(cols, 1, nfilt - 2);
    }
    // data = [1+a(2)         a(3:nfilt)        ones(1,nfilt-2)    -ones(1,nfilt-2)];

    auto klen = rows.size();
    vectord data;
    data.resize(klen);
    data[0] = 1 + A[1];  int j = 1;
    if (nfilt > 2)
    {
        for (int i = 2; i < nfilt; i++)
            data[j++] = A[i];
        for (int i = 0; i < nfilt - 2; i++)
            data[j++] = 1.0;
        for (int i = 0; i < nfilt - 2; i++)
            data[j++] = -1.0;
    }

    vectord leftpad = subvector_reverse(X, nfact, 1);
    double _2x0 = 2 * X[0];
    std::transform(leftpad.begin(), leftpad.end(), leftpad.begin(), [_2x0](double val) {return _2x0 - val; });

    vectord rightpad = subvector_reverse(X, len - 2, len - nfact - 1);
    double _2xl = 2 * X[len-1];
    std::transform(rightpad.begin(), rightpad.end(), rightpad.begin(), [_2xl](double val) {return _2xl - val; });

    double y0;
    vectord signal1, signal2, zi;

    signal1.reserve(leftpad.size() + X.size() + rightpad.size());
    append_vector(signal1, leftpad);
    append_vector(signal1, X);
    append_vector(signal1, rightpad);

    // Calculate initial conditions
    MatrixXd sp = MatrixXd::Zero(max_val(rows) + 1, max_val(cols) + 1);
    for (size_t k = 0; k < klen; ++k)
    {
        sp(rows[k], cols[k]) = data[k];
    }
    auto bb = VectorXd::Map(B.data(), B.size());
    auto aa = VectorXd::Map(A.data(), A.size());
    MatrixXd zzi = (sp.inverse() * (bb.segment(1, nfilt - 1) - (bb(0) * aa.segment(1, nfilt - 1))));
    zi.resize(zzi.size());

    // Do the forward and backward filtering
    y0 = signal1[0];
    std::transform(zzi.data(), zzi.data() + zzi.size(), zi.begin(), [y0](double val){ return val*y0; });
    filter(B, A, signal1, signal2, zi);
    std::reverse(signal2.begin(), signal2.end());
    y0 = signal2[0];
    std::transform(zzi.data(), zzi.data() + zzi.size(), zi.begin(), [y0](double val){ return val*y0; });
    filter(B, A, signal2, signal1, zi);
    Y = subvector_reverse(signal1, signal1.size() - nfact - 1, nfact);
}

// = end of filtfilt implementation = //

void plotwindow::zerophaseinit(int lcut, int hcut, int order, int sampr)
{
    butterfiltcoefs(lcut, hcut, order, sampr);
    acoeffs = vector<double>(order*2+1);
    for(int k = 0; k<order*2+1; k++)
    {
    //   cout<<DenC[k]<<" ";
        acoeffs[k]=DenC[k];
    }
    // cout<<endl;
    bcoeffs = vector<double>(order*2+1);
    for(int k = 0; k<order*2+1; k++)
    {
    //    cout<<NumC[k]<<" ";
        bcoeffs[k]=NumC[k];
    }
    // cout<<endl;
}

void plotwindow::zerophasefilt(int posstim, int length)
{
    vectord input_signal = vector<double>(length);
    double mean=0;
    for(int k = 0; k<length; k++)
    {
        input_signal[k]=arrc.amp1[posstim+k];
     //   mean+=input_signal[k];
    }
  //  mean/=length;
    vectord y_filter_out, y_filtfilt_out; vectord zi = {0};
    filter(bcoeffs, acoeffs, input_signal, y_filter_out, zi);
    filtfilt(bcoeffs, acoeffs, input_signal, y_filtfilt_out);
    for(int k = 0; k<length; k++)
        arrc.amp1[posstim+k]=y_filtfilt_out[k]*5;
}

void plotwindow::zerophasefiltar(int posstim, int length)
{
    vectord input_signal = vector<double>(length);
    double mean=0;
    for(int k = 0; k<length; k++)
    {
        input_signal[k]=arrc.of2[posstim+k];
      //  mean+=input_signal[k];
    }
 //   mean/=length;
    vectord y_filter_out, y_filtfilt_out; vectord zi = {0};
    filter(bcoeffs, acoeffs, input_signal, y_filter_out, zi);
    filtfilt(bcoeffs, acoeffs, input_signal, y_filtfilt_out);
    for(int k = 0; k<length; k++)
        arrc.of2[posstim+k]=(int)mean+y_filtfilt_out[k];
}

void plotwindow::zerophasefiltzc(int posstim, int length)
{
    if (!usefiltering)
    {
        vectord input_signal = vector<double>(length);
        for(int k = 0; k<length; k++)
            input_signal[k]=arrc.amp2[posstim+k];
        vectord y_filter_out, y_filtfilt_out; vectord zi = {0};
        filter(bcoeffs, acoeffs, input_signal, y_filter_out, zi);
        filtfilt(bcoeffs, acoeffs, input_signal, y_filtfilt_out);
        for(int k = 0; k<length; k++)
            arrc.amp2[posstim+k]=y_filtfilt_out[k];
    } else
    for(int k = 0; k<length; k++)
        arrc.amp2[posstim+k]=arrc.amp1[posstim+k];
}

void plotwindow::zerophasefiltssd(int posstim, int length)
{
    vectord input_signal = vector<double>(length);
    double mean=0;
    for(int k = 0; k<length; k++)
    {
        input_signal[k]=arrc.amp2[posstim+k];
        mean+=input_signal[k];
    }
    mean/=length;
    vectord y_filter_out, y_filtfilt_out; vectord zi = {0};
    filter(bcoeffs, acoeffs, input_signal, y_filter_out, zi);
    filtfilt(bcoeffs, acoeffs, input_signal, y_filtfilt_out);
    for(int k = 0; k<length; k++)
        arrc.amp2[posstim+k]=(int)mean+y_filtfilt_out[k];
}

void plotwindow::testfiltfilt()
{
    int length=samplesize;
    double* input = new double[length];
    double* filteredd = new double[length];

    QFile inputFile("C:/EEGdat/inp1.dat");
    inputFile.open(QIODevice::ReadOnly);
    QTextStream fin(&inputFile);
    QString line;
    QStringList sl;
    line = fin.readLine();
    sl = line.split(QRegExp("\\,"), QString::SkipEmptyParts);
    for (int k=0; k < length; k++)
        input[k]=sl[k].toDouble();
    inputFile.close();

    zerophaseinit(8, 12, 4, 200);
    vectord inps = vector<double>(length);
    for(int k = 0; k<length; k++)
        inps[k]=input[k];

    vectord y_filter_out, y_filtfilt_out; vectord zi = {0};
    filter(bcoeffs, acoeffs, inps, y_filter_out, zi);
    filtfilt(bcoeffs, acoeffs, inps, y_filtfilt_out);
    for(int k = 0; k<length; k++)
        filteredd[k]=y_filtfilt_out[k];

    QFile inputFil("C:/EEGdat/filtinp1.dat");
    inputFil.open(QIODevice::ReadOnly);
    QTextStream fi(&inputFil);
    line = fi.readLine();
    sl = line.split(QRegExp("\\,"), QString::SkipEmptyParts);
    for (int k=0; k < length; k++)
        input[k]=sl[k].toDouble();
    inputFil.close();

    double diff=0;
    for(int k = 0; k<length; k++)
        diff+=abs(filteredd[k]-input[k]);
    cout<<diff<<endl;
}

void plotwindow::getcoeffs(vectord& acc, vectord& bcc, int ord)
{
    for(int k = 0; k<ord*2+1; k++)
    {
        acc[k]=acoeffs[k];
        bcc[k]=bcoeffs[k];
    }
 //   for(int k = 0; k<ord*2+1; k++)
 //       cout<<acc[k]<<" ";
 //   cout<<endl;
//    for(int k = 0; k<ord*2+1; k++)
 //       cout<<bcc[k]<<" ";
//    cout<<endl;
//    cout<<endl;
}
MatrixXf plotwindow::covarian(MatrixXf mat)
{
    MatrixXf centered = mat.rowwise() - mat.colwise().mean();
    MatrixXf cov = (centered.adjoint() * centered) / double(mat.rows() - 1);
    return cov;
}

MatrixXf plotwindow::readdata(int snums, int ssize, QString fname)
{
    MatrixXf inp(snums,ssize); // points, channels
    QFile inputFile(fname); // C:/EEGdat/alexeeg.dat");
    inputFile.open(QIODevice::ReadOnly);
    QTextStream fin(&inputFile);
    QString line;
    QStringList sl;
    for (int j=0; j<ssize; j++)
    {
        line = fin.readLine();
        sl = line.split(QRegExp("\\s+"), QString::SkipEmptyParts);
        for (int i=0; i < snums; i++)
            inp(i,j)=sl[i].toDouble();
    }
    inputFile.close();
    return inp;
}

MatrixXf plotwindow::readdatafromeeg(int t)
{
    samplenums=t;
    samplesize=chnums;
    MatrixXf inp(samplenums,samplesize);
    for (int j=0; j<samplenums; j++)
        for (int i=0; i < samplesize; i++)
            inp(j,i)=rawdata[i][j];
    return inp;
}

void plotwindow::savessdtofile(MatrixXf xssd, int dim, QString fname)
{
    QFile outputFile(fname);
    outputFile.open(QIODevice::WriteOnly);
    QTextStream fout(&outputFile);
    for (int i=0; i<samplesize; i++)
    {
        for (int j=0; j<dim; j++)
            fout << xssd(i,j) << ",";
        fout<<endl;
    }
    outputFile.close();
}

MatrixXf plotwindow::filtfiltmatr(vectord b, vectord a, MatrixXf inp, int snums, int ssize)
{
    MatrixXf out(snums,ssize);
    vectord inpv = vector<double>(snums);
    for (int i=0; i<ssize; i++)
    {        
        for(int j = 0; j<snums; j++)
            inpv[j]=inp(j,i);        
        vectord y_filter_out, y_filtfilt_out; vectord zi = {0};
        filter(b, a, inpv, y_filter_out, zi);
        filtfilt(b, a, inpv, y_filtfilt_out);
        for(int j = 0; j<snums; j++)
            out(j,i)=y_filtfilt_out[j];
    } 
    return out;
}

MatrixXf plotwindow::filtfiltmatrbytrials(vectord b, vectord a, MatrixXf inp, int trialnums, int ssize)
{
    int trl = 750; // length of trial
    MatrixXf out(trialnums*trl,ssize);
    vectord inpv = vector<double>(trl);
    for (int i=0; i<ssize; i++)
    {
        for (int p = 0; p<trialnums; p++)
        {
            for (int j = p*trl; j<(p+1)*trl; j++)
                inpv[j-p*trl]=inp(j,i);
            vectord y_filter_out, y_filtfilt_out; vectord zi = {0};
            filter(b, a, inpv, y_filter_out, zi);
            filtfilt(b, a, inpv, y_filtfilt_out);
            for(int j = p*trl; j<(p+1)*trl; j++)
                out(j,i)=y_filtfilt_out[j-p*trl];
        }
    }
    return out;
}

void plotwindow::ssd(QString fname, int trialnum)
{
    int ord=butterord; int srf=500;//hnt->srfr;

    bc = vector<double>(ord*2+1);
    ac = vector<double>(ord*2+1);
    b_f = vector<double>(ord*2+1);
    a_f = vector<double>(ord*2+1);
    b_s = vector<double>(ord*2+1);
    a_s = vector<double>(ord*2+1);

  //  timer.restart();
    //noibs.lf=8; noibs.hf=13;

    zerophaseinit(sigb.lf,sigb.hf,ord,srf);
    getcoeffs(ac,bc,ord);
    zerophaseinit(noibp.lf,noibp.hf,ord,srf);
    getcoeffs(a_f,b_f,ord);
    //zerophaseinit(sigb.lf,sigb.hf,ord,srf);
    //getcoeffs(a_s,b_s,ord); 

    MatrixXf xinp;

    if (extractingssd)
        xinp = readdatafromeeg(hnt->srfr*ssdtime);
    else
        xinp = inpssd;  //readdata(samplenums, samplesize, fname); // offline case

    MatrixXf xfilt = filtfiltmatrbytrials(bc,ac,xinp,trialnum,samplesize);
    MatrixXf C_s = covarian(xfilt);

  //  cout << C_s;

    MatrixXf xflank = filtfiltmatrbytrials(b_f,a_f,xinp,trialnum,samplesize);
 // MatrixXf xfl = filtfiltmatr(b_s,a_s,xflank); // for case if have bandstop koeffs
    MatrixXf xfl = xflank - xfilt;
    MatrixXf C_n = covarian(xfl);

  //  cout << C_n;

    SelfAdjointEigenSolver<MatrixXf> es(C_s);
    MatrixXf V = es.eigenvectors();
    VectorXf eigv = es.eigenvalues();    

    std::sort(eigv.data(),eigv.data()+eigv.size(),std::greater<float>());
    MatrixXf rV = V.rowwise().reverse();

    double tol = eigv[0]*1e-10; int c=0;
    for (int i=0; i<samplesize; i++)
        if (eigv[i]>tol)
            c++;
    MatrixXf M(samplesize,c);
    if (c<samplesize)
    {
        cout<<"Input data does not have full rank. Only "<<c<<" components can be computed."<<endl<<endl;
        for (int i=0; i<c; i++)
            eigv[i]=1/qSqrt(eigv[i]);
        for (int i=0; i<samplesize; i++)
            for (int j=0; j<c; j++)
             M(i,j)=eigv[j]*rV(i,j);
        maxssdcomp=c;
    } else
    {
        DiagonalMatrix<float, Dynamic> Dl(samplesize);
        M=Dl;
        maxssdcomp=samplesize;
    }        

    MatrixXf C_s_r =  M.transpose() * C_s * M;  // different order with matlab, here ~ e-006, in matlab ~ e-015, but C_s = C_s (matlab)
    MatrixXf C_n_r =  M.transpose() * C_n * M;  // equal values with matlab, but different signs

  //  cout << C_s_r << endl << endl;
  //  cout << C_n_r << endl << endl;

    GeneralizedSelfAdjointEigenSolver<MatrixXf> ges(C_s_r,C_s_r+C_n_r);
    MatrixXf W = ges.eigenvectors();
   // cout << W << endl << endl;

    eigv = ges.eigenvalues();
 //   cout << eigv << endl << endl;

    std::sort(eigv.data(),eigv.data()+eigv.size(),std::greater<float>());
    MatrixXf rW = W.rowwise().reverse();  

    wW = M * rW;                   // the de-mixing matrix. Each column is a spatial filter and the
                                      // timecourse of the SSD components is extracted with X * W
   // cout<< wW << endl << endl;
  /*  MatrixXf t = wW.transpose() * C_s * wW;
    MatrixXf A = C_s * wW * t.inverse();    //  the spatial patterns (also called mixing matrix) with the i'th column
  */                                          //  corrresponding to the i'th SSD component
    MatrixXf invF = wW;
    QFile outputFile("C:/EEGdat/analysis/SSD_subj"+QString::number(currentsubj)+".dat");
    outputFile.open(QIODevice::WriteOnly);
    QTextStream fout(&outputFile);
    for (int i=0; i<invF.rows(); i++)
    {
        for (int j=0; j<invF.cols(); j++)
        {
            fout<<invF(i,j);
            if (j<invF.cols()-1)
                fout<<",";
        }
        fout<<endl;
    }
    outputFile.close();

   //   Xssd = xfilt * wW;

  // cout << "SSD time calculation (msec): " << timer.nsecsElapsed() / 1000000.0 << endl;

  //  cout << Xssd;

   // savessdtofile(wW,c,"C:/EEGdat/Wssdonline.dat");

  //  MatrixXf rXssd = Xssd.transpose();
  //  samplenums=1;
  //  savessdtofile(rXssd,1000,"C:/EEGdat/C1c++.dat");
}

void plotwindow::getssd(int comp, int length, QString fname, int blnum)
{
    if (applyssd)
    {
        MatrixXf xinp = readdatafromeeg(hnt->imlength);
        MatrixXf xfilt = filtfiltmatrbytrials(bc,ac,xinp,length/500,samplesize);
        MatrixXf Xres = xfilt * wW;
        // TO DO: Xres(i,t) -> arrc.amp1(i)
    }
    else
    {  
        spatial_fitering=true;

        QFile inputFile(fname);
        inputFile.open(QIODevice::ReadOnly);
        QTextStream fin(&inputFile);
        int lc=0;
        while (!fin.atEnd())
        {
            fin.readLine();
            lc++;
        }
        inputFile.close();

        samplesize=lc;

        int lenint = 25000*blnum; int chn=23; int int_len=500;
        MatrixXf xpre(lenint,samplesize);
        MatrixXf xpost(lenint,samplesize);

        QFile inpFile(fname);
        inpFile.open(QIODevice::ReadOnly);
        QTextStream finp(&inpFile);
        QString line;
        QStringList sl;

        line = finp.readLine();
        sl = line.split(QRegExp("\\s+"), QString::SkipEmptyParts);
        int length = sl.length();

        MatrixXf xinp(length,samplesize); // points, channels

        for (int i=0; i < length; i++)
            xinp(i,0)=sl[i].toDouble();

        for (int j=1; j<samplesize; j++)
        {
            line = finp.readLine();
            sl = line.split(QRegExp("\\s+"), QString::SkipEmptyParts);
            for (int i=0; i < length; i++)
                xinp(i,j)=sl[i].toDouble();
        }

        int lt = sl.length();
        int allst=blnum*50; int bordvalue=1000; int jump=1250;
        int c=0; int stnum=0; int sh1=20; int sh2=40;

        while ((c<lt) && (stnum<allst))  // extraction of pre / post data for all channels
        {
            c++;
            if (abs(xinp(c,chn))>=bordvalue)
            {
                for (int t=0; t<samplesize; t++)
                {
                    for (int i=c-int_len-sh1; i<c-sh1; i++)
                        xpre(stnum*int_len+i-c+int_len+sh1,t)=xinp(i,t);
                    for (int i=c+int_len+sh2; i<c+int_len*2+sh2; i++)
                        xpost(stnum*int_len-c-int_len-sh2+i,t)=xinp(i,t);
                }
                stnum++;
                c+=jump;
            }
        }   

        int inittrials=100;
        samplenums=2*inittrials*500;                     // 2 blocks (EC+EO), in each (50+50)*500
        MatrixXf inssd(samplenums,samplesize);

        for (int t=0; t<samplesize; t++)
        {
            for (int p=0; p<inittrials; p++)
            for (int i=0; i<int_len; i++)
            {
                inssd(int_len*(p*2)+i,t)=xpre(int_len*p+i);
                inssd(int_len*(p*2+1)+i,t)=xpost(int_len*p+i);
            }
        }

        inpssd=inssd;        
        ssd("",2*inittrials);

        xpre = filtfiltmatrbytrials(bc,ac,xpre,allst,samplesize);    // band-pass filtering
        xpost = filtfiltmatrbytrials(bc,ac,xpost,allst,samplesize);

        MatrixXf Xrespre = xpre * wW; // Xres: length:K, xfilt: length:31, wW: 31:K
        MatrixXf Xrespost = xpost * wW;

        for (int p=0; p<blnum; p++)
        for (int i=0; i<50; i++)
        {
            for (int j=0; j<int_len; j++)
            {
                arrc.amp1[blockstarts[p]+i*1500+j]=Xrespre(p*25000+i*500+j,comp-1)/5000000;
                arrc.amp1[blockstarts[p]+i*1500+1000+j]=Xrespost(p*25000+i*500+j,comp-1)/5000000;
            }
        }
        ui->widget->graph(0)->setData(arrc.xc, arrc.amp1);
        ui->widget->replot();
    }
}

VectorXf plotwindow::getmpart(VectorXf inp, int p)
{
    int int_len=500;
    VectorXf out(int_len);
    for (int i=(p-1)*int_len; i<p*int_len; i++)
        out(i-(p-1)*int_len)=inp(i);
    return out;
}

void plotwindow::detrend(double *y, double *x, int m)
{
    double xmean, ymean;
    int i;
    double temp, Sxy, Sxx, grad, yint;
    /********************************
       Set the X axis Liner Values
       *********************************/
    for (i = 0; i < m; i++)
        x[i] = i;

    /********************************
       Calculate the mean of x and y
       *********************************/
    xmean = 0;
    ymean = 0;
    for (i = 0; i < m; i++)
    {
        xmean += x[i];
        ymean += y[i];
    }
    xmean /= m;
    ymean /= m;

    /********************************
       Calculate Covariance
       *********************************/
    temp = 0;
    for (i = 0; i < m; i++)
        temp += x[i] * y[i];
    Sxy = temp / m - xmean * ymean;

    temp = 0;
    for (i = 0; i < m; i++)
        temp += x[i] * x[i];
    Sxx = temp / m - xmean * xmean;

    /********************************
       Calculate Gradient and Y intercept
       *********************************/
    grad = Sxy / Sxx;
    yint = -grad * xmean + ymean;

    /********************************
       Removing Linear Trend
       *********************************/
    for (i = 0; i < m; i++)
        y[i] = y[i] - (grad * i + yint);
}

void plotwindow::flipsign()
{
    QFile inresFile("E:/EEGdat/analysis/CSPtopo/o800_subj18.dat");
    inresFile.open(QIODevice::ReadOnly);
    QTextStream fresin(&inresFile);

    QFile outresFile("C:/EEGdat/analysis/CSP_topo800/CSPtopo800_subj18_.dat");
    outresFile.open(QIODevice::WriteOnly);
    QTextStream fresout(&outresFile);
    QStringList sl;
    int numbl;
    QString st;

    for (int i=0; i<25; i++)
    {
        st = fresin.readLine();
        sl = st.split(QRegExp("\\,"), QString::SkipEmptyParts);
        numbl = sl.length();
        for (int j=0; j<numbl; j++)
        if (j<numbl-1)
            fresout<<-1*sl[j].toDouble()<<",";
        else
            fresout<<-1*sl[j].toDouble();
        fresout<<endl;
    }
    inresFile.close();
    outresFile.close();
}


void plotwindow::sortreferencelist()
{
    QFile inputFile("D:/Androxim/refers.txt");
    inputFile.open(QIODevice::ReadOnly);
    QTextStream fin(&inputFile);
    QStringList sl;
    QString line;

    for (int i=0; i<68; i++)
    {
        line = fin.readLine();
        sl.append(line);
    }
    //inputFile.close();

    qSort(sl.begin(),sl.end());

    QFile outputFile("D:/Androxim/refsort.txt");
    outputFile.open(QIODevice::Append);
    QTextStream fout(&outputFile);

    for ( QStringList::Iterator it = sl.begin(); it != sl.end(); ++it )
                    fout << *it << "\n";

    outputFile.close();
    //for (int i=0; i<182; i++)
//        cout<<sl[i].<<endl;

}

MatrixXf plotwindow::getpart(MatrixXf inp, int p)
{
   // int bllen=25000;
    int int_len=750;
    MatrixXf out(int_len,samplesize);
    for (int t=0; t<samplesize; t++)
        for (int i=(p-1)*int_len; i<p*int_len; i++)
            out(i-(p-1)*int_len,t)=inp(i,t);
    return out;
}

MatrixXf plotwindow::cutpart(MatrixXf inp, int p, int trnum, int int_len)
{        
    MatrixXf out((trnum-1)*int_len,samplesize);
    for (int t=0; t<samplesize; t++)
        for (int i=0; i<(p-1)*int_len; i++)
            out(i,t)=inp(i,t);
    if (p!=trnum)
    {
        for (int t=0; t<samplesize; t++)
            for (int i = p*int_len; i<trnum*int_len; i++)
                out(i-int_len,t)=inp(i,t);
    }
    return out;
}

void plotwindow::fill_result_in_vector(VectorXf Xcsp_C1, VectorXf Xcsp_C2, int p)
{
    int int_len=750;
    for (int j=0; j<int_len; j++)
    {
        arrc.amp2[2*p*int_len+j]=Xcsp_C1(j);
        arrc.amp2[2*(p+1)*int_len-int_len+j]=Xcsp_C2(j);
    }
}

void plotwindow::meandatafromopregion(QString fname, int blnum, int currsubj)
{
    QFile inputFile(fname);                 // determine number of channels
    inputFile.open(QIODevice::ReadOnly);
    QTextStream fin(&inputFile);
    int lc=0;
    while (!fin.atEnd())
    {
        fin.readLine();
        lc++;
    }
    inputFile.close();
    samplesize=lc;

    int bllen = 25000; int chn=23; int int_len=500;
    int lenint = bllen*blnum;
    MatrixXf xpre(lenint,samplesize);
    MatrixXf xpost(lenint,samplesize);
    VectorXf meanxpre(lenint);
    VectorXf meanxpost(lenint);

    QFile inpFile(fname);                                        // read data from file
    inpFile.open(QIODevice::ReadOnly);
    QTextStream finp(&inpFile);
    QString line; QStringList sl;
    line = finp.readLine();
    sl = line.split(QRegExp("\\s+"), QString::SkipEmptyParts);
    int length = sl.length();
    MatrixXf xinp(length,samplesize);                    // points, channels
    for (int i=0; i < length; i++)
        xinp(i,0)=sl[i].toDouble();
    for (int j=1; j<samplesize; j++)
    {
        line = finp.readLine();
        sl = line.split(QRegExp("\\s+"), QString::SkipEmptyParts);
        for (int i=0; i < length; i++)
            xinp(i,j)=sl[i].toDouble();
    }
    inpFile.close();

    int lt = sl.length();                                       // extraction of pre / post data for all channels
    int allst=blnum*50; int bordvalue=1000; int jump=1250;
    int c=0; int stnum=0; int sh1=20; int sh2=40;
    while ((c<lt) && (stnum<allst))
    {
        c++;
        if (abs(xinp(c,chn))>=bordvalue)
        {
            for (int t=0; t<samplesize; t++)
            {
                for (int i=c-int_len-sh1; i<c-sh1; i++)
                    xpre(stnum*int_len+i-c+int_len+sh1,t)=xinp(i,t);
                for (int i=c+int_len+sh2; i<c+int_len*2+sh2; i++)
                    xpost(stnum*int_len-c-int_len-sh2+i,t)=xinp(i,t);
            }
            stnum++;
            c+=jump;
        }
    }

   lcutoff = 5; hcutoff=41;
   zerophaseinit(lcutoff, hcutoff, butterord, hnt->srfr);
   getcoeffs(ac,bc,butterord);
   xpre = filtfiltmatrbytrials(bc,ac,xpre,allst,samplesize); // band-pass filtering
   xpost = filtfiltmatrbytrials(bc,ac,xpost,allst,samplesize);


   int scale=10;

   double* valarrpre = new double[12];
   double* valarrpost = new double[12];
   for (int i=0; i<12; i++)
   {
       valarrpre[i] = 0;
       valarrpost[i] = 0;
   }

   double meanvpre;
   double meanvpost;

   for (int i=0; i<lenint; i++)
   {
       meanvpre = 0;
       meanvpost = 0;
       for (int j=0; j<po_list[currsubj]; j++) // j<po_list[currsubj]
       {
           meanvpre+=abs(xpre(i,parietal_occip_list[currsubj][j]));
           meanvpost+=abs(xpost(i,parietal_occip_list[currsubj][j]));
           valarrpre[j]=abs(xpre(i,parietal_occip_list[currsubj][j]));
           valarrpost[j]=abs(xpost(i,parietal_occip_list[currsubj][j]));
       }
       meanvpre/=po_list[currsubj];  // po_list[currsubj];
       meanvpost/=po_list[currsubj]; // po_list[currsubj];
       double stdvf = getstdvalue(valarrpre,po_list[currsubj]);
       meanxpre[i]=meanvpre/stdvf*scale;
       stdvf = getstdvalue(valarrpost,po_list[currsubj]);
       meanxpost[i]=meanvpost/stdvf*scale;
   }

   VectorXf xpre_data(int_len);
   VectorXf xpost_data(int_len);

   for (int i=1; i<=blnum*50; i++)
   {
       xpre_data = getmpart(meanxpre,i);
       xpost_data = getmpart(meanxpost,i);

       fill_result_in_vector(xpre_data, xpost_data, i);
   }

   ui->widget->graph(0)->setData(arrc.xc, arrc.amp1);
   ui->widget->replot();
   spatial_fitering = true;
}

void plotwindow::meandatafromopregionDT(QString fname, int blnum, int currsubj)
{
    QFile inputFile(fname);                 // determine number of channels
    inputFile.open(QIODevice::ReadOnly);
    QTextStream fin(&inputFile);
    int lc=0;
    while (!fin.atEnd())
    {
        fin.readLine();
        lc++;
    }
    inputFile.close();
    samplesize=lc;

    int bllen = 25000; int chn=23; int int_len=500;
    int lenint = bllen*blnum;
    MatrixXf xpre(lenint,samplesize);
    MatrixXf xpost(lenint,samplesize);
    MatrixXf xpre_data(int_len,samplesize);
    MatrixXf xpost_data(int_len,samplesize);

    QFile inpFile(fname);                                        // read data from file
    inpFile.open(QIODevice::ReadOnly);
    QTextStream finp(&inpFile);
    QString line; QStringList sl;
    line = finp.readLine();
    sl = line.split(QRegExp("\\s+"), QString::SkipEmptyParts);
    int length = sl.length();
    MatrixXf xinp(length,samplesize);                    // points, channels
    for (int i=0; i < length; i++)
        xinp(i,0)=sl[i].toDouble();
    for (int j=1; j<samplesize; j++)
    {
        line = finp.readLine();
        sl = line.split(QRegExp("\\s+"), QString::SkipEmptyParts);
        for (int i=0; i < length; i++)
            xinp(i,j)=sl[i].toDouble();
    }
    inpFile.close();

    int lt = sl.length();                                       // extraction of pre / post data for all channels
    int allst=blnum*50; int bordvalue=1000; int jump=1250;
    int c=0; int stnum=0; int sh1=20; int sh2=40;
    while ((c<lt) && (stnum<allst))
    {
        c++;
        if (abs(xinp(c,chn))>=bordvalue)
        {
            for (int t=0; t<samplesize; t++)
            {
                for (int i=c-int_len-sh1; i<c-sh1; i++)
                    xpre(stnum*int_len+i-c+int_len+sh1,t)=xinp(i,t);
                for (int i=c+int_len+sh2; i<c+int_len*2+sh2; i++)
                    xpost(stnum*int_len-c-int_len-sh2+i,t)=xinp(i,t);
            }
            stnum++;
            c+=jump;
        }
    }

   lcutoff = 5; hcutoff=41;
   zerophaseinit(lcutoff, hcutoff, butterord, hnt->srfr);
   getcoeffs(ac,bc,butterord);
   xpre = filtfiltmatrbytrials(bc,ac,xpre,allst,samplesize); // band-pass filtering
   xpost = filtfiltmatrbytrials(bc,ac,xpost,allst,samplesize);

  // int scale=10;

   double meanvpre;
   double meanvpost;
   double tpre;
   double tpost;

   double stdvf;
   double* valarrpre = new double[12];
   double* valarrpost = new double[12];
   for (int i=0; i<12; i++)
   {
       valarrpre[i] = 0;
       valarrpost[i] = 0;
   }

   for (int i=1; i<=blnum*50; i++)
   {
       xpre_data = getpart(xpre,i);
       xpost_data = getpart(xpost,i);
       meanvpre=0; meanvpost=0;
       for (int j=0; j<po_list[currsubj]; j++) // j<po_list[currsubj]
       {
           tpre = alphafromchannel(xpre_data,parietal_occip_list[currsubj][j],50,131);
           tpost = alphafromchannel(xpost_data,parietal_occip_list[currsubj][j],200,131);
           valarrpre[j]=tpre;
           valarrpost[j]=tpost;
           meanvpre+=tpre;
           meanvpost+=tpost;
       }
       stdvf = getstdvalue(valarrpre,po_list[currsubj]);
       meanvpre/=po_list[currsubj];
     //  meanvpre/=stdvf;
       stdvf = getstdvalue(valarrpost,po_list[currsubj]);
       meanvpost/=po_list[currsubj];
     //  meanvpost/=stdvf;
       arrpowerpre[i-1]=meanvpre;
       arrpowerpost[i-1]=meanvpost;
   }
}

void plotwindow::spectfromchannel(MatrixXf inp, int chnum, int sh, int zp, bool fl)
{
    for (int i=0; i<500; i++)
        tx[i]=0;
    for (int i=0; i<500; i++)
        ty[i]=inp(i,chnum);
    detrend(ty,tx,500);

    int start = sh;
    int len = 512;

    for (int i=start; i<start+len-zp*2; i++)
        fftarr[i-start].real(ty[i]);
    for (int i=0; i<len; i++)
        fftarr[i].imag(0);
    for (int i=0; i<zp*2; i++)  // zero-padding
        fftarr[len-zp*2+i].real(0);

    cdata = CArray(fftarr,len);
    hnt->fft(cdata);
    for (int i=0; i<len; i++)
    {
        fftframpl[i]=2*sqrt(cdata[i].real()*cdata[i].real()+cdata[i].imag()*cdata[i].imag());
        fftframpl[i]/=len;
    }

    if (fl)
    for (int p=lcutoff; p<hcutoff; p++)
        specttemp_pre[p-lcutoff]+=fftframpl[p-lcutoff];
    else
    for (int p=lcutoff; p<hcutoff; p++)
        specttemp_post[p-lcutoff]+=fftframpl[p-lcutoff];
}

void plotwindow::spectfromchannel2(MatrixXf inp, int chnum, int sh, bool fl)
{
 //   for (int i=0; i<750; i++)                   // pre: sh=150, post: sh=60
  //      tx[i]=0;
    for (int i=0; i<750; i++)
        ty[i]=inp(i,chnum);
  //  detrend(ty,tx,750);

    int lentr = 250;
    int len = 512;
    int zp = 262;

    for (int i=sh; i<sh+lentr; i++)
        fftarr[i-sh].real(ty[i-sh]);
    for (int i=0; i<len; i++)
        fftarr[i].imag(0);
    for (int i=0; i<zp; i++)  // zero-padding
        fftarr[lentr+i].real(0);

    cdata = CArray(fftarr,len);
    hnt->fft(cdata);
    for (int i=0; i<len; i++)
    {
        fftframpl[i]=2*sqrt(cdata[i].real()*cdata[i].real()+cdata[i].imag()*cdata[i].imag());
        fftframpl[i]/=len;
    }

    if (fl)
    for (int p=lcutoff; p<hcutoff; p++)
        specttemp_pre[p-lcutoff]+=fftframpl[p];
    else
    for (int p=lcutoff; p<hcutoff; p++)
        specttemp_post[p-lcutoff]+=fftframpl[p];
}

void plotwindow::meanspectfromopregion(QString fname, int blnum, int currsubj)
{
    QFile inputFile(fname);                 // determine number of channels
    inputFile.open(QIODevice::ReadOnly);
    QTextStream fin(&inputFile);
    int lc=0;
    while (!fin.atEnd())
    {
        fin.readLine();
        lc++;
    }
    inputFile.close();
    samplesize=lc;

    int bllen = 25000; int chn=23; int int_len=500;
    int lenint = bllen*blnum;
    MatrixXf xpre(lenint,samplesize);
    MatrixXf xpost(lenint,samplesize);
    MatrixXf xpre_data(int_len,samplesize);
    MatrixXf xpost_data(int_len,samplesize);

    QFile inpFile(fname);                                        // read data from file
    inpFile.open(QIODevice::ReadOnly);
    QTextStream finp(&inpFile);
    QString line; QStringList sl;
    line = finp.readLine();
    sl = line.split(QRegExp("\\s+"), QString::SkipEmptyParts);
    int length = sl.length();
    MatrixXf xinp(length,samplesize);                    // points, channels
    for (int i=0; i < length; i++)
        xinp(i,0)=sl[i].toDouble();
    for (int j=1; j<samplesize; j++)
    {
        line = finp.readLine();
        sl = line.split(QRegExp("\\s+"), QString::SkipEmptyParts);
        for (int i=0; i < length; i++)
            xinp(i,j)=sl[i].toDouble();
    }
    inpFile.close();

    int lt = sl.length();                                       // extraction of pre / post data for all channels
    int allst=blnum*50; int bordvalue=1000; int jump=1250;
    int c=0; int stnum=0; int sh1=20; int sh2=40;
    while ((c<lt) && (stnum<allst))
    {
        c++;
        if (abs(xinp(c,chn))>=bordvalue)
        {
            for (int t=0; t<samplesize; t++)
            {
                for (int i=c-int_len-sh1; i<c-sh1; i++)
                    xpre(stnum*int_len+i-c+int_len+sh1,t)=xinp(i,t);
                for (int i=c+int_len+sh2; i<c+int_len*2+sh2; i++)
                    xpost(stnum*int_len-c-int_len-sh2+i,t)=xinp(i,t);
            }
            stnum++;
            c+=jump;
        }
    }

   lcutoff = 5; hcutoff=41;
   zerophaseinit(lcutoff, hcutoff, butterord, hnt->srfr);
   getcoeffs(ac,bc,butterord);
   xpre = filtfiltmatrbytrials(bc,ac,xpre,allst,samplesize); // band-pass filtering
   xpost = filtfiltmatrbytrials(bc,ac,xpost,allst,samplesize);

   spectECpre = new double[hcutoff-lcutoff];
   spectECpost = new double[hcutoff-lcutoff];
   spectEOpre = new double[hcutoff-lcutoff];
   spectEOpost = new double[hcutoff-lcutoff];
   specttemp_pre = new double[hcutoff-lcutoff];
   specttemp_post = new double[hcutoff-lcutoff];

   for (int i=lcutoff; i<hcutoff; i++)
   {
       spectECpre[i-lcutoff]=0;
       spectECpost[i-lcutoff]=0;
       spectEOpre[i-lcutoff]=0;
       spectEOpost[i-lcutoff]=0;
   }

   for (int i=1; i<=blnum*50; i++)
   {

       for (int p=lcutoff; p<hcutoff; p++)
       {
           specttemp_pre[p-lcutoff]=0;
           specttemp_post[p-lcutoff]=0;
       }
       int nb = (i-1) / 50;

       xpre_data = getpart(xpre,i);
       xpost_data = getpart(xpost,i);

       for (int j=0; j<po_list[currsubj]; j++)
       {
           spectfromchannel(xpre_data,parietal_occip_list[currsubj][j],50,131,true);
           spectfromchannel(xpost_data,parietal_occip_list[currsubj][j],50,131,false);
       }

       for (int p=lcutoff; p<hcutoff; p++)
       {
           specttemp_pre[p-lcutoff]/=po_list[currsubj];
           specttemp_post[p-lcutoff]/=po_list[currsubj];
       }

       for (int p=lcutoff; p<hcutoff; p++)
           if (eyescondition[currsubj+1][nb]==1)
           {
               spectECpre[p-lcutoff]+=specttemp_pre[p-lcutoff];
               spectECpost[p-lcutoff]+=specttemp_post[p-lcutoff];
           }
           else if (eyescondition[currsubj+1][nb]==-1)
           {
               spectEOpre[p-lcutoff]+=specttemp_pre[p-lcutoff];
               spectEOpost[p-lcutoff]+=specttemp_post[p-lcutoff];
           }
   }

   for (int i=lcutoff; i<hcutoff; i++)
   {
       spectECpre[i-lcutoff]/=250;
       spectECpost[i-lcutoff]/=250;
       spectEOpre[i-lcutoff]/=250;
       spectEOpost[i-lcutoff]/=250;
   }
}

void plotwindow::meanspectfromopregion2(QString fname, int trnum, int blocki, int currsubj)
{
    // determine number of channels
    QFile inputFile(fname);
    inputFile.open(QIODevice::ReadOnly);
    QTextStream fin(&inputFile);
    int lc=0;
    while (!fin.atEnd())
    {
        fin.readLine();
        lc++;
    }
    inputFile.close();
    samplesize=lc;

    // read data from file in initial input matrix
    QFile inpFile(fname);
    inpFile.open(QIODevice::ReadOnly);
    QTextStream finp(&inpFile);
    QString line; QStringList sl;
    line = finp.readLine();                                     // to determine length of line
    sl = line.split(QRegExp("\\s+"), QString::SkipEmptyParts);
    int length = sl.length()-1;
    MatrixXf xinp(length,samplesize);                    // points, channels
    for (int i=0; i < length; i++)
        xinp(i,0)=sl[i+1].toDouble();
    for (int j=1; j<samplesize; j++)
    {
        line = finp.readLine();
        sl = line.split(QRegExp("\\s+"), QString::SkipEmptyParts);
        for (int i=0; i < length; i++)
            xinp(i,j)=sl[i+1].toDouble();
    }
    inpFile.close();

    // extraction of pre / post data for all channels
    int bordvalue=1500; int int_len=750; int jump=1525;
    int bllen = trnum*int_len; int chn=20;
    MatrixXf xpre(bllen,samplesize);
    MatrixXf xpost(bllen,samplesize);
    MatrixXf xpre_data(int_len,samplesize);
    MatrixXf xpost_data(int_len,samplesize);

    int c=100; int stnum=0; int sh1=20; int sh2=40;
    while (c<length-1)
    {
        c++;
        if (abs(xinp(c,chn)>=bordvalue))
        {
            for (int t=0; t<samplesize; t++)
            {
                for (int i=c-int_len-sh1; i<c-sh1; i++)
                    xpre(stnum*int_len+i-c+int_len+sh1,t)=xinp(i,t);
                for (int i=c+int_len*2+sh2; i<c+int_len*3+sh2; i++)
                    xpost(stnum*int_len+i-c-int_len*2-sh2,t)=xinp(i,t);
            }
            stnum++;
            c+=jump;
          //  qDebug()<<stnum;
        }
    }

   // filter data with broad band-pass
   lcutoff = 5; hcutoff=41;
   zerophaseinit(lcutoff, hcutoff, butterord, hnt->srfr);
   getcoeffs(ac,bc,butterord);
   xpre = filtfiltmatrbytrials(bc,ac,xpre,trnum,samplesize);
   xpost = filtfiltmatrbytrials(bc,ac,xpost,trnum,samplesize);

   spectECpre = new double[hcutoff-lcutoff];
   spectECpost = new double[hcutoff-lcutoff];
   specttemp_pre = new double[hcutoff-lcutoff];
   specttemp_post = new double[hcutoff-lcutoff];

   for (int i=lcutoff; i<hcutoff; i++)
   {
       spectECpre[i-lcutoff]=0;
       spectECpost[i-lcutoff]=0;
   }

   for (int i=1; i<=trnum; i++)
   {
       for (int p=lcutoff; p<hcutoff; p++)
       {
           specttemp_pre[p-lcutoff]=0;
           specttemp_post[p-lcutoff]=0;
       }

       xpre_data = getpart(xpre,i);
       xpost_data = getpart(xpost,i);

       for (int j=0; j<9; j++)
       {
           spectfromchannel2(xpre_data,channelsarr[currsubj][blocki][j]-1,400,true);
           spectfromchannel2(xpost_data,channelsarr[currsubj][blocki][j]-1,60,false);
       }

       for (int p=lcutoff; p<hcutoff; p++)
       {
           specttemp_pre[p-lcutoff]/=9;
           specttemp_post[p-lcutoff]/=9;
       }

       for (int p=lcutoff; p<hcutoff; p++)
       {
           spectECpre[p-lcutoff]+=specttemp_pre[p-lcutoff];
           spectECpost[p-lcutoff]+=specttemp_post[p-lcutoff];
       }

   }

   for (int i=lcutoff; i<hcutoff; i++)
   {
       spectECpre[i-lcutoff]/=trnum;
       spectECpost[i-lcutoff]/=trnum;
   }
}

void plotwindow::meanpowerfromopregion2(QString fname, int trnum, int blocki, int currsubj)
{
    // determine number of channels
    QFile inputFile(fname);
    inputFile.open(QIODevice::ReadOnly);
    QTextStream fin(&inputFile);
    int lc=0;
    while (!fin.atEnd())
    {
        fin.readLine();
        lc++;
    }
    inputFile.close();
    samplesize=lc;

    // read data from file in initial input matrix
    QFile inpFile(fname);
    inpFile.open(QIODevice::ReadOnly);
    QTextStream finp(&inpFile);
    QString line; QStringList sl;
    line = finp.readLine();                                     // to determine length of line
    sl = line.split(QRegExp("\\s+"), QString::SkipEmptyParts);
    int length = sl.length()-1;
    MatrixXf xinp(length,samplesize);                    // points, channels
    for (int i=0; i < length; i++)
        xinp(i,0)=sl[i+1].toDouble();
    for (int j=1; j<samplesize; j++)
    {
        line = finp.readLine();
        sl = line.split(QRegExp("\\s+"), QString::SkipEmptyParts);
        for (int i=0; i < length; i++)
            xinp(i,j)=sl[i+1].toDouble();
    }
    inpFile.close();

    // extraction of pre / post data for all channels
    int bordvalue=1500; int int_len=750; int jump=1525;
    int bllen = trnum*int_len; int chn=20;
    MatrixXf xpre(bllen,samplesize);
    MatrixXf xpost(bllen,samplesize);
    MatrixXf xpre_data(int_len,samplesize);
    MatrixXf xpost_data(int_len,samplesize);

    int c=100; int stnum=0; int sh1=20; int sh2=40;
    while (c<length-1)
    {
        c++;
        if (abs(xinp(c,chn)>=bordvalue))
        {
            for (int t=0; t<samplesize; t++)
            {
                for (int i=c-int_len-sh1; i<c-sh1; i++)
                    xpre(stnum*int_len+i-c+int_len+sh1,t)=xinp(i,t);
                for (int i=c+int_len*2+sh2; i<c+int_len*3+sh2; i++)
                    xpost(stnum*int_len+i-c-int_len*2-sh2,t)=xinp(i,t);
            }
            stnum++;
            c+=jump;
          //  qDebug()<<stnum;
        }
    }

   // filter data with broad band-pass
   lcutoff = 5; hcutoff=41;
   zerophaseinit(lcutoff, hcutoff, butterord, hnt->srfr);
   getcoeffs(ac,bc,butterord);
   xpre = filtfiltmatrbytrials(bc,ac,xpre,trnum,samplesize);
   xpost = filtfiltmatrbytrials(bc,ac,xpost,trnum,samplesize);

   double meanvpre;
   double meanvpost;
   double tpre;
   double tpost;

   double stdvf;
   double* valarrpre = new double[9];
   double* valarrpost = new double[9];
   for (int i=0; i<9; i++)
   {
       valarrpre[i] = 0;
       valarrpost[i] = 0;
   }

   for (int i=1; i<=trnum; i++)
   {
       xpre_data = getpart(xpre,i);
       xpost_data = getpart(xpost,i);
       meanvpre=0; meanvpost=0;
       for (int j=0; j<9; j++) // j<po_list[currsubj]
       {
           tpre = alphafromchannel(xpre_data,channelsarr[currsubj][blocki][j]-1,200,262);
           tpost = alphafromchannel(xpost_data,channelsarr[currsubj][blocki][j]-1,100,262);
           valarrpre[j]=tpre;
           valarrpost[j]=tpost;
           meanvpre+=tpre;
           meanvpost+=tpost;
       }
       stdvf = getstdvalue(valarrpre,po_list[currsubj]);
       meanvpre/=po_list[currsubj];
     //  meanvpre/=stdvf;
       stdvf = getstdvalue(valarrpost,po_list[currsubj]);
       meanvpost/=po_list[currsubj];
     //  meanvpost/=stdvf;
       arrpowerpre[i-1]=meanvpre;
       arrpowerpost[i-1]=meanvpost;
   }
}

void plotwindow::parselrtcres(QString stim)
{    
    QString fname;

    QFile outresFile1;
    QFile outresFile2;

    if (stim=="in")
    {
        outresFile1.setFileName("E:/EEGdat/analysis/LRTC/aamppre-in.dat");
        outresFile2.setFileName("E:/EEGdat/analysis/LRTC/aamppost-in.dat");
    }
    else if (stim=="anti")
    {
        outresFile1.setFileName("E:/EEGdat/analysis/LRTC/aamppre-anti.dat");
        outresFile2.setFileName("E:/EEGdat/analysis/LRTC/aamppost-anti.dat");
    }
    else if (stim=="rand")
    {
        outresFile1.setFileName("E:/EEGdat/analysis/LRTC/aamppre-rand.dat");
        outresFile2.setFileName("E:/EEGdat/analysis/LRTC/aamppost-rand.dat");
    }
    else if (stim=="noise")
    {
        outresFile1.setFileName("E:/EEGdat/analysis/LRTC/aamppre-noise.dat");
        outresFile2.setFileName("E:/EEGdat/analysis/LRTC/aamppost-noise.dat");
    }

    outresFile1.open(QIODevice::WriteOnly);
    outresFile2.open(QIODevice::WriteOnly);
    QTextStream fresout1(&outresFile1);
    QTextStream fresout2(&outresFile2);
    for (int i=0; i<19; i++)
    {
        if (i<18)
        {
            fresout1<<"s"+QString::number(i+1)<<",";
            fresout2<<"s"+QString::number(i+1)<<",";
        }
        else
        {
            fresout1<<"s"+QString::number(i+1);
            fresout2<<"s"+QString::number(i+1);
        }
    }
    fresout1<<endl;
    fresout2<<endl;

    for (int t=0; t<31; t++)
    {
        fname="E:/EEGdat/analysis/LRTC/amparr"+QString::number(t+1)+".dat";

        QFile inputFile(fname);
        inputFile.open(QIODevice::ReadOnly);
        QTextStream finp(&inputFile);
        QString line; QStringList sl;
        MatrixXf exparr(19,6);
        double arrin[19][2];
        double arranti[19][2];
        double arrrand[19][2];
        double arrnoise[19][2];
        for (int k=0; k<19; k++)
        {
            arrin[k][0]=0; arrin[k][1]=0;
            arranti[k][0]=0; arranti[k][1]=0;
            arrrand[k][0]=0; arrrand[k][1]=0;
            arrnoise[k][0]=0; arrnoise[k][1]=0;
        }
        int i1, i2;
        for (int i=0; i<19; i++)
        {
            line = finp.readLine();
            sl = line.split(QRegExp("\\,"), QString::SkipEmptyParts);
            for (int j=0; j<sl.length(); j++)
                exparr(i,j)=sl[j].toDouble();
            for (int j=0; j<4; j++)
            {
                if (j<2)
                {
                    i1=j;
                    i2=j+1;
                } else
                {
                    i1=j+1;
                    i2=j+2;
                }
                switch (blocksarr[i][j])
                {
                case 0:
                    //   qDebug()<<i1<<i2<<exparr(i,i1)<<exparr(i,i2);
                    arrnoise[i][0]=exparr(i,i1);
                    arrnoise[i][1]=exparr(i,i2);
                    break;
                case 1:
                    arrin[i][0]=exparr(i,i1);
                    arrin[i][1]=exparr(i,i2);
                    break;
                case -1:
                    arranti[i][0]=exparr(i,i1);
                    arranti[i][1]=exparr(i,i2);
                    break;
                case 3:
                    arrrand[i][0]=exparr(i,i1);
                    arrrand[i][1]=exparr(i,i2);
                    break;
                }
            }
        }
        for (int i=0; i<19; i++)
        {
            if (i<18)
            {
                if (stim=="in")
                {
                    fresout1<<QString::number(arrin[i][0],'f',3)<<",";
                    fresout2<<QString::number(arrin[i][1],'f',3)<<",";
                }
                else if (stim=="anti")
                {
                    fresout1<<QString::number(arranti[i][0],'f',3)<<",";
                    fresout2<<QString::number(arranti[i][1],'f',3)<<",";
                } else if (stim=="rand")
                {
                    fresout1<<QString::number(arrrand[i][0],'f',3)<<",";
                    fresout2<<QString::number(arrrand[i][1],'f',3)<<",";
                } else if (stim=="noise")
                {
                    fresout1<<QString::number(arrnoise[i][0],'f',3)<<",";
                    fresout2<<QString::number(arrnoise[i][1],'f',3)<<",";
                }
            } else
            {
                if (stim=="in")
                {
                    fresout1<<QString::number(arrin[i][0],'f',3);
                    fresout2<<QString::number(arrin[i][1],'f',3);
                }
                else if (stim=="anti")
                {
                    fresout1<<QString::number(arranti[i][0],'f',3);
                    fresout2<<QString::number(arranti[i][1],'f',3);
                } else if (stim=="rand")
                {
                    fresout1<<QString::number(arrrand[i][0],'f',3);
                    fresout2<<QString::number(arrrand[i][1],'f',3);
                } else if (stim=="noise")
                {
                    fresout1<<QString::number(arrnoise[i][0],'f',3);
                    fresout2<<QString::number(arrnoise[i][1],'f',3);
                }
            }
        }
        inputFile.close();
        fresout1<<endl;
        fresout2<<endl;
    }
    outresFile1.close();
    outresFile2.close();
   // cout<<exparr;
}

void plotwindow::meanplvfromopregionDT(QString fname, int blnum, int currsubj)
{
    QFile inputFile(fname);                 // determine number of channels
    inputFile.open(QIODevice::ReadOnly);
    QTextStream fin(&inputFile);
    int lc=0;
    while (!fin.atEnd())
    {
        fin.readLine();
        lc++;
    }
    inputFile.close();
    samplesize=lc;

    int bllen = 25000; int chn=23; int int_len=500;
    int lenint = bllen*blnum;
    MatrixXf xpre(lenint,samplesize);
    MatrixXf xpost(lenint,samplesize);
    MatrixXf xpre_data(int_len,samplesize);
    MatrixXf xpost_data(int_len,samplesize);

    QFile inpFile(fname);                                        // read data from file
    inpFile.open(QIODevice::ReadOnly);
    QTextStream finp(&inpFile);
    QString line; QStringList sl;
    line = finp.readLine();
    sl = line.split(QRegExp("\\s+"), QString::SkipEmptyParts);
    int length = sl.length();
    MatrixXf xinp(length,samplesize);                    // points, channels
    for (int i=0; i < length; i++)
        xinp(i,0)=sl[i].toDouble();
    for (int j=1; j<samplesize; j++)
    {
        line = finp.readLine();
        sl = line.split(QRegExp("\\s+"), QString::SkipEmptyParts);
        for (int i=0; i < length; i++)
            xinp(i,j)=sl[i].toDouble();
    }
    inpFile.close();

    int lt = sl.length();                                       // extraction of pre / post data for all channels
    int allst=blnum*50; int bordvalue=1000; int jump=1250;
    int c=0; int stnum=0; int sh1=20; int sh2=40;
    while ((c<lt) && (stnum<allst))
    {
        c++;
        if (abs(xinp(c,chn))>=bordvalue)
        {
            for (int t=0; t<samplesize; t++)
            {
                for (int i=c-int_len-sh1; i<c-sh1; i++)
                    xpre(stnum*int_len+i-c+int_len+sh1,t)=xinp(i,t);
                for (int i=c+int_len+sh2; i<c+int_len*2+sh2; i++)
                    xpost(stnum*int_len-c-int_len-sh2+i,t)=xinp(i,t);
            }
            stnum++;
            c+=jump;
        }
    }

   lcutoff = 5; hcutoff=41;
   zerophaseinit(lcutoff, hcutoff, butterord, hnt->srfr);
   getcoeffs(ac,bc,butterord);
   xpre = filtfiltmatrbytrials(bc,ac,xpre,allst,samplesize); // band-pass filtering
   xpost = filtfiltmatrbytrials(bc,ac,xpost,allst,samplesize);

  // int scale=10;

   double meanvpre;
   double meanvpost;
   double tpre;
   double tpost;

   double stdvf;
   double* valarrpre = new double[12];
   double* valarrpost = new double[12];
   for (int i=0; i<12; i++)
   {
       valarrpre[i] = 0;
       valarrpost[i] = 0;
   }

   int blocknum, trialnum;
   int trnum=50; int len = 500; int p1=500;

 //  qDebug()<<hnt->lfilt;

   for (int i=1; i<=blnum*50; i++)
   {
       blocknum = (i-1) / trnum;
       trialnum = i-1 - blocknum*trnum;
       xpre_data = getpart(xpre,i);
       xpost_data = getpart(xpost,i);
       meanvpre=0; meanvpost=0;

       extractoptphase=true;
       minphasediff((len*3)*trialnum+blockstarts[blocknum]+len,p1);
       fillstim((len*3)*trialnum+blockstarts[blocknum]+len-p1,p1); // extend stim to the left
       fillstim((len*3)*trialnum+blockstarts[blocknum]+len*2,p1); // extend stim to the right
       extractoptphase=false;

       for (int j=0; j<po_list[currsubj]; j++)
       {
           tpre = plvfromchannel(xpre_data,parietal_occip_list[currsubj][j],blocknum,trialnum,true);
           tpost = plvfromchannel(xpost_data,parietal_occip_list[currsubj][j],blocknum,trialnum,false);
           valarrpre[j]=tpre;
           valarrpost[j]=tpost;
           meanvpre+=tpre;
           meanvpost+=tpost;
       }
     //  stdvf = getstdvalue(valarrpre,po_list[currsubj]);
       meanvpre/=po_list[currsubj];
    //   meanvpre/=stdvf;
   //    stdvf = getstdvalue(valarrpost,po_list[currsubj]);
       meanvpost/=po_list[currsubj];
   //    meanvpost/=stdvf;
       prestim[i-1]=meanvpre;
       poststim[i-1]=meanvpost;
   }
   ui->widget->graph(3)->setData(arrc.xc, arrc.of1);
   ui->widget->replot();
}

double plotwindow::alphafromchannel(MatrixXf inp, int chnum, int sh, int zp)
{
    for (int i=0; i<750; i++)
        tx[i]=0;
    for (int i=0; i<750; i++)
        ty[i]=inp(i,chnum);
    detrend(ty,tx,750);
    double apow=getalphapower(ty, sh, zp);
    return apow;
}

double plotwindow::plvfromchannel(MatrixXf inp, int chnum, int blnum, int trialnum, bool pre)
{
    for (int i=0; i<500; i++)
        tx[i]=0;
    for (int i=0; i<500; i++)
        ty[i]=inp(i,chnum);
    detrend(ty,tx,500);

    int p1=300; int len=500;
    int artshift = 50;
    double plvv;

    if (pre)
    {
        double st1[p1-artshift];                                 // PLV for pre-stim
        for (int k=0; k<p1-artshift; k++)
            st1[k]=arrc.of1[(len*3)*trialnum+blockstarts[blnum]+len-p1+k];
        double st2[p1];
        for (int k=0; k<p1-artshift; k++)
            st2[k]=ty[len-p1+k];
        hnt->phsdif = phasediff(st1, st2, p1 - artshift);
        plvv = plv;
    } else
    {
        double s1[p1-artshift];                                  // PLV for post-stim
        for (int k=artshift; k<p1; k++)
            s1[k-artshift]=arrc.of1[(len*3)*trialnum+blockstarts[blnum]+len*2+k+200];
        double s2[p1];
        for (int k=artshift; k<p1; k++)
            s2[k-artshift]=ty[k+200];
        hnt->phsdif = phasediff(s1, s2, p1 - artshift);
        plvv = plv;
    }
    return plvv;
}

void plotwindow::applycsphalf(QString fname, int blnum, int comp, QString condition, QString state)
{
    QFile inputFile(fname);                 // determine number of channels
    inputFile.open(QIODevice::ReadOnly);
    QTextStream fin(&inputFile);
    int lc=0;
    while (!fin.atEnd())
    {
        fin.readLine();
        lc++;
    }
    inputFile.close();
    samplesize=lc;

    int bllen = 25000; int chn=23; int int_len=500;
    int lenint = bllen*blnum;
    MatrixXf xpre(lenint,samplesize);
    MatrixXf xpost(lenint,samplesize);
    MatrixXf xin_post(lenint/2,samplesize);
    MatrixXf xanti_post(lenint/2,samplesize);
    MatrixXf xin_pre(lenint/2,samplesize);
    MatrixXf xanti_pre(lenint/2,samplesize);

    QFile inpFile(fname);                                        // read data from file
    inpFile.open(QIODevice::ReadOnly);
    QTextStream finp(&inpFile);
    QString line; QStringList sl;
    line = finp.readLine();
    sl = line.split(QRegExp("\\s+"), QString::SkipEmptyParts);
    int length = sl.length();
    MatrixXf xinp(length,samplesize);                    // points, channels
    for (int i=0; i < length; i++)
        xinp(i,0)=sl[i].toDouble();
    for (int j=1; j<samplesize; j++)
    {
        line = finp.readLine();
        sl = line.split(QRegExp("\\s+"), QString::SkipEmptyParts);
        for (int i=0; i < length; i++)
            xinp(i,j)=sl[i].toDouble();
    }
    inpFile.close();

    int lt = sl.length();                                       // extraction of pre / post data for all channels
    int allst=blnum*50; int bordvalue=1000; int jump=1250;
    int c=0; int stnum=0; int sh1=20; int sh2=40;
    while ((c<lt) && (stnum<allst))
    {
        c++;
        if (abs(xinp(c,chn))>=bordvalue)
        {
            for (int t=0; t<samplesize; t++)
            {
                for (int i=c-int_len-sh1; i<c-sh1; i++)
                    xpre(stnum*int_len+i-c+int_len+sh1,t)=xinp(i,t);
                for (int i=c+int_len+sh2; i<c+int_len*2+sh2; i++)
                    xpost(stnum*int_len-c-int_len-sh2+i,t)=xinp(i,t);
            }
            stnum++;
            c+=jump;
        }
    }

   lcutoff = 5; hcutoff=41; // 5 41
   zerophaseinit(lcutoff, hcutoff, butterord, hnt->srfr);
   getcoeffs(ac,bc,butterord);
   MatrixXf xprebroad = filtfiltmatrbytrials(bc,ac,xpre,allst,samplesize);
   MatrixXf xpostbroad = filtfiltmatrbytrials(bc,ac,xpost,allst,samplesize);

   lcutoff = 7; hcutoff=14;
   zerophaseinit(lcutoff, hcutoff, butterord, hnt->srfr);
   getcoeffs(ac,bc,butterord);
   xpre = filtfiltmatrbytrials(bc,ac,xpre,allst,samplesize); // band-pass filtering
   xpost = filtfiltmatrbytrials(bc,ac,xpost,allst,samplesize);
   MatrixXf xprehalf(12500*blnum,samplesize); // 20000*blnum for 800ms, 12500* for 500ms
   MatrixXf xposthalf(12500*blnum,samplesize);
   int tp=250; // 400 for 800ms, 250 for 500
   for (int p=0; p<samplesize; p++)
   {
       for (int i = 0; i<allst; i++)
       for (int j = 0; j<tp; j++)
       {
           xprehalf(i*tp+j,p)=xpre(i*500+200+j,p); // +50 for 800 ms, +200 for 500 ms
           xposthalf(i*tp+j,p)=xpost(i*500+50+j,p);
       }
   }

   int in_c = -1; int anti_c = -1;                       // separation of "in" and "anti" intervals for "in/anti" condition
   if (condition=="in/anti")
   for (int k=0; k<stnum; k++)
   {
       if (relations[k]==1)
       {
           in_c++;
           for (int t=0; t<samplesize; t++)
             for (int j=0; j<int_len; j++)
             {
               xin_post(in_c*int_len+j,t)=xpost(k*int_len+j,t);
               xin_pre(in_c*int_len+j,t)=xpre(k*int_len+j,t);
             }
       }
       if (relations[k]==-1)
       {
           anti_c++;
           for (int t=0; t<samplesize; t++)
             for (int j=0; j<int_len; j++)
             {
               xanti_post(anti_c*int_len+j,t)=xpost(k*int_len+j,t);
               xanti_pre(anti_c*int_len+j,t)=xpre(k*int_len+j,t);
             }
       }
   }

   MatrixXf xpre_data(int_len,samplesize);
   MatrixXf xpost_data(int_len,samplesize);
   int_len=250; // 400 for 800ms, 250 for 500 ms
   MatrixXf xpre_filter((blnum*50-1)*int_len,samplesize);
   MatrixXf xpost_filter((blnum*50-1)*int_len,samplesize);
   MatrixXf xin_filter((blnum*25-1)*int_len,samplesize);
   MatrixXf xanti_filter((blnum*25-1)*int_len,samplesize);
   MatrixXf C_1, C_2;
   VectorXf Xcsp_C1_pre, Xcsp_C2_pre, Xcsp_C1_post, Xcsp_C2_post;
   VectorXf Xres_pre(500);
   VectorXf Xres_post(500);
   VectorXf eigvC1, eigvC2;

   if (use_mean_trial)
   for (int i=0; i<500; i++)
   {
       meanpretrial[i]=0;
       meanposttrial[i]=0;
   }

   in_c = 0; anti_c = 0;
   for (int i=1; i<=blnum*50; i++)                         // extraction and application of CSPs in similar to cross validation procedure
   {
       xpre_data = getpart(xprebroad,i);  // xprebroad for spectrums    // extract one block / trial for application
       xpost_data = getpart(xpostbroad,i);

       if (condition=="pre/post")
       {
           xpre_filter = cutpart(xprehalf,i,blnum,tp);  // extract all other blocks / trials for filter creation
           xpost_filter = cutpart(xposthalf,i,blnum,tp);
           C_1 = covarian(xpre_filter)/(blnum*50-1);
           C_2 = covarian(xpost_filter)/(blnum*50-1);
       } else
       if (condition=="in/anti")
       {
           if (relations[i-1]==1)
           {
               in_c++;
               xin_filter = cutpart(xin_post,in_c,blnum,500);
               C_1 = covarian(xin_filter)/(blnum*25-1);
               C_2 = covarian(xanti_post)/(blnum*25);
           } else
           if (relations[i-1]==-1)
           {
               anti_c++;
               xanti_filter = cutpart(xanti_post,anti_c,blnum,500);
               C_1 = covarian(xin_post)/(blnum*25);
               C_2 = covarian(xanti_filter)/(blnum*25-1);
           }
       }

       GeneralizedSelfAdjointEigenSolver<MatrixXf> ges(C_1,C_1+C_2);        // filters for C1
       //    eigvC1 = ges.eigenvalues();
       VectorXf CSPfilter1(samplesize);
       CSPfilter1=ges.eigenvectors().col(samplesize-comp);          // highest eigenvalue corresponds to last vector: col[samplesize-1]

     //  GeneralizedSelfAdjointEigenSolver<MatrixXf> gs(C_2,C_1+C_2);     // filters for C2
     //  eigvC2 = gs.eigenvalues();
       VectorXf CSPfilter2(samplesize);
       CSPfilter2=ges.eigenvectors().col(comp-1);

       double* filter1 = new double[samplesize];
       for (int t=0; t<samplesize; t++)
           filter1[t]=CSPfilter1[t];
       double stdvf = getstdvalue(filter1,samplesize);
       for (int t=0; t<samplesize; t++)
           CSPfilter1[t]/=stdvf;

       double* filter2 = new double[samplesize];
       for (int t=0; t<samplesize; t++)
           filter2[t]=CSPfilter2[t];
       stdvf = getstdvalue(filter2,samplesize);
       for (int t=0; t<samplesize; t++)
           CSPfilter2[t]/=stdvf;

       Xcsp_C1_pre = CSPfilter1.transpose() * xpre_data.transpose();
       Xcsp_C1_post = CSPfilter1.transpose() * xpost_data.transpose();
       Xcsp_C2_pre = CSPfilter2.transpose() * xpre_data.transpose();
       Xcsp_C2_post = CSPfilter2.transpose() * xpost_data.transpose();

    /*   double* filter3 = new double[1000];
       for (int t=0; t<1000; t++)
           if (t<500)
               filter3[t]=Xcsp_C1_pre[t];
           else
               filter3[t]=Xcsp_C1_post[t-500];
       stdvf = getstdvalue(filter3,1000);
       for (int t=0; t<500; t++)
       {
           Xcsp_C1_pre[t]/=stdvf;
           Xcsp_C1_post[t]/=stdvf;
       } */

       if (state=="S1")
       {
           Xres_pre = Xcsp_C1_pre;
           Xres_post = Xcsp_C1_post;
       } else
       if (state=="S2")
       {
           Xres_pre = Xcsp_C2_pre;
           Xres_post = Xcsp_C2_post;
       }

       fill_result_in_vector(Xres_pre, Xres_post, i);

       if (use_mean_trial)
           getmeantrials(Xres_pre,Xres_post);
   }

   if (use_mean_trial)
   for (int i=0; i<500; i++)
   {
       meanpretrial[i]/=allst;
       meanposttrial[i]/=allst;
   }

 /*  cout<<"Eigv C1: ";
   for (int i=0; i<4; i++)
       cout<<fixed<<setprecision(3)<<eigvC1(i)<<" ";
   cout<<"... ... ";
   for (int i=samplesize-4; i<samplesize; i++)
       cout<<fixed<<setprecision(3)<<eigvC1(i)<<" ";
   cout<<endl;

   cout<<"Eigv C2: ";
   for (int i=0; i<4; i++)
       cout<<fixed<<setprecision(3)<<eigvC2(i)<<" ";
   cout<<"... ... ";
   for (int i=samplesize-4; i<samplesize; i++)
       cout<<fixed<<setprecision(3)<<eigvC2(i)<<" ";
   cout<<endl; */

    ui->widget->graph(0)->setData(arrc.xc, arrc.amp1);
    ui->widget->graph(3)->setData(arrc.xc, arrc.of1);
  //  ui->widget->graph(5)->setData(arrc.xc, arrc.of4);
  //  ui->widget->graph(6)->setData(arrc.xc, arrc.of6);
    ui->widget->replot();
    spatial_fitering = true;
}

void plotwindow::applycsp2ndexp(QString fname, int trnum, int comp, QString state)
{
    // determine number of channels
    QFile inputFile(fname);
    inputFile.open(QIODevice::ReadOnly);
    QTextStream fin(&inputFile);
    int lc=0;
    while (!fin.atEnd())
    {
        fin.readLine();
        lc++;
    }
    inputFile.close();
    samplesize=lc;

    // read data from file in initial input matrix
    QFile inpFile(fname);
    inpFile.open(QIODevice::ReadOnly);
    QTextStream finp(&inpFile);
    QString line; QStringList sl;
    line = finp.readLine();                                     // to determine length of line
    sl = line.split(QRegExp("\\s+"), QString::SkipEmptyParts);
    int length = sl.length()-1;
    MatrixXf xinp(length,samplesize);                    // points, channels
    for (int i=0; i < length; i++)
        xinp(i,0)=sl[i+1].toDouble();
    for (int j=1; j<samplesize; j++)
    {
        line = finp.readLine();
        sl = line.split(QRegExp("\\s+"), QString::SkipEmptyParts);
        for (int i=0; i < length; i++)
            xinp(i,j)=sl[i+1].toDouble();
    }
    inpFile.close();

    // extraction of pre / post data for all channels
    int bllen = trnum*750; int chn=20;
    MatrixXf xpre(bllen,samplesize);
    MatrixXf xpost(bllen,samplesize);
    int bordvalue=1500; int int_len=750; int jump=1525;
    int c=100; int stnum=0; int sh1=20; int sh2=40;
    while (c<length-1)
    {
        c++;
        if (abs(xinp(c,chn)>=bordvalue))
        {
            for (int t=0; t<samplesize; t++)
            {
                for (int i=c-int_len-sh1; i<c-sh1; i++)
                    xpre(stnum*int_len+i-c+int_len+sh1,t)=xinp(i,t);
                for (int i=c+int_len*2+sh2; i<c+int_len*3+sh2; i++)
                    xpost(stnum*int_len+i-c-int_len*2-sh2,t)=xinp(i,t);
            }
            stnum++;
            c+=jump;
          //  qDebug()<<stnum;
        }
    }

   // filter data with broad band-pass for projecting CSP on it
   lcutoff = 5; hcutoff=41;
   zerophaseinit(lcutoff, hcutoff, butterord, hnt->srfr);
   getcoeffs(ac,bc,butterord);
   MatrixXf xprebroad = filtfiltmatrbytrials(bc,ac,xpre,trnum,samplesize);
   MatrixXf xpostbroad = filtfiltmatrbytrials(bc,ac,xpost,trnum,samplesize);

   // filter data with narrow band-pass for CSP contruction
   lcutoff = 7; hcutoff=14;
   zerophaseinit(lcutoff, hcutoff, butterord, hnt->srfr);
   getcoeffs(ac,bc,butterord);
   xpre = filtfiltmatrbytrials(bc,ac,xpre,trnum,samplesize);
   xpost = filtfiltmatrbytrials(bc,ac,xpost,trnum,samplesize);

   int xsl = 500; // length of csp-trial after cutting filtering edges

   MatrixXf xspre(xsl*trnum,samplesize);
   MatrixXf xspost(xsl*trnum,samplesize);
   for (int p=0; p<samplesize; p++)
   {
       for (int i = 0; i<trnum; i++)
       for (int j = 0; j<xsl; j++)
       {
           xspre(i*xsl+j,p)=xpre(i*int_len+200+j,p);        //  cut 100 ms before end of pre-stim
           xspost(i*xsl+j,p)=xpost(i*int_len+50+j,p);       //  cut 100 ms after end of stim
       }
   }

   MatrixXf xpre_data(int_len,samplesize);
   MatrixXf xpost_data(int_len,samplesize);
   MatrixXf xpre_filter((trnum-1)*xsl,samplesize);
   MatrixXf xpost_filter((trnum-1)*xsl,samplesize);
   MatrixXf C_1, C_2;
   VectorXf Xcsp_C1_pre, Xcsp_C2_pre, Xcsp_C1_post, Xcsp_C2_post;
   VectorXf Xres_pre(int_len);
   VectorXf Xres_post(int_len);

   for (int i=1; i<=trnum; i++)     // extraction and application of CSPs in similar to cross validation procedure
   {
       xpre_data = getpart(xprebroad,i);  // extract one block / trial for application
       xpost_data = getpart(xpostbroad,i);

       xpre_filter = cutpart(xspre,i,trnum,xsl);  // extract all other blocks / trials for filter creation
       xpost_filter = cutpart(xspost,i,trnum,xsl);

       C_1 = covarian(xpre_filter)/(trnum-1);
       C_2 = covarian(xpost_filter)/(trnum-1);

       GeneralizedSelfAdjointEigenSolver<MatrixXf> ges(C_1,C_1+C_2);        // filters for C1
       //    eigvC1 = ges.eigenvalues();
       VectorXf CSPfilter1(samplesize);
       CSPfilter1=ges.eigenvectors().col(samplesize-comp);          // highest eigenvalue corresponds to last vector: col[samplesize-1]

       VectorXf CSPfilter2(samplesize);
       CSPfilter2=ges.eigenvectors().col(comp-1);

     /*  double* filter1 = new double[samplesize];
       for (int t=0; t<samplesize; t++)
           filter1[t]=CSPfilter1[t];
       double stdvf = getstdvalue(filter1,samplesize);
       for (int t=0; t<samplesize; t++)
           CSPfilter1[t]/=stdvf;

       double* filter2 = new double[samplesize];
       for (int t=0; t<samplesize; t++)
           filter2[t]=CSPfilter2[t];
       stdvf = getstdvalue(filter2,samplesize);
       for (int t=0; t<samplesize; t++)
           CSPfilter2[t]/=stdvf; */

       Xcsp_C1_pre = CSPfilter1.transpose() * xpre_data.transpose();
       Xcsp_C1_post = CSPfilter1.transpose() * xpost_data.transpose();
       Xcsp_C2_pre = CSPfilter2.transpose() * xpre_data.transpose();
       Xcsp_C2_post = CSPfilter2.transpose() * xpost_data.transpose();


       if (state=="S1")
       {
           Xres_pre = Xcsp_C1_pre;
           Xres_post = Xcsp_C1_post;
       } else
       if (state=="S2")
       {
           Xres_pre = Xcsp_C2_pre;
           Xres_post = Xcsp_C2_post;
       }

       fill_result_in_vector(Xres_pre, Xres_post, i-1);

   }

    ui->widget->graph(1)->setData(arrc.xc, arrc.amp2);
    ui->widget->replot();
    spatial_fitering = true;
}

void plotwindow::applycsp(QString fname, int blnum, int comp, QString condition, QString state)
{
    QFile inputFile(fname);                 // determine number of channels
    inputFile.open(QIODevice::ReadOnly);
    QTextStream fin(&inputFile);
    int lc=0;
    while (!fin.atEnd())
    {
        fin.readLine();
        lc++;
    }
    inputFile.close();
    samplesize=lc;

    int bllen = 25000; int chn=23; int int_len=500;
    int lenint = bllen*blnum;
    MatrixXf xpre(lenint,samplesize);
    MatrixXf xpost(lenint,samplesize);
    MatrixXf xin_post(lenint/2,samplesize);
    MatrixXf xanti_post(lenint/2,samplesize);
    MatrixXf xin_pre(lenint/2,samplesize);
    MatrixXf xanti_pre(lenint/2,samplesize);

    QFile inpFile(fname);                                        // read data from file
    inpFile.open(QIODevice::ReadOnly);
    QTextStream finp(&inpFile);
    QString line; QStringList sl;
    line = finp.readLine();
    sl = line.split(QRegExp("\\s+"), QString::SkipEmptyParts);
    int length = sl.length();
    MatrixXf xinp(length,samplesize);                    // points, channels
    for (int i=0; i < length; i++)
        xinp(i,0)=sl[i].toDouble();    
    for (int j=1; j<samplesize; j++)
    {        
        line = finp.readLine();
        sl = line.split(QRegExp("\\s+"), QString::SkipEmptyParts);
        for (int i=0; i < length; i++)
            xinp(i,j)=sl[i].toDouble();
    } 
    inpFile.close();

    int lt = sl.length();                                       // extraction of pre / post data for all channels
    int allst=blnum*50; int bordvalue=1000; int jump=1250;
    int c=0; int stnum=0; int sh1=20; int sh2=40;
    while ((c<lt) && (stnum<allst))
    {
        c++;
        if (abs(xinp(c,chn))>=bordvalue)
        {
            for (int t=0; t<samplesize; t++)
            {
                for (int i=c-int_len-sh1; i<c-sh1; i++)
                    xpre(stnum*int_len+i-c+int_len+sh1,t)=xinp(i,t);
                for (int i=c+int_len+sh2; i<c+int_len*2+sh2; i++)
                    xpost(stnum*int_len-c-int_len-sh2+i,t)=xinp(i,t);
            }
            stnum++;
            c+=jump;
        }       
    }

   zerophaseinit(lcutoff, hcutoff, butterord, hnt->srfr);
   getcoeffs(ac,bc,butterord);
   xpre = filtfiltmatrbytrials(bc,ac,xpre,allst,samplesize); // band-pass filtering
   xpost = filtfiltmatrbytrials(bc,ac,xpost,allst,samplesize);

   int in_c = -1; int anti_c = -1;                       // separation of "in" and "anti" intervals for "in/anti" condition
   for (int k=0; k<stnum; k++)
   {
       if (relations[k]==1)
       {
           in_c++;
           for (int t=0; t<samplesize; t++)
             for (int j=0; j<int_len; j++)
             {
               xin_post(in_c*int_len+j,t)=xpost(k*int_len+j,t);
               xin_pre(in_c*int_len+j,t)=xpre(k*int_len+j,t);
             }
       }
       if (relations[k]==-1)
       {
           anti_c++;
           for (int t=0; t<samplesize; t++)
             for (int j=0; j<int_len; j++)
             {
               xanti_post(anti_c*int_len+j,t)=xpost(k*int_len+j,t);
               xanti_pre(anti_c*int_len+j,t)=xpre(k*int_len+j,t);
             }
       }
   }

   MatrixXf xpre_data(int_len,samplesize);
   MatrixXf xpost_data(int_len,samplesize);
   MatrixXf xpre_filter((blnum*50-1)*int_len,samplesize);
   MatrixXf xpost_filter((blnum*50-1)*int_len,samplesize);
   MatrixXf xin_filter((blnum*25-1)*int_len,samplesize);
   MatrixXf xanti_filter((blnum*25-1)*int_len,samplesize);
   MatrixXf C_1, C_2;
   VectorXf Xcsp_C1_pre, Xcsp_C2_pre, Xcsp_C1_post, Xcsp_C2_post;
   VectorXf Xres_pre(int_len);
   VectorXf Xres_post(int_len);
   VectorXf eigvC1, eigvC2;

   if (use_mean_trial)
   for (int i=0; i<500; i++)
   {
       meanpretrial[i]=0;
       meanposttrial[i]=0;
   }

   in_c = 0; anti_c = 0;
   for (int i=1; i<=blnum*50; i++)                         // extraction and application of CSPs in similar to cross validation procedure
   {
       xpre_data = getpart(xpre,i);       // extract one block / trial for application
       xpost_data = getpart(xpost,i);

       if (condition=="pre/post")
       {
           xpre_filter = cutpart(xpre,i,blnum,500);  // extract all other blocks / trials for filter creation
           xpost_filter = cutpart(xpost,i,blnum,500);
           C_1 = covarian(xpre_filter)/(blnum*50-1);
           C_2 = covarian(xpost_filter)/(blnum*50-1);
       } else
       if (condition=="in/anti")
       {
           if (relations[i-1]==1)
           {
               in_c++;
               xin_filter = cutpart(xin_post,in_c,blnum,500);
               C_1 = covarian(xin_filter)/(blnum*25-1);
               C_2 = covarian(xanti_post)/(blnum*25);
           } else
           if (relations[i-1]==-1)
           {
               anti_c++;
               xanti_filter = cutpart(xanti_post,anti_c,blnum,500);
               C_1 = covarian(xin_post)/(blnum*25);
               C_2 = covarian(xanti_filter)/(blnum*25-1);
           }
       }

       GeneralizedSelfAdjointEigenSolver<MatrixXf> ges(C_1,C_1+C_2);        // filters for C1
       //    eigvC1 = ges.eigenvalues();
       VectorXf CSPfilter1(samplesize);
       CSPfilter1=ges.eigenvectors().col(samplesize-comp);          // highest eigenvalue corresponds to last vector: col[samplesize-1]

     //  GeneralizedSelfAdjointEigenSolver<MatrixXf> gs(C_2,C_1+C_2);     // filters for C2
     //  eigvC2 = gs.eigenvalues();
       VectorXf CSPfilter2(samplesize);
       CSPfilter2=ges.eigenvectors().col(comp-1);

       Xcsp_C1_pre = CSPfilter1.transpose() * xpre_data.transpose();
       Xcsp_C1_post = CSPfilter1.transpose() * xpost_data.transpose();
       Xcsp_C2_pre = CSPfilter2.transpose() * xpre_data.transpose();
       Xcsp_C2_post = CSPfilter2.transpose() * xpost_data.transpose();

       if (state=="S1")
       {
           Xres_pre = Xcsp_C1_pre;
           Xres_post = Xcsp_C1_post;
       } else
       if (state=="S2")
       {
           Xres_pre = Xcsp_C2_pre;
           Xres_post = Xcsp_C2_post;
       }

       fill_result_in_vector(Xres_pre, Xres_post, i);

       if (use_mean_trial)
           getmeantrials(Xres_pre,Xres_post);
   }   

   if (use_mean_trial)
   for (int i=0; i<500; i++)
   {
       meanpretrial[i]/=allst;
       meanposttrial[i]/=allst;
   }

 /*  cout<<"Eigv C1: ";
   for (int i=0; i<4; i++)
       cout<<fixed<<setprecision(3)<<eigvC1(i)<<" ";
   cout<<"... ... ";
   for (int i=samplesize-4; i<samplesize; i++)
       cout<<fixed<<setprecision(3)<<eigvC1(i)<<" ";
   cout<<endl;

   cout<<"Eigv C2: ";
   for (int i=0; i<4; i++)
       cout<<fixed<<setprecision(3)<<eigvC2(i)<<" ";
   cout<<"... ... ";
   for (int i=samplesize-4; i<samplesize; i++)
       cout<<fixed<<setprecision(3)<<eigvC2(i)<<" ";
   cout<<endl; */

    ui->widget->graph(0)->setData(arrc.xc, arrc.amp1);
  //  ui->widget->graph(5)->setData(arrc.xc, arrc.of4);
  //  ui->widget->graph(6)->setData(arrc.xc, arrc.of6);
    ui->widget->replot();
    spatial_fitering = true;
}

void plotwindow::getmeantrials(VectorXf datapre, VectorXf datapost)
{
    for (int i=0; i<500; i++)
    {
        meanpretrial[i]+=abs(datapre[i]);
        meanposttrial[i]+=abs(datapost[i]);
    }
   // qDebug()<<"!";
}

void plotwindow::getcsptopo(QString fname, int blnum)
{
    QFile inputFile(fname);                 // determine number of channels
    inputFile.open(QIODevice::ReadOnly);
    QTextStream fin(&inputFile);
    int lc=0;
    while (!fin.atEnd())
    {
        fin.readLine();
        lc++;
    }
    inputFile.close();
    samplesize=lc;

    int bllen = 25000; int chn=23; int int_len=500;
    int lenint = bllen*blnum;
    MatrixXf xpre(lenint,samplesize);
    MatrixXf xpost(lenint,samplesize);

    QFile inpFile(fname);                                        // read data from file
    inpFile.open(QIODevice::ReadOnly);
    QTextStream finp(&inpFile);
    QString line; QStringList sl;
    line = finp.readLine();
    sl = line.split(QRegExp("\\s+"), QString::SkipEmptyParts);
    int length = sl.length();
    MatrixXf xinp(length,samplesize);                    // points, channels
    for (int i=0; i < length; i++)
        xinp(i,0)=sl[i].toDouble();
    for (int j=1; j<samplesize; j++)
    {
        line = finp.readLine();
        sl = line.split(QRegExp("\\s+"), QString::SkipEmptyParts);
        for (int i=0; i < length; i++)
            xinp(i,j)=sl[i].toDouble();
    }
    inpFile.close();

    int lt = sl.length();                                       // extraction of pre / post data for all channels
    int allst=blnum*50; int bordvalue=1000; int jump=1250;
    int c=0; int stnum=0; int sh1=20; int sh2=40;
    while ((c<lt) && (stnum<allst))
    {
        c++;
        if (abs(xinp(c,chn))>=bordvalue)
        {
            for (int t=0; t<samplesize; t++)
            {
                for (int i=c-int_len-sh1; i<c-sh1; i++)
                    xpre(stnum*int_len+i-c+int_len+sh1,t)=xinp(i,t);
                for (int i=c+int_len+sh2; i<c+int_len*2+sh2; i++)
                    xpost(stnum*int_len-c-int_len-sh2+i,t)=xinp(i,t);
            }
            stnum++;
            c+=jump;
        }
    }

   zerophaseinit(lcutoff, hcutoff, butterord, hnt->srfr);
   getcoeffs(ac,bc,butterord);
   xpre = filtfiltmatrbytrials(bc,ac,xpre,allst,samplesize);    // band-pass filtering
   xpost = filtfiltmatrbytrials(bc,ac,xpost,allst,samplesize);

   MatrixXf xprehalf(20000*blnum,samplesize);
   MatrixXf xposthalf(20000*blnum,samplesize);
   for (int p=0; p<samplesize; p++)
   {
       for (int i = 0; i<allst; i++)
       for (int j = 0; j<400; j++)
       {
           xprehalf(i*400+j,p)=xpre(i*500+50+j,p);   // for 512 ms: xprehalf(i*250+j,p)=xpre((i+1)*500-300+j,p);
           xposthalf(i*400+j,p)=xpost(i*500+50+j,p);      //             xposthalf(i*250+j,p)=xpost(i*500+50+j,p);
       }
   }

   MatrixXf C_1, C_2;
   VectorXf eigvC1, eigvC2;

   C_1 = covarian(xprehalf)/(blnum*50);
   C_2 = covarian(xposthalf)/(blnum*50);

   GeneralizedSelfAdjointEigenSolver<MatrixXf> ges(C_1,C_1+C_2);        // filters for C1
   //    eigvC1 = ges.eigenvalues();
   // VectorXf CSPfilter1(samplesize);
   // CSPfilter1=ges.eigenvectors().col(samplesize-comp);

   MatrixXf CSPfilters=ges.eigenvectors();
   MatrixXf invF = CSPfilters.inverse();

   QFile outputFile("C:/EEGdat/analysis/CSPtopo800_subj"+QString::number(currentsubj)+".dat");
   outputFile.open(QIODevice::WriteOnly);
   QTextStream fout(&outputFile);

   for (int i=0; i<invF.rows(); i++)
   {
       for (int j=0; j<invF.cols(); j++)
       {
           fout<<invF(i,j);
           if (j<invF.cols()-1)
               fout<<",";
       }
       fout<<endl;
   }
   outputFile.close();

   //  GeneralizedSelfAdjointEigenSolver<MatrixXf> gs(C_2,C_1+C_2);     // filters for C2
   //  eigvC2 = gs.eigenvalues();
   //VectorXf CSPfilter2(samplesize);
   //CSPfilter2=ges.eigenvectors().col(comp-1);

}

void plotwindow::getcsptopo2(QString fname, int trnum)
{
    QFile inputFile(fname);                 // determine number of channels
    inputFile.open(QIODevice::ReadOnly);
    QTextStream fin(&inputFile);
    int lc=0;
    while (!fin.atEnd())
    {
        fin.readLine();
        lc++;
    }
    inputFile.close();
    samplesize=lc;

    int bllen = trnum*750; int chn=20;
    MatrixXf xpre(bllen,samplesize);
    MatrixXf xpost(bllen,samplesize);

    QFile inpFile(fname);                                        // read data from file
    inpFile.open(QIODevice::ReadOnly);
    QTextStream finp(&inpFile);
    QString line; QStringList sl;
    line = finp.readLine();                                     // to determine length of line
    sl = line.split(QRegExp("\\s+"), QString::SkipEmptyParts);
    int length = sl.length()-1;
    MatrixXf xinp(length,samplesize);                    // points, channels
    for (int i=0; i < length; i++)
        xinp(i,0)=sl[i+1].toDouble();
    for (int j=1; j<samplesize; j++)
    {
        line = finp.readLine();
        sl = line.split(QRegExp("\\s+"), QString::SkipEmptyParts);
        for (int i=0; i < length; i++)
            xinp(i,j)=sl[i+1].toDouble();
    }
    inpFile.close();


    // extraction of pre / post data for all channels
    int bordvalue=1500; int int_len=750; int jump=1525;
    int c=100; int stnum=0; int sh1=20; int sh2=40;
    while (c<length-1)
    {
        c++;
        if (abs(xinp(c,chn)>=bordvalue))
        {
            for (int t=0; t<samplesize; t++)
            {
                for (int i=c-int_len-sh1; i<c-sh1; i++)
                    xpre(stnum*int_len+i-c+int_len+sh1,t)=xinp(i,t);
                for (int i=c+int_len*2+sh2; i<c+int_len*3+sh2; i++)
                    xpost(stnum*int_len+i-c-int_len*2-sh2,t)=xinp(i,t);
            }
            stnum++;
            c+=jump;
          //  qDebug()<<stnum;
        }
    }

   lcutoff=7; hcutoff=14; butterord=2;
   zerophaseinit(lcutoff, hcutoff, butterord, hnt->srfr);
   getcoeffs(ac,bc,butterord);
   xpre = filtfiltmatrbytrials(bc,ac,xpre,trnum,samplesize);    // band-pass filtering
   xpost = filtfiltmatrbytrials(bc,ac,xpost,trnum,samplesize);

   int xsl = 500; // length of csp-trial after cutting filtering edges

   MatrixXf xspre(xsl*trnum,samplesize);
   MatrixXf xspost(xsl*trnum,samplesize);
   for (int p=0; p<samplesize; p++)
   {
       for (int i = 0; i<trnum; i++)
       for (int j = 0; j<xsl; j++)
       {
           xspre(i*xsl+j,p)=xpre(i*int_len+200+j,p);        //  cut 100 ms before end of pre-stim
           xspost(i*xsl+j,p)=xpost(i*int_len+50+j,p);       //  cut 100 ms after end of stim
       }
   }

   MatrixXf C_1, C_2;
   VectorXf eigvC1, eigvC2;

   C_1 = covarian(xspre)/(trnum);
   C_2 = covarian(xspost)/(trnum);

   GeneralizedSelfAdjointEigenSolver<MatrixXf> ges(C_1,C_1+C_2);        // filters for C1
   //    eigvC1 = ges.eigenvalues();
   // VectorXf CSPfilter1(samplesize);
   // CSPfilter1=ges.eigenvectors().col(samplesize-comp);

   MatrixXf CSPfilters=ges.eigenvectors();
   MatrixXf invF = CSPfilters.inverse();

   QFile outputFile("E:/EEGdat/analysis/CSPtopo/subj"+QString::number(currentsubj)+"_block"+QString::number(currentblock)+".dat");
   outputFile.open(QIODevice::WriteOnly);
   QTextStream fout(&outputFile);

   for (int i=0; i<invF.rows(); i++)
   {
       for (int j=0; j<invF.cols(); j++)
       {
           fout<<invF(i,j);
           if (j<invF.cols()-1)
               fout<<",";
       }
       fout<<endl;
   }
   outputFile.close();

   //  GeneralizedSelfAdjointEigenSolver<MatrixXf> gs(C_2,C_1+C_2);     // filters for C2
   //  eigvC2 = gs.eigenvalues();
   //VectorXf CSPfilter2(samplesize);
   //CSPfilter2=ges.eigenvectors().col(comp-1);

}

void plotwindow::readxsddfromfile(QString fname)
{
    QFile inputFile(fname);
    inputFile.open(QIODevice::ReadOnly);
    QTextStream fin(&inputFile);
    QStringList sl;

    QString line = fin.readLine();
    line = fin.readLine();
    if (line.length()>0)
    {
        sl = line.split(QRegExp("\\,"), QString::SkipEmptyParts);
        for (int k=1; k < 1000; k++)
            arrc.amp1[hnt->posstim+k-1]=sl[k].toDouble()*20;
    }

    line = fin.readLine();
    if (line.length()>0)
    {
        sl = line.split(QRegExp("\\,"), QString::SkipEmptyParts);
        for (int k=1; k < 1000; k++)
            arrc.of1[hnt->posstim+k-1]=sl[k].toDouble()*20;
    }

    line = fin.readLine();
    if (line.length()>0)
    {
        sl = line.split(QRegExp("\\,"), QString::SkipEmptyParts);
        for (int k=1; k < 1000; k++)
            arrc.of3[hnt->posstim+k-1-10]=sl[k].toDouble()*10; // compensate phase shift
    }
    inputFile.close();
    ui->widget->graph(0)->setData(arrc.xc, arrc.amp1);
    ui->widget->graph(3)->setData(arrc.xc, arrc.of1);
    ui->widget->graph(4)->setData(arrc.xc, arrc.of3);
    ui->widget->replot();
}

void plotwindow::savedataforplv(QString fname)
{
    QFile outputFile(fname);
    outputFile.open(QIODevice::WriteOnly);
    QTextStream fout(&outputFile);
    for (int j=0; j<hnt->stlength; j++)
        fout << arrc.amp1[hnt->posstim+j] << ",";
    fout<<endl;
    for (int j=0; j<hnt->stlength; j++)
        fout << arrc.of1[hnt->posstim+j] << ",";
    fout<<endl;
    for (int j=0; j<hnt->stlength; j++)
        fout << arrc.of3[hnt->posstim+j] << ",";
    fout<<endl;
    for (int j=0; j<hnt->stlength; j++)
        fout << arrc.of4[hnt->posstim+j] << ",";
    fout<<endl;
    outputFile.close();
}

void plotwindow::extractandsaveresult(QString fname, int blockstartpos) // make parameter control for block extraction
{
    // add parameter to settings
    QFile outputFile(fname);
    outputFile.open(QIODevice::WriteOnly);
    QTextStream fout(&outputFile);
    int meanvalue=250;
    for (int i=0; i<50; i++)
        for (int j=blockstartpos+i*1500; j<blockstartpos+(i+1)*1500; j++)
            if ((((j>blockstartpos+i*1500+980) && (j<blockstartpos+i*1500+1045)) || ((j>blockstartpos+i*1500+480) && (j<blockstartpos+i*1500+510)))
                    && (abs(arrc.amp1[j])>meanvalue))
                resblock[i][j-blockstartpos-i*1500]=0;
            else
            if ((j>=blockstartpos+i*1500+500) && (j<blockstartpos+i*1500+1000))
                resblock[i][j-blockstartpos-i*1500]=arrc.of1[j];
            else
                resblock[i][j-blockstartpos-i*1500]=arrc.amp1[j];
    for (int i=0; i<50; i++)
    {
        for (int j=0; j<1500; j++)
            fout << resblock[i][j] << " ";
        fout << endl;
    }
    outputFile.close();
}

void plotwindow::loadresultblock(QString fname, int blockstartpos)
{
    QFile inputFile(fname);
    inputFile.open(QIODevice::ReadOnly);
    QTextStream fin(&inputFile);
    QStringList sl;
    QString line;

    for (int i=0; i<50; i++)
    {
        line = fin.readLine();
        sl = line.split(QRegExp("\\s+"), QString::SkipEmptyParts);
        for (int k=0; k < 1500; k++)
        {
            if ((k<500) || (k>=999))
                arrc.amp1[blockstartpos+i*1500+k]=sl[k].toDouble();
            if ((k>=500) && (k<=999))
                arrc.of1[blockstartpos+i*1500+k]=sl[k].toDouble();
        }
    }
    inputFile.close();
}


void plotwindow::sortblocks(QString fname, int blnum, bool draw)
{
    QString st = fname.mid(10,6).toLower();
    QString sst;
    for (int i=0; i<blnum; i++)
    {
        sst = fname+st+"block"+QString::number(i+1)+".dat";
        sortblock(sst,i+1,50);
        loadresultblock(sst+"sort",blockstarts[i]);
    }
    if (draw)
    {
        ui->widget->graph(0)->setData(arrc.xc, arrc.amp1);
        ui->widget->graph(3)->setData(arrc.xc, arrc.of1);
        ui->widget->replot();
    }
}

void plotwindow::sortblock(QString fname, int blocknum, int blocklength)
{
    QFile inputFile(fname);
    inputFile.open(QIODevice::ReadOnly);
    QTextStream fin(&inputFile);
    QStringList sl;
    QString line;

    QFile outputFile(fname+"sort");
    outputFile.open(QIODevice::WriteOnly);
    QTextStream fout(&outputFile);

    for (int i=0; i<blocklength; i++)
    {
        line = fin.readLine();
        sl.push_back(line);
    }
    int b = (blocknum-1)*blocklength;
    for (int i=b; i<b+blocklength; i++)
        if (relations[i]==1)
            fout<<sl.at(i-b)<<endl;
    for (int i=b; i<b+blocklength; i++)
        if (relations[i]==-1)
            fout<<sl.at(i-b)<<endl;
    inputFile.close();
    outputFile.close();
}

void plotwindow::blockinline(QString fname, int blocklength)
{
    QFile inputFile(fname);
    inputFile.open(QIODevice::ReadOnly);
    QTextStream fin(&inputFile);    
    QString line;

    QFile outputFile1(fname+"inp");
    QFile outputFile2(fname+"out");

    outputFile1.open(QIODevice::WriteOnly);
    outputFile2.open(QIODevice::WriteOnly);

    QTextStream fout1(&outputFile1);
    QTextStream fout2(&outputFile2);

    for (int i=0; i<blocklength; i++)
    {
        line = fin.readLine();
        fout1<<line;
    }

    for (int i=blocklength; i<blocklength*2; i++)
    {
        line = fin.readLine();
        fout2<<line;
    }

    inputFile.close();
    outputFile1.close();
    outputFile2.close();
}

void plotwindow::testsyntdata()
{
    stimshift=-70;
    hnt->noise=0;
    while (hnt->noise<=80)
    {
        ui->spinBox_4->setValue(hnt->noise);        
        generatesignals(hnt->posstim,hnt->imlength*hnt->numst*2,4.3,9.0,20,32.2,0,hnt->noise);
        on_pushButton_clicked();
        hnt->noise+=5;
    }

}

void plotwindow::analyzestimdata(int predperiods, int postperiods, int numbl, bool draw)
{
    // go through blocks and trials   

 //   lcutoff=5; hcutoff=41;
 //   zerophaseinit(lcutoff, hcutoff, butterord, hnt->srfr);

    int len = 500; int bl=50; hnt->lfilt=50;  // lfilt 50 for 1/2 sec intervals, 24 - for 1 period intervals
    int p1 = predperiods * hnt->srfr/hnt->osfr;
    int p2 = postperiods * hnt->srfr/hnt->osfr;
    p1 = p2 = 120; //hnt->srfr/2;
    int t1,t2,ph1,ph2,ph3;

    double meanaccsh, shiftperc;

    for (int j=0; j<numbl; j++)
    {
        meanaccsh = 0;
      //  if (eyescondition[currentsubj-1][j]==-1)
        for (int i=0; i<bl; i++)
        {
          //  zerophasefilt((len*3)*i+blockstarts[j],len);
          //  zerophasefilt((len*3)*i+len*2+blockstarts[j],len);

         /*   if (relations[j*bl+i]==1)                       // extract phase for pre-stim
                minphasediff((len*3)*i+blockstarts[j]+len-p1,p1);
            else
                maxphasediff((len*3)*i+blockstarts[j]+len-p1,p1);

            double diff = hnt->diff; // diff for optimal phase

            ph1=hnt->phase;                                 // fill prediction to the end of pre-stim interval
            extractoptphase=true;
            filloptstim((len*3)*i+blockstarts[j]+len-p1,p1,hnt->phase, 6); */

            minphasediff((len*3)*i+blockstarts[j]+len,p1);   // extend stim interval on P1 to the left
            ph2=hnt->phase;
            fillstim((len*3)*i+blockstarts[j]+len-p1,p1);

         /*   if (relations[j*bl+i]==1)                       // average relative accuracy for pre-stim
                diff=diff/hnt->phsdif;
            else
                diff=hnt->phsdif/diff;
            relativeacc[j*bl+i]=diff; */

            fillstim((len*3)*i+blockstarts[j]+len*2-1,p2+200);  // extend stim interval on P2 to the right
          //  extractoptphase=false;

            if (relations[j*bl+i]==1)                       // find optimal prediction of phase in post-stim interval
                minphasediff((len*3)*i+blockstarts[j]+len*2,p2);
            else
                maxphasediff((len*3)*i+blockstarts[j]+len*2,p2);
            ph3=hnt->phase;
            hnt->phase=ph2;                                 // fill again right extension of stim interval
            fillstim((len*3)*i+blockstarts[j]+len*2-1,p2);

            extractoptphase=true;                           // fill optimal prediction of phase in post-stim interval
            filloptstim((len*3)*i+blockstarts[j]+len*2,p2,ph3, 6);
            extractoptphase=false;

            // t1=abs(ph1-ph2);
            t2=abs(ph2-ph3);               // phase shift differences for pre and post-stim
           // if (t1>(hnt->srfr/hnt->osfr)/2)
           //     t1=hnt->srfr/hnt->osfr-t1;
            if (t2>(hnt->srfr/hnt->osfr)/2)
                t2=hnt->srfr/hnt->osfr-t2;

         //   phaseshiftdiffs[j*bl+i][0]=t1;
            phaseshiftdiffs[j*bl+i][1]=abs(t2);
            cout<<abs(t2)<<" ";
          /*  int artshift = 50;

            double st1[p1-artshift];                                 // PLV and Synch for pre-stim
            for (int k=0; k<p1-artshift; k++)
                st1[k]=arrc.of1[(len*3)*i+blockstarts[j]+len-p1+k];
            double st2[p1];
            for (int k=0; k<p1-artshift; k++)
                st2[k]=arrc.amp1[(len*3)*i+blockstarts[j]+len-p1+k];
            hnt->phsdif = phasediff(st1, st2, p1 - artshift);
            prestim[j*bl+i]=plv;
         //   synchprestim[j*bl+i]=gethindex();

            double s1[p2-artshift];                                  // PLV and Synch for post-stim
            for (int k=artshift; k<p2; k++)
                s1[k-artshift]=arrc.of1[(len*3)*i+blockstarts[j]+len*2+k+150];
            double s2[p2];
            for (int k=artshift; k<p2; k++)
                s2[k-artshift]=arrc.amp1[(len*3)*i+blockstarts[j]+len*2+k+150];
            hnt->phsdif = phasediff(s1, s2, p2 - artshift);
            poststim[j*bl+i]=plv;       */
          //  synchpoststim[j*bl+i]=gethindex();

           //if (i<bl-1)
          //      cout<<fixed<<setprecision(3)<<plvpoststim[j*bl+i]<<",";
          //  else
           //     cout<<fixed<<setprecision(3)<<plvprestim[j*bl+i];

          //  cout<<plv<<" "; */
          //  cout<<t1<<","<<t2<<" ";
          //  cout<<t2<<","<<plv<<" ";
       //   shiftperc = 1 - (double)t1/(double)(hnt->srfr/hnt->osfr)/2;
       //   relativeacc[j*bl+i] = shiftperc;
         // meanaccsh+=shiftperc;
        //  cout<<j*bl+i+1<<" ";
        //  cout<<relations[j*bl+i]<<" ";
     //     cout<<fixed<<setprecision(2)<<shiftperc<<" ";
        //  cout<<fixed<<setprecision(2)<<averagerelacc[j*bl+i]<<endl;
        }
        cout<<endl;
      //  cout<<meanaccsh/50<<endl;
      //  cout<<endl<<endl;
    }
   // qDebug()<<"1";
    //cout<<endl;

 /*   double* meanphasediff = new double [numbl];
    double* meanpoststimphasediff = new double [numbl];
    double* meanrelativeacc = new double [numbl];

    for (int i=0; i<numbl; i++)
    {
        meanphasediff[i] = 0;
        meanpoststimphasediff[i] = 0;
        meanrelativeacc[i] = 0;
    }

    for (int i=0; i<numbl; i++)
        for (int j=0; j<bl; j++)
        {
            meanphasediff[i]+=phaseshiftdiffs[i*bl+j][0];
            meanpoststimphasediff[i]+=phaseshiftdiffs[i*bl+j][1];
            meanrelativeacc[i]+=averagerelacc[i*bl+j];
        }

    for (int i=0; i<numbl; i++)
    {
        meanphasediff[i]/=bl;
        meanpoststimphasediff[i]/=bl;
        meanrelativeacc[i]/=bl;
    }

    // do for all subjs

//    double step = 360/(hnt->srfr/hnt->osfr);

    for (int i=0; i<numbl; i++)
        cout<<"Mean phase prediction diff / post-stim phase diff (block "<<i+1<<"):   "<<(int)(meanphasediff[i])<<" "<<(int)(meanpoststimphasediff[i])<<endl;
    cout<<endl;

  //  for (int i=0; i<numbl; i++)
  //      cout<<"Mean post-stim phase diff degree (block "<<i+1<<"):   "<<(int)(meanpoststimphasediff[i])<<endl;
  //  cout<<endl;

    for (int i=0; i<numbl; i++)
        cout<<"Mean relative accuracy (block "<<i+1<<"):   "<<meanrelativeacc[i]<<endl;
    cout<<endl;
*/

 /*   meanpreplvs = new double [numbl];              // Mean PLVs for blocks
    meanpostplvs = new double [numbl];
    for (int i=0; i<numbl; i++)
    {
        meanpreplvs[i] = 0;
        meanpostplvs[i] = 0;
        meanpre_in[i] = 0;
        meanpost_in[i] = 0;
        meanpre_anti[i] = 0;
        meanpost_anti[i] = 0;
    }

    for (int i=0; i<numbl; i++)
        for (int j=0; j<bl; j++)
        {
            meanpreplvs[i]+=prestim[i*bl+j];
            meanpostplvs[i]+=poststim[i*bl+j];
            if (relations[i*bl+j]==1)
            {
                meanpre_in[i]+=prestim[i*bl+j];
                meanpost_in[i]+=poststim[i*bl+j];
            } else
            if (relations[i*bl+j]==-1)
            {
                meanpre_anti[i]+=prestim[i*bl+j];
                meanpost_anti[i]+=poststim[i*bl+j];
            }
        } */

  /*  meanpresynchs = new double [numbl];             // Mean Synchs for blocks
    for (int i=0; i<numbl; i++)
        meanpresynchs[i] = 0;
    for (int i=0; i<numbl; i++)
        for (int j=0; j<bl; j++)
             meanpresynchs[i]+=synchprestim[i*bl+j];
    meanpostsynchs = new double [numbl];
    for (int i=0; i<numbl; i++)
        meanpostsynchs[i] = 0;
    for (int i=0; i<numbl; i++)
        for (int j=0; j<bl; j++)
             meanpostsynchs[i]+=synchpoststim[i*bl+j]; */

  /*  for (int i=0; i<numbl; i++)
    {
        meanpreplvs[i]/=bl;
        meanpostplvs[i]/=bl;
        meanpre_in[i]=meanpre_in[i]/(bl/2);
        meanpost_in[i]=meanpost_in[i]/(bl/2);
        meanpre_anti[i]=meanpre_anti[i]/(bl/2);
        meanpost_anti[i]=meanpost_anti[i]/(bl/2);
      //  meanpresynchs[i]/=bl;
     //   meanpostsynchs[i]/=bl;

      //  cout<<"Mean pre / post PLV for "<<i+1<<" block:   "<<meanpreplvs[i]<<"  "<<meanpostplvs[i]<<endl;
       // cout<<"Diff pre / post PLV for "<<i+1<<" block:   "<<meanpreplvs[i]-meanpostplvs[i]<<endl;
      //  cout<<"Mean pre / post Synchrony for "<<i+1<<" block:   "<<meanpresynchs[i]<<"  "<<meanpostsynchs[i]<<endl;
        //cout<<"Diff pre / post Synchrony for "<<i+1<<" block:   "<<meanpresynchs[i]-meanpostsynchs[i]<<endl;
    } */

  /*  cout<<"Mean pre PLV: ";
    for (int i=0; i<numbl; i++)
        cout<<meanpreplvs[i]<<',';
    cout<<endl;
    cout<<"Mean post PLV: ";
    for (int i=0; i<numbl; i++)
        cout<<meanpostplvs[i]<<',';
    cout<<endl;
    cout<<endl;

    for (int i=0; i<numbl; i++)
        cout<<meanpre_in[i]<<',';
    cout<<endl;
    for (int i=0; i<numbl; i++)
        cout<<meanpre_anti[i]<<',';
    cout<<endl;
    for (int i=0; i<numbl; i++)
        cout<<meanpost_in[i]<<',';
    cout<<endl;
    for (int i=0; i<numbl; i++)
        cout<<meanpost_anti[i]<<',';
    cout<<endl; */

    /*
    cout<<"Mean pre Synchs: ";
    for (int i=0; i<numbl; i++)
        cout<<meanpresynchs[i]<<',';
    cout<<endl;
    cout<<"Mean post Synchs: ";
    for (int i=0; i<numbl; i++)
        cout<<meanpostsynchs[i]<<',';
    cout<<endl; */

  if (draw)
  {
      ui->widget->graph(0)->setData(arrc.xc, arrc.amp1);
      ui->widget->graph(3)->setData(arrc.xc, arrc.of1);
      ui->widget->graph(7)->setData(arrc.xc, arrc.of5);
      ui->widget->replot();
  }
}

void plotwindow::fillmeantrials(int numbl)
{
    for (int p=0; p<1000; p++)
    {
        meantrEOin[p]=0;
        meantrEOanti[p]=0;
        meantrECin[p]=0;
        meantrECanti[p]=0;
    }

    zerophaseinit(lcutoff, hcutoff, butterord, hnt->srfr);
    int bl = 50;
    int length = 500;
    int start = 0;

    for (int j=0; j<numbl; j++)
    {
        for (int i=0; i<bl; i++)
        {
            start=blockstarts[j]+(length*3)*i;
            zerophasefilt(start,length);
            zerophasefilt(start+length*2,length);
            if ((eyescondition[currentsubj-1][j]==1) && (relations[j*bl+i]==1))
            for (int t=0; t<length; t++)
            {
                meantrECin[t]+=abs(arrc.amp1[start+t]);
                meantrECin[t+length]+=abs(arrc.amp1[start+length*2+t]);
            } else
            if ((eyescondition[currentsubj-1][j]==1) && (relations[j*bl+i]==-1))
            for (int t=0; t<length; t++)
            {
                meantrECanti[t]+=abs(arrc.amp1[start+t]);
                meantrECanti[t+length]+=abs(arrc.amp1[start+length*2+t]);
            } else
            if ((eyescondition[currentsubj-1][j]==-1) && (relations[j*bl+i]==1))
            for (int t=0; t<length; t++)
            {
                meantrEOin[t]+=abs(arrc.amp1[start+t]);
                meantrEOin[t+length]+=abs(arrc.amp1[start+length*2+t]);
            } else
            if ((eyescondition[currentsubj-1][j]==-1) && (relations[j*bl+i]==-1))
            for (int t=0; t<length; t++)
            {
                meantrEOanti[t]+=abs(arrc.amp1[start+t]);
                meantrEOanti[t+length]+=abs(arrc.amp1[start+length*2+t]);
            }
        }
    }

    int trnum = numbl*50;

    for (int t=0; t<1000; t++)
    {
        meantrEOin[t]/=trnum;
        meantrEOanti[t]/=trnum;
        meantrECin[t]/=trnum;
        meantrECanti[t]/=trnum;
    }
}

void plotwindow::analyzeallsubjsmean()
{
    int totalsubj=5;

    lcutoff = 7; hcutoff=14;

    use_mean_trial=true;

    QString fname1="C:/EEGdat/analysis/POz_EOin.dat";
    QString fname2="C:/EEGdat/analysis/POz_EOanti.dat";
    QString fname3="C:/EEGdat/analysis/POz_ECin.dat";
    QString fname4="C:/EEGdat/analysis/POz_ECanti.dat";

    QFile outFile1(fname1);
    QFile outFile2(fname2);
    QFile outFile3(fname3);
    QFile outFile4(fname4);

   /* QFile outFile_m1(fname1.left(fname1.length()-4)+"_m.dat");
    QFile outFile_m2(fname2.left(fname2.length()-4)+"_m.dat");
    double m1=0; double m2=0; */

    outFile1.open(QIODevice::WriteOnly);
    outFile2.open(QIODevice::WriteOnly);
    outFile3.open(QIODevice::WriteOnly);
    outFile4.open(QIODevice::WriteOnly);
  //  outFile_m1.open(QIODevice::WriteOnly);
  //  outFile_m2.open(QIODevice::WriteOnly);

    QTextStream fout1(&outFile1);
    QTextStream fout2(&outFile2);
    QTextStream fout3(&outFile3);
    QTextStream fout4(&outFile4);
   // QTextStream out1m(&outFile_m1);
   // QTextStream out2m(&outFile_m2);

    int nump=1000;
    for (int i=0; i<nump; i++)
    {
        fout1<<i;
        fout2<<i;
        fout3<<i;
        fout4<<i;
        if (i<nump)
        {
            fout1<<","; fout2<<",";
            fout3<<","; fout4<<",";
        }
    }
    fout1<<endl;  fout2<<endl; fout3<<endl;  fout4<<endl;

   /* for (int i=2; i<totalsubj; i++)
    {
        out1m<<"subj"<<QString::number(i);
        out2m<<"subj"<<QString::number(i);
        if (i<totalsubj-1)
        {
            out1m<<",";
            out2m<<",";
        }
    }
    out1m<<endl; out2m<<endl; */

    QString filename;
    for (int i=1; i<totalsubj; i++)
    {
        currentsubj = i;// m1=0; m2=0;
        cout<<"Subject "<<i<<endl;
        if (i<10)
            filename = "C:/EEGdat/Subj0"+QString::number(i)+"/rawdata.dat";
        else
            filename = "C:/EEGdat/Subj"+QString::number(i)+"/rawdata.dat";
        loadstimfromfile(filename);

        fsubjshort = filename.left(17);
        extractstimdata(fsubjshort,10,false);

        int numbl = 10;
        if ((i==3) || (i==5))
            numbl = 9;

        fillmeantrials(numbl);
     /*   for (int i=0; i<500; i++)
        {
            meanpretrial[i]=0;
            meanposttrial[i]=0;
        } */

       /* if ((i==3) || (i==5))
            applycsphalf(fsubjshort+"10blocks.dat",9,1,"pre/post","S2");
        else
            applycsphalf(fsubjshort+"10blocks.dat",10,1,"pre/post","S2"); */

        int st = 0; int end = 1000;

        for (int i=st; i<end; i++)
        {
            fout1<<meantrEOin[i];
            fout2<<meantrEOanti[i];
            fout3<<meantrECin[i];
            fout4<<meantrECanti[i];
            if (i<end-1)
            {
                fout1<<","; fout2<<",";
                fout3<<","; fout4<<",";
            }
        }

      /*  for (int i=st; i<end; i++)
        {
            fout1<<meanpretrial[i];
            fout2<<meanposttrial[i];
            if (i<end-1)
            {
                fout1<<",";
                fout2<<",";
            }
            m1+=meanpretrial[i];
            m2+=meanposttrial[i];
        }

        out1m<<m1/(double)(end-st);
        out2m<<m2/(double)(end-st);
        if (i<totalsubj-1)
        {
            out1m<<",";
            out2m<<",";
        } */

        fout1<<endl;  fout2<<endl; fout3<<endl;  fout4<<endl;
    }

    outFile1.close();
    outFile2.close();
    outFile3.close();
    outFile4.close();
    //outFile_m1.close();
   // outFile_m2.close();

}

void plotwindow::mergeresults(QString prefix)
{
    int totalsubjs=23;

    QString fname1, fname2, fname3, fresname;

    for (int i=2; i<totalsubjs; i++)
    {
        currentsubj = i;

        cout<<"Subject "<<i<<endl;
        if (i<10)
        {
            fname1 = "C:/EEGdat/analysis/group_res_"+prefix+"/group_res_POz/subj0"+QString::number(i)+".dat";
            fname2 = "C:/EEGdat/analysis/group_res_"+prefix+"/group_res_CSPpre/subj0"+QString::number(i)+".dat";
            fname3 = "C:/EEGdat/analysis/group_res_"+prefix+"/group_res_CSPpost/subj0"+QString::number(i)+".dat";
            fresname = "C:/EEGdat/analysis/group_res_"+prefix+"/subj0"+QString::number(i)+".dat";
        }
        else
        {
            fname1 = "C:/EEGdat/analysis/group_res_"+prefix+"/group_res_POz/subj"+QString::number(i)+".dat";
            fname2 = "C:/EEGdat/analysis/group_res_"+prefix+"/group_res_CSPpre/subj"+QString::number(i)+".dat";
            fname3 = "C:/EEGdat/analysis/group_res_"+prefix+"/group_res_CSPpost/subj"+QString::number(i)+".dat";
            fresname = "C:/EEGdat/analysis/group_res_"+prefix+"/subj"+QString::number(i)+".dat";
        }
        QFile f1(fname1);
        QFile f2(fname2);
        QFile f3(fname3);
        QFile fres(fresname);
        f1.open(QIODevice::ReadOnly);
        f2.open(QIODevice::ReadOnly);
        f3.open(QIODevice::ReadOnly);
        fres.open(QIODevice::WriteOnly);
        QTextStream ft1(&f1);
        QTextStream ft2(&f2);
        QTextStream ft3(&f3);
        QTextStream ftres(&fres);
        QStringList sl;
        QString line;

        int numbl = 10;
        if ((i==3) || (i==5))
            numbl = 9;

        for (int j=0; j<numbl*50; j++)
        {
            line = ft1.readLine();
            ftres<<line;
            line = ft2.readLine();
            sl = line.split(QRegExp("\\,"), QString::SkipEmptyParts);
            ftres<<","<<sl[4]<<","<<sl[5];
            line = ft3.readLine();
            sl = line.split(QRegExp("\\,"), QString::SkipEmptyParts);
            ftres<<","<<sl[4]<<","<<sl[5]<<endl;
        }

        f1.close();
        f2.close();
        f3.close();
        fres.close();
    }
}

void plotwindow::fillresfile(QString fresname, int numbl)
{
    QFile outresFile(fresname);
    outresFile.open(QIODevice::WriteOnly);
    QTextStream fresout(&outresFile);
    int bl;
    for (int i=0; i<50*numbl; i++)
    {
        bl = (i / 50) + 1;
        fresout<<"block "<<QString::number(bl)<<",";
        if (eyescondition[currentsubj-1][bl-1]==1)
            fresout<<"EC,";
        else
            fresout<<"EO,";
        fresout<<relations[i]<<","<<relativeacc[i]<<","<<arrpowerpre[i]<<","<<arrpowerpost[i]<<endl;    // for power
      //  fresout<<relations[i]<<","<<relativeacc[i]<<","<<prestim[i]<<","<<poststim[i]<<endl;                // for PLV
    }
    outresFile.close();
}

void plotwindow::sortresultblocks(QString prefix)
{
    QFile inresFile("C:/EEGdat/analysis/"+prefix+".dat");
    inresFile.open(QIODevice::ReadOnly);
    QTextStream fresin(&inresFile);

    QFile outresFile("C:/EEGdat/analysis/a"+prefix+".dat");
    outresFile.open(QIODevice::WriteOnly);
    QTextStream fresout(&outresFile);
    QStringList sl;
    int numbl;

    QString st = fresin.readLine();
    //fresout<<st<<endl;
    fresout<<"bl1(EC),bl2(EC),bl3(EC),bl4(EC),bl5(EC),bl1(EO),bl2(EO),bl3(EO),bl4(EO),bl5(EO)"<<endl;

    for (int i=0; i<21; i++)
    {
        st = fresin.readLine();
        sl = st.split(QRegExp("\\,"), QString::SkipEmptyParts);
        numbl = sl.length();
        for (int j=0; j<numbl; j++)
            if (eyescondition[i+1][j]==1)
                fresout<<sl[j]<<",";
        for (int j=0; j<numbl; j++)
            if (eyescondition[i+1][j]==-1)
                if (j<numbl-2)
                    fresout<<sl[j]<<",";
                else
                    fresout<<sl[j];
        fresout<<endl;
    }

    inresFile.close();
    outresFile.close();
}

void plotwindow::alignresultblocks(QString prefix)
{
    QFile inresFile("C:/EEGdat/analysis/"+prefix);
    inresFile.open(QIODevice::ReadOnly);
    QTextStream fresin(&inresFile);

    QFile outresFile("C:/EEGdat/analysis/st"+prefix);
    outresFile.open(QIODevice::WriteOnly);
    QTextStream fresout(&outresFile);
    QStringList sl;
    int numbl;

    QString st = fresin.readLine();
    fresout<<st<<endl;

    for (int i=0; i<21; i++)
    {
        st = fresin.readLine();
        if (eyescondition[i+1][0]==-1)
            fresout<<st;
        else
        {
            sl = st.split(QRegExp("\\,"), QString::SkipEmptyParts);
            numbl = sl.length();
            for (int j=1; j<numbl; j++)
                fresout<<sl[j]<<",";
            fresout<<sl[0];
        }
        fresout<<endl;
    }
    inresFile.close();
    outresFile.close();
}

void plotwindow::alignallresultblocks(QString prefix)
{
    alignresultblocks(prefix+"_pre_in.dat");
    alignresultblocks(prefix+"_post_in.dat");
    alignresultblocks(prefix+"_pre_anti.dat");
    alignresultblocks(prefix+"_post_anti.dat");
}

void plotwindow::combineresultblocks()
{
    QFile inresFile("C:/EEGdat/analysis/stPOz800_post_anti.dat");
    inresFile.open(QIODevice::ReadOnly);
    QTextStream fresin(&inresFile);

    QFile outresFile("C:/EEGdat/analysis/stPOz800.dat");
    outresFile.open(QIODevice::Append);
    QTextStream fresout(&outresFile);
    QStringList sl;
    int numbl;

    QString st = fresin.readLine();

    for (int i=0; i<22; i++)
    {
        st = fresin.readLine();
        sl = st.split(QRegExp("\\,"), QString::SkipEmptyParts);
        numbl = sl.length();
        for (int j=0; j<numbl; j++)
            fresout<<"subj"+QString::number(i+1)<<","<<QString::number(j+1)<<","<<"post"<<","<<"anti"<<","<<sl[j]<<endl;
    }
    inresFile.close();
    outresFile.close();
}

void plotwindow::getchanlabels()
{
    int totalsubjs=23;
    QString fresname;
    for (int i=2; i<totalsubjs; i++)
    {
        currentsubj = i;
        fresname = "C:/EEGdat/analysis/channelslist/subj"+QString::number(i)+"channels.txt";
        QFile inputFile(fresname);
        inputFile.open(QIODevice::ReadOnly);
        QTextStream fin(&inputFile);
        QStringList sl;
        QString line;
        int c = -1; int k = 0;
        while (!fin.atEnd())
        {
            line=fin.readLine();
            c++;
            if ((parietal_occip_list[currentsubj-2][k]==c) && (k<9))
            {
                sl=line.split(QRegExp("\\s+"), QString::SkipEmptyParts);
                cout<<sl[2].toStdString()<<" ";
                k++;
            }
        }
        cout<<endl;
    }
}

void plotwindow::geteceopower()
{
    int totalsubjs=23;
    QString fresname;

    ECin = 0; ECanti = 0; EOin = 0; EOanti = 0;
    cECin = 0; cECanti = 0; cEOin = 0; cEOanti = 0;

    for (int i=2; i<totalsubjs; i++)
    {
        currentsubj = i;

        if (i<10)
            fresname = "C:/EEGdat/analysis/data/res_500_1_POC/subj0"+QString::number(i)+".dat";
        else
            fresname = "C:/EEGdat/analysis/data/res_500_1_POC/subj"+QString::number(i)+".dat";

        int numbl = 10;
        if ((i==3) || (i==5))
            numbl = 9;

        QFile inputFile(fresname);
        inputFile.open(QIODevice::ReadOnly);
        QTextStream fin(&inputFile);
        QStringList sl;
        QString line;

        for (int k=0; k<numbl*50; k++)
        {
            line = fin.readLine();
            sl = line.split(QRegExp("\\,"), QString::SkipEmptyParts);
            if ((sl[1]=="EC") && (sl[2]=="1"))
            {
                ECin+=sl[4].toDouble()+sl[5].toDouble();
                cECin++;
            }
            else
            if ((sl[1]=="EC") && (sl[2]=="-1"))
            {
                ECanti+=sl[4].toDouble()+sl[5].toDouble();
                cECanti++;
            }
            else
            if ((sl[1]=="EO") && (sl[2]=="1"))
            {
                EOin+=sl[4].toDouble()+sl[5].toDouble();
                cEOin++;
            }
            else
            if ((sl[1]=="EO") && (sl[2]=="-1"))
            {
                EOanti+=sl[4].toDouble()+sl[5].toDouble();
                cEOanti++;
            }
        }
    }
    ECin/=cECin; ECanti/=cECanti; EOin/=cEOin; EOanti/=cEOanti;
    cout<<"EC_in: "<<ECin<<" EC_anti: "<<ECanti<<" EO_in: "<<EOin<<" EO_anti: "<<EOanti<<endl;
}

void plotwindow::analyzeblockpower(int numtr, int timew)
{
    for (int i=0; i<numtr; i++)
    {
        arrpowerpre[i]=0;
        arrpowerpost[i]=0;
    }
    lcutoff=7; hcutoff=14; butterord=2;
    zerophaseinit(lcutoff, hcutoff, butterord, hnt->srfr);

    int length = 750;
    int zp;
    int start = 0;
    int sh1;
    int sh2 = 60;  // shift from start of post-stim
    int len = 512;
    int lentr;

    if (timew==1000)
    {
        lentr = 500;
        sh1 = 150;
        zp = 12;
    }
    else if (timew==500)
    {
        lentr = 250;
        sh1 = 400;
        zp = 262;
    }
    double* frampl;
    int fr = (int)round(hnt->osfr); // IAF

    for (int j=0; j<numtr; j++)
    {
        start=j*(length*2);
        Complex t[len];
        if (!spatial_fitering)
            zerophasefilt(start,length);
        for (int i=start+sh1; i<start+sh1+lentr; i++)
            t[i-start-sh1].real(arrc.amp1[i]);

        for (int i=0; i<len; i++)
            t[i].imag(0);
        for (int i=0; i<zp; i++)  // zero-padding
            t[lentr+i].real(0);

        cdata = CArray(t,len);
        hnt->fft(cdata);
        frampl = new double[len];
        for (int i=0; i<len; i++)
        {
            frampl[i]=2*sqrt(cdata[i].real()*cdata[i].real()+cdata[i].imag()*cdata[i].imag());
            frampl[i]/=len;
        }
        arrpowerpre[j]=(frampl[fr-1]+frampl[fr]+frampl[fr+1]);

        start=length+j*(length*2);
        if (!spatial_fitering)
            zerophasefilt(start,length);
        for (int i=start+sh2; i<start+sh2+lentr; i++)
            t[i-start-sh2].real(arrc.amp1[i]);

        for (int i=0; i<len; i++)
            t[i].imag(0);
        for (int i=0; i<zp; i++)  // zero-padding
            t[lentr+i].real(0);

        cdata = CArray(t,len);
        hnt->fft(cdata);
        frampl = new double[len];
        for (int i=0; i<len; i++)
        {
            frampl[i]=2*sqrt(cdata[i].real()*cdata[i].real()+cdata[i].imag()*cdata[i].imag());
            frampl[i]/=len;
        }
        arrpowerpost[j]=(frampl[fr-1]+frampl[fr]+frampl[fr+1]);

      //  cout<<fixed<<setprecision(2)<<arrpowerpost[j]-arrpowerpre[j]<<" ";
    }

  //  for (int i=0; i<numtr; i++)
  //      cout<<i+1<<",";
  //  cout<<endl;
 /*   for (int i=0; i<numtr; i++)
        cout<<fixed<<setprecision(2)<<arrpowerpre[i]<<",";
    cout<<endl;
    for (int i=0; i<numtr; i++)
        cout<<fixed<<setprecision(2)<<arrpowerpost[i]<<",";
    cout<<endl; */
}

int plotwindow::extractblock(QString fname)
{
    QFile inpFile(fname);                                        // read data from file
    inpFile.open(QIODevice::ReadOnly);
    QTextStream finp(&inpFile);
    QString line; QStringList sl;
    int pozindex = 0;
    while (!finp.atEnd())
    {
        pozindex++;
        line = finp.readLine();
        sl = line.split(QRegExp("\\s+"), QString::SkipEmptyParts);
        line = sl[0];
     //   cout<<pozindex<<" "<<line.toStdString()<<" ";
        if (line=="POz")
            break;
    }
    cout<<endl;
    inpFile.close();
  //  qDebug()<<pozindex;
    inpFile.open(QIODevice::ReadOnly);
    for (int i=0; i<pozindex-1; i++)
        line = finp.readLine();
    line = finp.readLine();
    sl = line.split(QRegExp("\\s+"), QString::SkipEmptyParts);
    int length = sl.length()-1;
    //qDebug()<<length;
    VectorXi xinp(length);
    for (int i=0; i < length; i++)
    {
        xinp(i)=sl[i].toDouble();
   //     arrc.amp1[i]=xinp(i);
    }
    inpFile.close();

    // extraction of pre / post data for POz channels
    int bordvalue=1500; int int_len=750; int jump=1525;
    int c=100; int stnum=0; int sh1=20; int sh2=40;
    while (c<length-1)
    {
        c++;
        if (abs(xinp(c)>=bordvalue))
        {
            for (int i=c-int_len-sh1; i<c-sh1; i++)
               arrc.amp1[stnum*(2*int_len)+i-c+int_len+sh1]=xinp(i);
            for (int i=c+int_len*2+sh2; i<c+int_len*3+sh2; i++)
               arrc.amp1[int_len+stnum*(2*int_len)+i-c-int_len*2-sh2]=xinp(i);
            stnum++;
            c+=jump;
          //  qDebug()<<stnum;
        }
    }
    ui->widget->graph(0)->setData(arrc.xc, arrc.amp1);
    ui->widget->replot();
    return stnum;
}

void plotwindow::analyzeallsubjs()
{
  //  use_mean_trial=true;

    lcutoff = 5; hcutoff=41;

    int totalsubjs=23;

    QString fname1="C:/EEGdat/analysis/phasediff_EC_in.dat";
    QString fname2="C:/EEGdat/analysis/phasediff_EC_anti.dat";
    QString fname3="C:/EEGdat/analysis/phasediff_EO_in.dat";
    QString fname4="C:/EEGdat/analysis/phasediff_EO_anti.dat";

    QFile outFile1(fname1);
    QFile outFile2(fname2);
    QFile outFile3(fname3);
    QFile outFile4(fname4);

 //   QFile outFile_m1(fname1.left(fname1.length()-4)+"_m.dat");
 //   QFile outFile_m2(fname2.left(fname2.length()-4)+"_m.dat");
 //   double m1=0; double m2=0;

    outFile1.open(QIODevice::WriteOnly);
    outFile2.open(QIODevice::WriteOnly);
    outFile3.open(QIODevice::WriteOnly);
    outFile4.open(QIODevice::WriteOnly);

  //  outFile_m1.open(QIODevice::WriteOnly);
  //  outFile_m2.open(QIODevice::WriteOnly);

    QTextStream fout1(&outFile1);
    QTextStream fout2(&outFile2);
    QTextStream fout3(&outFile3);
    QTextStream fout4(&outFile4);

  //  QTextStream out1m(&outFile_m1);
  //  QTextStream out2m(&outFile_m2);

  /* for (int i=lcutoff; i<hcutoff; i++)
    {
        fout1<<i; fout2<<i; fout3<<i; fout4<<i;
        if (i<hcutoff-1)
        {
            fout1<<","; fout2<<","; fout3<<","; fout4<<",";
        }
    }
    fout1<<endl; fout2<<endl; fout3<<endl; fout4<<endl; */

  /*  for (int i=0; i<500; i++)
    {
        fout1<<i+1;
        fout2<<i+1;
        if (i<499)
        {
            fout1<<",";
            fout2<<",";
        }
    }
    fout1<<endl;
    fout2<<endl; */

    fout1<<"block1,block2,block3,block4,block5,block6,block7,block8,block9,block10"<<endl;
    fout2<<"block1,block2,block3,block4,block5,block6,block7,block8,block9,block10"<<endl;
    fout3<<"block1,block2,block3,block4,block5,block6,block7,block8,block9,block10"<<endl;
    fout4<<"block1,block2,block3,block4,block5,block6,block7,block8,block9,block10"<<endl;

  /*  for (int i=2; i<totalsubjs; i++)
    {
        out1m<<"subj"<<QString::number(i);
        out2m<<"subj"<<QString::number(i);
        if (i<totalsubjs-1)
        {
            out1m<<",";
            out2m<<",";
        }
    }

    out1m<<endl; out2m<<endl; */

    QString filename;
    QString fresname;

    for (int i=9; i<10; i++)
    {
        currentsubj = i; // m1=0; m2=0;

        cout<<"Subject "<<i<<endl;
        if (i<10)
        {
            filename = "C:/EEGdat/1stexp/Subj0"+QString::number(i)+"/rawdata.dat";
            fresname = "C:/EEGdat/analysis/subj0"+QString::number(i)+".dat";
        }
        else
        {
            filename = "C:/EEGdat/1stexp/Subj"+QString::number(i)+"/rawdata.dat";
            fresname = "C:/EEGdat/analysis/subj"+QString::number(i)+".dat";
        }

        loadstimfromfile(filename);
        fsubjshort = filename.left(24); // 17
      //  spatial_fitering = true;
      //  loadstimfromfile(fsubjshort+"/csppre.dat");
        extractstimdata(fsubjshort,10,true);

        int numbl = 10;
        if ((i==3) || (i==5))
            numbl = 9;

       // sortblocks(fsubjshort,numbl,false);
        analyzestimdata(1,1,numbl,true);

      //  applycsphalf(fsubjshort+"10blocks.dat",numbl,cspcomparr[currentsubj-2],"pre/post","S1");
      //  savedatatofile(fsubjshort+"csppost.dat");

       // getssd(1,samplenums,fsubjshort+"10blocks.dat",numbl);
      //  meandatafromopregion(fsubjshort+"10blocks.dat", numbl, currentsubj-2);
       //   meandatafromopregionDT(fsubjshort+"10blocks.dat", numbl, currentsubj-2);
      //    meanspectfromopregion(fsubjshort+"10blocks.dat", numbl, currentsubj-2);
      //  meanplvfromopregionDT(fsubjshort+"10blocks.dat", numbl, currentsubj-2);

       //  analyzepower(numbl);
      //  analyzegammapow(numbl);
       //   getpowerspectrum(numbl);
      //   fillresfile(fresname,numbl);
      //  getcsptopo(fsubjshort+"10blocks.dat",numbl);
      /*  for (int i=lcutoff; i<hcutoff; i++)
        {
            fout1<<spectECpre[i-lcutoff];
            fout2<<spectECpost[i-lcutoff];
            fout3<<spectEOpre[i-lcutoff];
            fout4<<spectEOpost[i-lcutoff];
            if (i<hcutoff-1)
            {
                fout1<<","; fout2<<",";
                fout3<<","; fout4<<",";
            }
        } */

      /*  for (int i=0; i<numbl; i++)
        {
            //  fout1<<meanpreplvs[i];
            //  fout2<<meanpostplvs[i];
            //  fout1<<meanpresynchs[i];
            //  fout2<<meanpostsynchs[i];
            // fout1<<meanprepower[i];
           //  fout2<<meanpostpower[i];
           //    fout1<<meanpre_in[i];
           //    fout2<<meanpre_anti[i];
           //    fout3<<meanpost_in[i];
           //    fout4<<meanpost_anti[i];
            fout1<<meanphsdiffEC_in[i];
            fout2<<meanphsdiffEC_anti[i];
            fout3<<meanphsdiffEO_in[i];
            fout4<<meanphsdiffEO_anti[i];
         //   m1+=meanprepower[i];
         //   m2+=meanpostpower[i];
            if (i<numbl-1)
            {
                fout1<<','; fout2<<',';
                fout3<<','; fout4<<',';
            }
        }
       /* out1m<<(m1/numbl);
        out2m<<(m2/numbl);
        if (i<totalsubjs-1)
        {
            out1m<<',';
            out2m<<',';
        } */

     /*   if (use_mean_trial)
        for (int i=0; i<500; i++)
        {
            fout1<<meanpretrial[i];
            fout2<<meanposttrial[i];
            if (i<499)
            {
                fout1<<",";
                fout2<<",";
            }
        } */

        fout1<<endl; fout2<<endl; fout3<<endl; fout4<<endl;
    }

    outFile1.close();
    outFile2.close();
    outFile3.close();
    outFile4.close();

  //  outFile_m1.close();
  //  outFile_m2.close();
  //  getstdvalue(allopttimes,opttimecount);
}

void plotwindow::extractstimdata(QString fname, int numbl, bool draw)
{
    QFile inputFile(fname+"rawdata.datlog");
    inputFile.open(QIODevice::ReadOnly);
    QTextStream fin(&inputFile);
    QString line;
    QStringList sl;
    blockstarts = new int[numbl];
    int j=numbl-1;

    while (!fin.atEnd())
    {
        line = fin.readLine();
        if (line.contains("Sequence start point:"))
        {
            sl = line.split(QRegExp("\\s+"), QString::SkipEmptyParts);
            blockstarts[j]=sl[3].toInt();
           // qDebug()<<blockstarts[j];
            j--;
        }
       /* if (line.contains("Optimization time (msec):"))
        {
            sl = line.split(QRegExp("\\s+"), QString::SkipEmptyParts);
            if (sl[0]=="Optimization")
                 allopttimes[opttimecount]=sl[3].toDouble();
            else
                 allopttimes[opttimecount]=sl[5].toDouble();
            opttimecount++;
        } */
    }
    inputFile.close();

   // for (int i=0; i<10; i++)
   //     cout<<starts[i]<<" ";
   // cout<<endl;    

 /*   QString fn = fname.left(23) + "subj" + fname.mid(14,2) +"block"; //+QString::number(1)+".dat";

    for (int i=0; i<numbl; i++)
    {
        extractandsaveresult(fn+QString::number(i+1)+".dat",blockstarts[i]);        
        loadresultblock(fn+QString::number(i+1)+".dat",blockstarts[i]);
    }
*/
    if (draw)
    {
        ui->widget->graph(0)->setData(arrc.xc, arrc.amp1);
        ui->widget->graph(3)->setData(arrc.xc, arrc.of1);
        ui->widget->replot();
    }    

    writerelationstofile(fname.left(23)+"relations.dat");
    QFile inpFile(fname+"params.txt");
    inpFile.open(QIODevice::ReadOnly);
    QTextStream fnt(&inpFile);
    line = fnt.readLine();
    sl = line.split(QRegExp("\\s+"), QString::SkipEmptyParts);
    hnt->osfr = sl[1].toDouble();
    ui->doubleSpinBox->setValue(hnt->osfr);
    inpFile.close();

}

double plotwindow::getphaseoninterval(int start, int length) // need to correct!!!
{
    zerophasefiltzc(start,length);

    Complex t[length];
    for (int i=start; i<start+length; i++)
    {
        t[i-start].real(arrc.amp2[i]);
        t[i-start].imag(0);
    }
    cdata = CArray(t,length);
    hnt->fft(cdata);
    double* frampl = new double[length];
    double* phs = new double[length];
    for (int i=0; i<length; i++)
    {
        frampl[i]=2*sqrt(cdata[i].real()*cdata[i].real()+cdata[i].imag()*cdata[i].imag());
        phs[i]=atan2(cdata[i].imag(),cdata[i].real());
    }
    int domfreq=0; double maxamp=0;
    for (int i=7; i<30; i++) // 6-29 Hz
        if (frampl[i]>maxamp)
        {
            maxamp=frampl[i];
            domfreq=i-1;
        }
    //   cout<<"Dominant freq: "<<domfreq<<" Hz with phase: "<<phs[domfreq+1]<<endl;
    return phs[domfreq+1];
}

void plotwindow::analyzegammapow(int numbl)
{
    lcutoff = 28; hcutoff = 60;
    zerophaseinit(lcutoff, hcutoff, butterord, hnt->srfr);
    int bl = 50;
    int length = 500;
    int zp = 31;
    int start = 0;
    int artsh = 25; // artifact shift
    int len = length - artsh*2 + zp*2;
    double* frampl;

    for (int i=0; i<numbl*50; i++)
    {
        arrpowerpre[i]=0;
        arrpowerpost[i]=0;
    }

    for (int j=0; j<numbl; j++)
    {
        for (int i=0; i<bl; i++)
        {
            start=blockstarts[j]+(length*3)*i+artsh; // to have same number of points as in post-stim
            Complex t[len];
            if (!spatial_fitering)
            {
                zerophasefiltzc(start,len);
                for (int i=start; i<start+len-zp*2; i++)
                {
                    t[zp+i-start].real(arrc.amp2[i]);
                    t[zp+i-start].imag(0);
                }
                for (int i=0; i<zp; i++)
                {
                    t[i].real(0); t[len-zp+i].real(0);
                    t[i].imag(0); t[len-zp+i].imag(0);
                }
            } else
            {
                for (int i=start; i<start+len-zp*2; i++)
                {
                    t[zp+i-start].real(arrc.amp1[i]);
                    t[zp+i-start].imag(0);
                }
                for (int i=0; i<zp; i++)
                {
                    t[i].real(0); t[len-zp+i].real(0);
                    t[i].imag(0); t[len-zp+i].imag(0);
                }
            }
            cdata = CArray(t,len);
            hnt->fft(cdata);
            frampl = new double[len];
            for (int i=0; i<len; i++)
            {
                frampl[i]=2*sqrt(cdata[i].real()*cdata[i].real()+cdata[i].imag()*cdata[i].imag());
                frampl[i]/=len;
            }

            for (int p=lcutoff; p<hcutoff; p++)
                arrpowerpre[j*bl+i]+=frampl[p];

            start=blockstarts[j]+(length*3)*i+length*2+artsh;

            if (!spatial_fitering)
            {
                zerophasefiltzc(start,len);
                for (int i=start; i<start+len-zp*2; i++)
                {
                    t[zp+i-start].real(arrc.amp2[i]);
                    t[zp+i-start].imag(0);
                }
                for (int i=0; i<zp; i++)
                {
                    t[i].real(0); t[len-zp+i].real(0);
                    t[i].imag(0); t[len-zp+i].imag(0);
                }
            } else
            {
                for (int i=start; i<start+len-zp*2; i++)
                {
                    t[zp+i-start].real(arrc.amp1[i]);
                    t[zp+i-start].imag(0);
                }
                for (int i=0; i<zp; i++)
                {
                    t[i].real(0); t[len-zp+i].real(0);
                    t[i].imag(0); t[len-zp+i].imag(0);
                }
            }
            cdata = CArray(t,len);
            hnt->fft(cdata);
            frampl = new double[len];
            for (int i=0; i<len; i++)
            {
                frampl[i]=2*sqrt(cdata[i].real()*cdata[i].real()+cdata[i].imag()*cdata[i].imag());
                frampl[i]/=len;
            }

            for (int p=lcutoff; p<hcutoff; p++)
                arrpowerpost[j*bl+i]+=frampl[p];

           // cout<<arrpowerpre[j*bl+i]<<" "<<arrpowerpost[j*bl+i]<<" ";
        }
      //  cout<<endl;
      //  cout<<endl<<endl;
    }

    for (int i=0; i<numbl; i++)
    {
        meanprepower[i]=0;
        meanpostpower[i]=0;
        meanpre_in[i]=0;
        meanpost_in[i]=0;
        meanpre_anti[i]=0;
        meanpost_anti[i]=0;
    }

    for (int i=0; i<numbl; i++)
        for (int j=0; j<bl; j++)
        {
            meanprepower[i]+=arrpowerpre[i*bl+j];
            meanpostpower[i]+=arrpowerpost[i*bl+j];
            if (relations[i*bl+j]==1)
            {
                meanpre_in[i]+=arrpowerpre[i*bl+j];
                meanpost_in[i]+=arrpowerpost[i*bl+j];
            } else
            if (relations[i*bl+j]==-1)
            {
                meanpre_anti[i]+=arrpowerpre[i*bl+j];
                meanpost_anti[i]+=arrpowerpost[i*bl+j];
            }
        }

    for (int i=0; i<numbl; i++)
    {
        meanprepower[i]/=bl;
        meanpostpower[i]/=bl;
        meanpre_in[i]=meanpre_in[i]/(bl/2);
        meanpost_in[i]=meanpost_in[i]/(bl/2);
        meanpre_anti[i]=meanpre_anti[i]/(bl/2);
        meanpost_anti[i]=meanpost_anti[i]/(bl/2);
    }

  /* for (int i=0; i<numbl; i++)
        cout<<meanprepower[i]<<",";
    cout<<endl;
    for (int i=0; i<numbl; i++)
        cout<<meanpostpower[i]<<",";
    cout<<endl;
    for (int i=0; i<numbl; i++)
        cout<<meanpre_in[i]<<",";
    cout<<endl;
    for (int i=0; i<numbl; i++)
        cout<<meanpre_anti[i]<<",";
    cout<<endl;
    for (int i=0; i<numbl; i++)
        cout<<meanpost_in[i]<<",";
    cout<<endl;
    for (int i=0; i<numbl; i++)
        cout<<meanpost_anti[i]<<",";
    cout<<endl; */

  //  for (int i=0; i<numbl; i++)
  //      cout<<"Mean pre / post power for "<<i+1<<" block:   "<<meanprepower[i]<<"  "<<meanpostpower[i]<<endl;
  //  cout<<endl;

  //  ui->widget->graph(1)->setData(arrc.xc, arrc.amp2);
  //  ui->widget->replot();

}

void plotwindow::getpowerspectrumrs()
{
    QString filename;
    lcutoff=5; hcutoff=41;
    zerophaseinit(lcutoff, hcutoff, butterord, hnt->srfr);

    QString fname1="C:/EEGdat/analysis/POz800_Spect_EC.dat";

    QFile outFile1(fname1);
    outFile1.open(QIODevice::WriteOnly);
    QTextStream fout1(&outFile1);

    for (int i=lcutoff; i<hcutoff; i++)
    {
        fout1<<i;
        if (i<hcutoff-1)
            fout1<<",";
    }
    fout1<<endl;

    for (int i=1; i<23; i++)
    {
        cout<<"Subject "<<i<<endl;
        currentsubj = i;
        if (i<10)
            filename = "C:/EEGdat/Subj0"+QString::number(i)+"/restingstateEC.dat";
        else
            filename = "C:/EEGdat/Subj"+QString::number(i)+"/restingstateEC.dat";
        loadrestingstate(filename);
        analyzerspowerspect();
        for (int i=lcutoff; i<hcutoff; i++)
        {
            fout1<<spectECpre[i-lcutoff];
            if (i<hcutoff-1)
                fout1<<",";
        }
        fout1<<endl;
    }
    outFile1.close();
}

void plotwindow::getpowerspectrum(int numbl)
{
    lcutoff = 5; hcutoff=41;
    zerophaseinit(lcutoff, hcutoff, butterord, hnt->srfr);
    int bl = 50;
    int length = 500;
    int zp = 56; // 56 for 800ms, 131 for 500ms
    int start = 0;
    int artsh = 50;                 // artifact shift
    int sh2 = 200; // 0                 // for CSP_half, to analyze 1/2 sec before and 1/2 sec after stim
    int len = 512; // 512, 256
    double* frampl;

    spectECpre = new double[hcutoff-lcutoff];
    spectECpost = new double[hcutoff-lcutoff];
    spectEOpre = new double[hcutoff-lcutoff];
    spectEOpost = new double[hcutoff-lcutoff];

    for (int i=lcutoff; i<hcutoff; i++)
    {
        spectECpre[i-lcutoff]=0;
        spectECpost[i-lcutoff]=0;
        spectEOpre[i-lcutoff]=0;
        spectEOpost[i-lcutoff]=0;
    }

    for (int j=0; j<numbl; j++)
    {
        for (int i=0; i<bl; i++)
        {
            start=blockstarts[j]+(length*3)*i+artsh; // to have same number of points as in post-stim
            Complex t[len];

            if (!spatial_fitering)
            {
                zerophasefiltzc(start-artsh,length);
                for (int i=start; i<start+len-zp*2; i++)
                    t[i-start].real(arrc.amp2[i]);
            } else
                for (int i=start; i<start+len-zp*2; i++)
                    t[i-start].real(arrc.amp1[i]);
            for (int i=0; i<len; i++)
                t[i].imag(0);
            for (int i=0; i<zp*2; i++)  // zero-padding
                t[len-zp*2+i].real(0);

            cdata = CArray(t,len);
            hnt->fft(cdata);
            frampl = new double[len];
            for (int i=0; i<len; i++)
            {
                frampl[i]=2*sqrt(cdata[i].real()*cdata[i].real()+cdata[i].imag()*cdata[i].imag());
                frampl[i]/=len;
            }

            for (int p=lcutoff; p<hcutoff; p++)
            if (eyescondition[currentsubj-1][j]==1)
                spectECpre[p-lcutoff]+=frampl[p];
            else if (eyescondition[currentsubj-1][j]==-1)
                spectEOpre[p-lcutoff]+=frampl[p];

            start=blockstarts[j]+(length*3)*i+length*2+artsh;

            if (!spatial_fitering)
            {
                zerophasefiltzc(start-artsh,length);
                for (int i=start; i<start+len-zp*2; i++)
                    t[i-start].real(arrc.amp2[i]);
            } else
                for (int i=start; i<start+len-zp*2; i++)
                    t[i-start].real(arrc.amp1[i]);
            for (int i=0; i<len; i++)
                t[i].imag(0);
            for (int i=0; i<zp*2; i++)  // zero-padding
                t[len-zp*2+i].real(0);

            cdata = CArray(t,len);
            hnt->fft(cdata);
            frampl = new double[len];
            for (int i=0; i<len; i++)
            {
                frampl[i]=2*sqrt(cdata[i].real()*cdata[i].real()+cdata[i].imag()*cdata[i].imag());
                frampl[i]/=len;
            }

            for (int p=lcutoff; p<hcutoff; p++)
            if (eyescondition[currentsubj-1][j]==1)
                spectECpost[p-lcutoff]+=frampl[p];
            else if (eyescondition[currentsubj-1][j]==-1)
                spectEOpost[p-lcutoff]+=frampl[p];
        }
    }

    for (int i=lcutoff; i<hcutoff; i++)
    {
        spectECpre[i-lcutoff]/=250;
        spectECpost[i-lcutoff]/=250;
        spectEOpre[i-lcutoff]/=250;
        spectEOpost[i-lcutoff]/=250;
    }

}

double plotwindow::getalphapower(double* arr, int sh, int zp)
{
    int start = sh;
    int len = 512;
    int fr = (int)round(hnt->osfr);

    for (int i=start; i<start+len-zp; i++)
        fftarr[i-start].real(arr[i]);
    for (int i=0; i<len; i++)
        fftarr[i].imag(0);
    for (int i=0; i<zp; i++)  // zero-padding
        fftarr[len-zp+i].real(0);

    cdata = CArray(fftarr,len);
    hnt->fft(cdata);
    for (int i=0; i<len; i++)
    {
        fftframpl[i]=2*sqrt(cdata[i].real()*cdata[i].real()+cdata[i].imag()*cdata[i].imag());
        fftframpl[i]/=len;
    }
    double alphapow = fftframpl[fr-1]+fftframpl[fr]+fftframpl[fr+1];
    return alphapow;
}

void plotwindow::analyzepower(int numbl)
{
    zerophaseinit(lcutoff, hcutoff, butterord, hnt->srfr);
    int bl = 50;
    int length = 500;
    int zp = 131;       // 56 - for 800 ms, 131 / 156 - for 500 ms;
    int start = 0;
    int artsh = 50;    // artifact shift
    int sh2 = 200;     // shift for pre-stim in 500 ms case
    int len = 512;
    double* frampl;
    int fr = (int)trunc(hnt->osfr); // IAF

    for (int j=0; j<numbl; j++)
    {
        for (int i=0; i<bl; i++)
        {
            start=blockstarts[j]+(length*3)*i+sh2; //sh2 // to have same number of points as in post-stim
            Complex t[len];
            if (!spatial_fitering)
            {
                zerophasefiltzc(start-sh2,length);
                for (int i=start; i<start+len-zp*2; i++)
                    t[i-start].real(arrc.amp2[i]);
            } else
                for (int i=start; i<start+len-zp*2; i++)
                    t[i-start].real(arrc.amp1[i]);

            for (int i=0; i<len; i++)
                t[i].imag(0);
            for (int i=0; i<zp*2; i++)  // zero-padding
                t[len-zp*2+i].real(0);

            cdata = CArray(t,len);
            hnt->fft(cdata);
            frampl = new double[len];
            for (int i=0; i<len; i++)
            {
                frampl[i]=2*sqrt(cdata[i].real()*cdata[i].real()+cdata[i].imag()*cdata[i].imag());
                frampl[i]/=len;
            }
            arrpowerpre[j*bl+i]=(frampl[fr-1]+frampl[fr]+frampl[fr+1]);

            start=blockstarts[j]+(length*3)*i+length*2+artsh*3;
            if (!spatial_fitering)
            {
                zerophasefiltzc(start-artsh,length);
                for (int i=start; i<start+len-zp*2; i++)
                    t[i-start].real(arrc.amp2[i]);
            } else
                for (int i=start; i<start+len-zp*2; i++)
                    t[i-start].real(arrc.amp1[i]);

            for (int i=0; i<len; i++)
                t[i].imag(0);
            for (int i=0; i<zp*2; i++)  // zero-padding
                t[len-zp*2+i].real(0);

            cdata = CArray(t,len);
            hnt->fft(cdata);
            frampl = new double[len];
            for (int i=0; i<len; i++)
            {
                frampl[i]=2*sqrt(cdata[i].real()*cdata[i].real()+cdata[i].imag()*cdata[i].imag());
                frampl[i]/=len;
            }
            arrpowerpost[j*bl+i]=(frampl[fr-1]+frampl[fr]+frampl[fr+1]);
        }
        //  cout<<endl;
        //  cout<<endl<<endl;
    }

    for (int i=0; i<numbl; i++)
    {
        meanprepower[i]=0;
        meanpostpower[i]=0;
        meanpre_in[i]=0;
        meanpost_in[i]=0;
        meanpre_anti[i]=0;
        meanpost_anti[i]=0;
    }

    for (int i=0; i<numbl; i++)
        for (int j=0; j<bl; j++)
        {
            meanprepower[i]+=arrpowerpre[i*bl+j];
            meanpostpower[i]+=arrpowerpost[i*bl+j];
            if (relations[i*bl+j]==1)
            {
                meanpre_in[i]+=arrpowerpre[i*bl+j];
                meanpost_in[i]+=arrpowerpost[i*bl+j];
            } else
            if (relations[i*bl+j]==-1)
            {
                meanpre_anti[i]+=arrpowerpre[i*bl+j];
                meanpost_anti[i]+=arrpowerpost[i*bl+j];
            }
        }

    for (int i=0; i<numbl; i++)
    {
        meanprepower[i]/=bl;
        meanpostpower[i]/=bl;
        meanpre_in[i]=meanpre_in[i]/(bl/2);
        meanpost_in[i]=meanpost_in[i]/(bl/2);
        meanpre_anti[i]=meanpre_anti[i]/(bl/2);
        meanpost_anti[i]=meanpost_anti[i]/(bl/2);
    }

    /* for (int i=0; i<numbl; i++)
          cout<<meanprepower[i]<<",";
      cout<<endl;
      for (int i=0; i<numbl; i++)
          cout<<meanpostpower[i]<<",";
      cout<<endl;
      for (int i=0; i<numbl; i++)
          cout<<meanpre_in[i]<<",";
      cout<<endl;
      for (int i=0; i<numbl; i++)
          cout<<meanpre_anti[i]<<",";
      cout<<endl;
      for (int i=0; i<numbl; i++)
          cout<<meanpost_in[i]<<",";
      cout<<endl;
      for (int i=0; i<numbl; i++)
          cout<<meanpost_anti[i]<<",";
      cout<<endl; */

    //  for (int i=0; i<numbl; i++)
    //      cout<<"Mean pre / post power for "<<i+1<<" block:   "<<meanprepower[i]<<"  "<<meanpostpower[i]<<endl;
    //  cout<<endl;
}


void plotwindow::analyzerspowerspect()
{
  //  lcutoff=5; hcutoff=41;
  //  zerophaseinit(lcutoff, hcutoff, butterord, hnt->srfr);
    int bl = 50;  
    int length = 500;
    int zp = 56; //  131 for 500ms
    int start = 0;
    int artsh = 50;                 // artifact shift
    int sh2 = 200;                  // for CSP_half, to analyze 1/2 sec before and 1/2 sec after stim
    int len = 512;

    double* frampl;

    spectECpre = new double[hcutoff-lcutoff];

    for (int i=lcutoff; i<hcutoff; i++)
        spectECpre[i-lcutoff]=0;

    for (int i=0; i<bl; i++)
    {
        start=(length*2)*i+artsh;
        zerophasefiltzc(start-artsh,length);
        Complex t[len];
        for (int i=start; i<start+len-zp*2; i++)
            t[i-start].real(arrc.amp2[i]);

        for (int i=0; i<len; i++)
            t[i].imag(0);
        for (int i=0; i<zp*2; i++)  // zero-padding
            t[len-zp*2+i].real(0);

        cdata = CArray(t,len);
        hnt->fft(cdata);
        frampl = new double[len];
        for (int i=0; i<len; i++)
        {
            frampl[i]=2*sqrt(cdata[i].real()*cdata[i].real()+cdata[i].imag()*cdata[i].imag());
            frampl[i]/=len;
        }
        for (int p=lcutoff; p<hcutoff; p++)
            spectECpre[p-lcutoff]+=frampl[p];

        start=(length*2)*i+length+artsh;
        zerophasefiltzc(start-artsh,length);
        for (int i=start; i<start+len-zp*2; i++)
            t[i-start].real(arrc.amp2[i]);

        for (int i=0; i<len; i++)
            t[i].imag(0);
        for (int i=0; i<zp*2; i++)  // zero-padding
            t[len-zp*2+i].real(0);

        cdata = CArray(t,len);
        hnt->fft(cdata);
        frampl = new double[len];
        for (int i=0; i<len; i++)
        {
            frampl[i]=2*sqrt(cdata[i].real()*cdata[i].real()+cdata[i].imag()*cdata[i].imag());
            frampl[i]/=len;
        }

        for (int p=lcutoff; p<hcutoff; p++)
            spectECpre[p-lcutoff]+=frampl[p];
    }

    for (int i=lcutoff; i<hcutoff; i++)
        spectECpre[i-lcutoff]/=100;
}

void plotwindow::analyzerspower()
{
    lcutoff=5; hcutoff=41;
    zerophaseinit(lcutoff, hcutoff, butterord, hnt->srfr);
    int bl = 50;
    int length = 500;
    int zp = 56; // 56 for 800ms
    int start = 0;    
    int artsh = 50; // artifact shift
    int sh2 = 200; // shift for pre-stim in 500ms case
    int len = 512;
    double* frampl;
    int fr = (int)trunc(hnt->osfr); // IAF

    for (int i=0; i<bl; i++)
    {
        start=(length*2)*i+artsh; // sh2 - in case of 500 ms, artsh - 800 ms
        zerophasefiltzc(start-artsh,length);
        Complex t[len];
        for (int i=start; i<start+len-zp*2; i++)
            t[i-start].real(arrc.amp2[i]);

        for (int i=0; i<len; i++)
            t[i].imag(0);
        for (int i=0; i<zp*2; i++)  // zero-padding
            t[len-zp*2+i].real(0);

        cdata = CArray(t,len);
        hnt->fft(cdata);
        frampl = new double[len];
        for (int i=0; i<len; i++)
        {
            frampl[i]=2*sqrt(cdata[i].real()*cdata[i].real()+cdata[i].imag()*cdata[i].imag());
            frampl[i]/=len;
        }
        arrpowerpre[i]=(frampl[fr-1]+frampl[fr]+frampl[fr+1]);

        start=(length*2)*i+length+artsh;
        zerophasefiltzc(start-artsh,length);
        for (int i=start; i<start+len-zp*2; i++)
            t[i-start].real(arrc.amp2[i]);

        for (int i=0; i<len; i++)
            t[i].imag(0);
        for (int i=0; i<zp*2; i++)  // zero-padding
            t[len-zp*2+i].real(0);

        cdata = CArray(t,len);
        hnt->fft(cdata);
        frampl = new double[len];
        for (int i=0; i<len; i++)
        {
            frampl[i]=2*sqrt(cdata[i].real()*cdata[i].real()+cdata[i].imag()*cdata[i].imag());
            frampl[i]/=len;
        }
        arrpowerpost[i]=(frampl[fr-1]+frampl[fr]+frampl[fr+1]);
    }
    //  cout<<endl;
    //  cout<<endl<<endl;

    meanprepower[0]=0;
    meanpostpower[0]=0;

    for (int i=0; i<bl; i++)
    {
        meanprepower[0]+=arrpowerpre[i];
        meanpostpower[0]+=arrpowerpost[i];
    }

    meanprepower[0]/=bl;
    meanpostpower[0]/=bl;

    double rspower = (meanprepower[0]+meanpostpower[0])/2.0;

    if ((currentsubj+1==3) || (currentsubj+1==5))
        for (int i=0; i<9; i++)
            cout<<rspower<<",";
    else
        for (int i=0; i<10; i++)
            cout<<rspower<<",";
    cout<<endl;

   // cout<<"Mean pre / post power for RS:   "<<meanprepower[0]<<"  "<<meanpostpower[0]<<endl;
  //  cout<<endl;

}

void plotwindow::phasepreservation(int blocknum, int nump)
{
    zerophaseinit(lcutoff, hcutoff, butterord, hnt->srfr);
    int p = hnt->srfr/hnt->osfr;
    int numbl=50; int arsh=20; int t=p;//2*p
    int presh=250; // shift closer to stim start in pre-stim to eliminate possible after-stim effect from previous trial
    int blstart=blockstarts[blocknum];
    for (int k=0; k<numbl; k++)
    {
        preinstphasediffs[k][0]=getphaseoninterval(blstart+presh+k*1500,2*t);
        postinstphasediffs[k][0]=getphaseoninterval(blstart+(k+1)*1500-500+arsh,2*t);
        for (int i=0; i<nump; i++)
        {
            preinstphasediffs[k][i+1]=preinstphasediffs[k][0]-getphaseoninterval(blstart+presh+k*1500+t+i*p,2*p);
            postinstphasediffs[k][i+1]=postinstphasediffs[k][0]-getphaseoninterval(blstart+(k+1)*1500-500+arsh+t+i*p,2*p);
        }
    }
   /* for (int k=0; k<numbl; k++)
    {
        cout<<"Trial "<<k+1<<":  ";
        for (int j=1; j<nump+1; j++)
            cout<<fixed<<setprecision(3)<<instphasediffs[k][j]<<"  ";
        cout<<endl;
    } */
    for (int i=0; i<nump; i++)
    {
        complex<double> postsumexp(0.0,0.0);
        complex<double> presumexp(0.0,0.0);
        for (int k=0; k<numbl; k++)
        {
            presumexp+=exp(I*preinstphasediffs[k][i+1]);
            postsumexp+=exp(I*postinstphasediffs[k][i+1]);
        }
        preppi[i]=fabs(presumexp)/numbl;
        postppi[i]=fabs(postsumexp)/numbl;
    }
   // cout<<endl;
  //  cout<<"PPI diffs for "<<blocknum+1<<" block: ";
    // for (int i=0; i<nump; i++)
     //   if (i==nump-1)
     //       cout<<fixed<<setprecision(3)<<postppi[i]-preppi[i];
      //  else
         //   cout<<fixed<<setprecision(3)<<preppi[1]<<',';
   // cout<<endl;
}

void plotwindow::ppiforallsubjs()
{
    QFile outFile1("C:/EEGdat/analysis/_preppi_1P.dat");
    QFile outFile2("C:/EEGdat/analysis/_postppi_1P.dat");
    QFile outFile3("C:/EEGdat/analysis/_preppi_2P.dat");
    QFile outFile4("C:/EEGdat/analysis/_postppi_2P.dat");
    outFile1.open(QIODevice::WriteOnly);
    outFile2.open(QIODevice::WriteOnly);
    outFile3.open(QIODevice::WriteOnly);
    outFile4.open(QIODevice::WriteOnly);
    QTextStream fout1(&outFile1);
    QTextStream fout2(&outFile2);
    QTextStream fout3(&outFile3);
    QTextStream fout4(&outFile4);
    QString filename;
    for (int i=1; i<21; i++)
    {
        cout<<"Subject "<<i<<endl;
        if (i<10)
            filename = "C:/EEGdat/Subj0"+QString::number(i)+"/rawdata.dat";
        else
            filename = "C:/EEGdat/Subj"+QString::number(i)+"/rawdata.dat";
        hnt->firstinit(filename,NMAX);
        loadstimfromfile(filename);
        for (int i=0; i<hnt->npt; i++)
        {
            arrc.xc[i] = i;
            arrc.amp1[i]=hnt->x[i];
        }
        fsubjshort = filename.left(17);
        extractstimdata(fsubjshort,10,false);
        if ((i==3) || (i==5))
        {
            for (int i=0; i<9; i++)
            {
                phasepreservation(i,2);
                fout1<<preppi[0]<<',';
                fout2<<postppi[0]<<',';
                fout3<<preppi[1]<<',';
                fout4<<postppi[1]<<',';
            }
        } else
        {
            for (int i=0; i<10; i++)
            {
                phasepreservation(i,2);
                fout1<<preppi[0]<<',';
                fout2<<postppi[0]<<',';
                fout3<<preppi[1]<<',';
                fout4<<postppi[1]<<',';
            }
        }
        fout1<<endl;  fout2<<endl;
        fout3<<endl;  fout4<<endl;
       // cout<<endl;
    }
}

void plotwindow::loadrestingstate(QString fname)
{
    QFile inpFile(fname.left(17)+"params.txt");
    inpFile.open(QIODevice::ReadOnly);
    QTextStream fnt(&inpFile);
    QString line;
    QStringList sl;
    line = fnt.readLine();
    sl = line.split(QRegExp("\\s+"), QString::SkipEmptyParts);
    hnt->osfr = sl[1].toDouble();
    ui->doubleSpinBox->setValue(hnt->osfr);
    inpFile.close();

    QFile inputFile(fname);
    inputFile.open(QIODevice::ReadOnly);
    QTextStream fin(&inputFile);

    line = fin.readLine();
    int leng=0;

    if (line.length()>0)
    {
        sl = line.split(QRegExp("\\s+"), QString::SkipEmptyParts);
        leng=sl.size();
        for (int k=0; k < leng; k++)
            arrc.amp1[k]=sl[k].toDouble();
    }

 //   for (int i=0; i<100; i++)
 //       for (int j=0; j<500; j++)
 //           rsinitarray[i][j]=arrc.amp1[i*500+j];

    inputFile.close();

    cleareeg(leng+1, NMAX);
    ui->widget->graph(0)->setData(arrc.xc, arrc.amp1);
    ui->widget->replot();
}

void plotwindow::getrandompermutation()
{
    int numbl=100; int len=500;
    random_shuffle(&indexesforperms[0],&indexesforperms[numbl]);
    for (int i=0; i<numbl; i++)
        for (int j=0; j<len; j++)
        arrc.amp1[i*len+j]=rsinitarray[indexesforperms[i]][j];
  //  ui->widget->graph(0)->setData(arrc.xc, arrc.amp1);
  //  ui->widget->replot();
}

double plotwindow::getpvalue(int size, double value, bool side)
{
    int c=0;
    if (side)
    {
        for (int i=0; i<size; i++)
            if (perms[i]>value)
                c++;
    }
    else
    {
        for (int i=0; i<size; i++)
            if (perms[i]<value)
                c++;
    }
    return (double)c/size;
}

void plotwindow::analyzerestingstate()
{
    int len = 250; int trials=50;
    int p1, p2;
    p1 = p2 = hnt->srfr/2;

    int eyescond = -1;
    int low = 0; int high = hnt->srfr/hnt->osfr;
    hnt->lfilt = 50;

    double meanperm = 0;
    int size=100; int numbl=10; // size - number of permutations
  //  for (int ip=0; ip<size; ip++)
  //  {
      //  getrandompermutation();

        for (int i=0; i<trials; i++)
        {
            double pres = 0;
            double posts = 0;
            int k = 50;
            for (int j=0; j<k; j++)
            {
                hnt->phase=qrand() % ((high + 1) - low) + low;

             //   minphasediff((len*4)*i,len);
                fillstim((len*4)*i,p1);
               // prestim[i]=gethindex(); // plv;
                pres+=plv;// gethindex(); //plv;        

                fillstim((len*4)*i+len*3,p2); // +len*3, bl=50
               // poststim[i]=gethindex(); //plv;
                posts+=plv;//gethindex(); //plv;
            }
            prestim[i]=pres/k;
            poststim[i]=posts/k;
            // if (i<bl-1)
             //  cout<<fixed<<setprecision(3)<<prestim[i]<<","<<poststim[i]<<",";
            //  else
            //  cout<<fixed<<setprecision(3)<<prestim[i]<<","<<poststim[i];
        }
        double meanpre = 0;              // Mean PLV
        for (int j=0; j<trials; j++)
            meanpre+=prestim[j];
        double meanpost = 0;
        meanpost = 0;
        for (int j=0; j<trials; j++)
            meanpost+=poststim[j];
        meanpre/=trials;
        meanpost/=trials;

        meanperm+=meanpre+meanpost;
      //  cout<<meanpre<<","<<meanpost<<" ";
    //    perms[ip]=meanpre; perms[ip+size]=meanpost;
  //  }
    for (int i=0; i<numbl; i++)
        cout<<meanperm/2<<","; // (size*2)
    cout<<endl;
 /*   for (int i=0; i<numbl; i++)                         // pre-stim / rs: pre-stim > rs, pre-stim < rs
        if (eyescondition[currentsubj][i]==eyescond)
             cout<<getpvalue(size*2,prestimsynchs[currentsubj][i],true)<<","<<getpvalue(size*2,prestimsynchs[currentsubj][i],false)<<" ";
    cout<<endl;
    for (int i=0; i<numbl; i++)
        if (eyescondition[currentsubj][i]==eyescond)    // post-stim / rs: post-stim > rs, post-stim < rs
             cout<<getpvalue(size*2,poststimsynchs[currentsubj][i],true)<<","<<getpvalue(size*2,poststimsynchs[currentsubj][i],false)<<" ";
    cout<<endl; */

 //   ui->widget->graph(3)->setData(arrc.xc, arrc.of1);
 //   ui->widget->replot();

   // for (int i=0; i<10; i++)
     //   if (i<9)
     //       cout<<setprecision(3)<<(meanpresynchs[0]+meanpostsynchs[0])/2.0<<',';
     //   else
     //       cout<<setprecision(3)<<(meanpresynchs[0]+meanpostsynchs[0])/2.0;
      //  cout<<setprecision(3)<<(meanpreplvs[0]+meanpostplvs[0])/2.0<<',';
   // cout<<setprecision(3)<<meanpreplvs[0]<<" "<<meanpostplvs[0];

    //cout<<setprecision(3)<<"Mean pre / post PLV, average: "<<meanpreplvs[0]<<"  "<<meanpostplvs[0]<<",  "<<(meanpreplvs[0]+meanpostplvs[0])/2.0<<endl;
    //cout<<setprecision(3)<<"Mean pre / post Synchrony, average: "<<meanpresynchs[0]<<"  "<<meanpostsynchs[0]<<",  "<<(meanpresynchs[0]+meanpostsynchs[0])/2.0<<endl;

 /*   zerophaseinit(lcutoff, hcutoff, butterord, hnt->srfr);
    int nump=3;
    int p = hnt->srfr/hnt->osfr;
    numbl=100; int t=p;//2*p
    len=300;
    for (int k=0; k<numbl; k++)
    {
        preinstphasediffs[k][0]=getphaseoninterval(k*len,2*t);
        postinstphasediffs[k][0]=getphaseoninterval((k+1)*len,2*t);
        for (int i=0; i<nump; i++)
        {
            preinstphasediffs[k][i+1]=preinstphasediffs[k][0]-getphaseoninterval(k*len+t+i*p,2*p);
            postinstphasediffs[k][i+1]=postinstphasediffs[k][0]-getphaseoninterval((k+1)*len+t+i*p,2*p);
        }
    } */
   /* for (int k=0; k<numbl; k++)
    {
        cout<<"Trial "<<k+1<<":  ";
        for (int j=1; j<nump+1; j++)
            cout<<fixed<<setprecision(3)<<instphasediffs[k][j]<<"  ";
        cout<<endl;
    } */
  /*  for (int i=0; i<nump; i++)
    {
        complex<double> postsumexp(0.0,0.0);
        complex<double> presumexp(0.0,0.0);
        for (int k=0; k<numbl; k++)
        {
            presumexp+=exp(I*preinstphasediffs[k][i+1]);
            postsumexp+=exp(I*postinstphasediffs[k][i+1]);
        }
        preppi[i]=fabs(presumexp)/numbl;
        postppi[i]=fabs(postsumexp)/numbl;
    }
   // cout<<endl;
    //cout<<"Mean pre/post PPIs for 1P, 2P, 3P : ";
    //for (int i=0; i<nump; i++)
     //   if (i==nump-1)
     //       cout<<fixed<<setprecision(3)<<postppi[i]-preppi[i];
      //  else
 //   for (int i=0; i<10; i++)  // [0] - 1P, [1] - 2P, [2] - 3P
 //   if (i<9)
 //       cout<<(preppi[1]+postppi[1])/2.0<<',';
 //   else
 //      cout<<(preppi[1]+postppi[1])/2.0;
    cout<<endl; */
}

void plotwindow::analyzeallrestingstate()
{
    QString filename;

    for (int i=1; i<23; i++)
    {
      //  cout<<"Subject "<<i<<endl;
        currentsubj = i-1; // because eyes condition array starts from 0 :)
        if (i<10)
            filename = "C:/EEGdat/Subj0"+QString::number(i)+"/restingstateEO.dat";
        else
            filename = "C:/EEGdat/Subj"+QString::number(i)+"/restingstateEO.dat";
        loadrestingstate(filename);
       // analyzerestingstate();
        analyzerspower();
    }
}

void plotwindow::readstimplvs()
{
    QFile inputFile1("C:/EEGdat/Analysis/preplvs.dat");
    inputFile1.open(QIODevice::ReadOnly);
    QTextStream fin1(&inputFile1);
    QFile inputFile2("C:/EEGdat/Analysis/postplvs.dat");
    inputFile2.open(QIODevice::ReadOnly);
    QTextStream fin2(&inputFile2);

    QString line1 = fin1.readLine();
    QString line2 = fin2.readLine();
    QStringList sl;
    int leng=0; int subjn=20;

    for (int i=0; i<subjn; i++)
    {
        line1 = fin1.readLine();
        line2 = fin2.readLine();
        sl = line1.split(QRegExp("\\,"), QString::SkipEmptyParts);
        leng=sl.size();
        for (int k=0; k < leng; k++)
            prestimplvs[i][k]=sl[k].toDouble();
        sl = line2.split(QRegExp("\\,"), QString::SkipEmptyParts);
        leng=sl.size();
        for (int k=0; k < leng; k++)
            poststimplvs[i][k]=sl[k].toDouble();
    }
    inputFile1.close();
    inputFile2.close();
  /*  for (int i=0; i<subjn; i++)
    {
        for (int k=0; k < 10; k++)
            cout<<prestimplvs[i][k]<<",";
        cout<<endl;
    }
    cout<<endl;
    for (int i=0; i<subjn; i++)
    {
        for (int k=0; k < 10; k++)
            cout<<poststimplvs[i][k]<<",";
        cout<<endl;
    } */
}

void plotwindow::readstimsynchs()
{
    QFile inputFile1("C:/EEGdat/Analysis/presynchs.dat");
    inputFile1.open(QIODevice::ReadOnly);
    QTextStream fin1(&inputFile1);
    QFile inputFile2("C:/EEGdat/Analysis/postsynchs.dat");
    inputFile2.open(QIODevice::ReadOnly);
    QTextStream fin2(&inputFile2);

    QString line1 = fin1.readLine();
    QString line2 = fin2.readLine();
    QStringList sl;
    int leng=0; int subjn=20;

    for (int i=0; i<subjn; i++)
    {
        line1 = fin1.readLine();
        line2 = fin2.readLine();
        sl = line1.split(QRegExp("\\,"), QString::SkipEmptyParts);
        leng=sl.size();
        for (int k=0; k < leng; k++)
            prestimsynchs[i][k]=sl[k].toDouble();
        sl = line2.split(QRegExp("\\,"), QString::SkipEmptyParts);
        leng=sl.size();
        for (int k=0; k < leng; k++)
            poststimsynchs[i][k]=sl[k].toDouble();
    }
    inputFile1.close();
    inputFile2.close();

   /*   for (int i=0; i<subjn; i++)
      {
          for (int k=0; k < 10; k++)
              cout<<prestimsynchs[i][k]<<",";
          cout<<endl;
      }
      cout<<endl;
      for (int i=0; i<subjn; i++)
      {
          for (int k=0; k < 10; k++)
              cout<<poststimsynchs[i][k]<<",";
          cout<<endl;
      } */
}

int plotwindow::getchnind(QString str)
{
    int t=0;
    for (int i=0; i<31; i++)
        if (channelsinfo[i].chnname==str)
            t=i;
    return t;
}

void plotwindow::addchannelvalues(channinfo* chnin, int chnnums, double stdv)
{
    int indexchn; double val;
    for (int i=0; i<chnnums; i++)
    {
        indexchn = getchnind(chnin[i].chnname);
        val = abs(chnin[i].value/stdv); // abs()
        channelsinfo[indexchn].value+=val;
    }
}

 double plotwindow::getstdvalue(double *arr, int sz)
 {
     vector<double> v(sz);
     for (int i=0; i<sz; i++)
         v[i]=arr[i];
     double sum = std::accumulate(v.begin(), v.end(), 0.0);
     double mean = sum / v.size();

     std::vector<double> diff(v.size());
     std::transform(v.begin(), v.end(), diff.begin(),std::bind2nd(std::minus<double>(), mean));
     double sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
     double stdev = std::sqrt(sq_sum / v.size());

    // qDebug()<<mean;
   //  qDebug()<<stdev;

     return stdev;
 }

void plotwindow::getmeantopo()
{
    QFile inputFile("E:/EEGdat/analysis/CSPtopo/allchannels.txt");
    inputFile.open(QIODevice::ReadOnly);
    QTextStream fin(&inputFile);
    QStringList sl;
    QString st;
    for (int i=0; i<31; i++)
    {
        st = fin.readLine();
        sl = st.split(QRegExp("\\s+"), QString::SkipEmptyParts);
        channelsinfo[i].chnname=sl[0];
        channelsinfo[i].value=0;
    }
    inputFile.close();

    int totalsubjs=20;
    int skiplines=0;
    int blocknum=-1;
    int chnum=0;

    for (int j=1; j<totalsubjs; j++)
    {
       // if ((j==13))
       //     continue;
        double arrchnvalue[30];  // getting values of channels
        blocknum=noisearr[j-1];
        QFile valueFile("E:/EEGdat/Analysis/CSPtopo/subj"+QString::number(j)+"_block"+QString::number(blocknum)+".dat");
        valueFile.open(QIODevice::ReadOnly);
        QTextStream finval(&valueFile);
        skiplines = cspcomparr2[j-1][blocknum-1];
       // for (int i=0; i<skiplines-1; i++) // for CSPpre
       //     st = finval.readLine();
       // st = finval.readLine();
        {                                  // for CSPpost
            st = finval.readLine();
            sl = st.split(QRegExp("\\,"), QString::SkipEmptyParts);
            chnum = sl.length();
            for (int i=0; i<chnum-skiplines-1; i++)
                st = finval.readLine();
            st = finval.readLine();
        }
        sl = st.split(QRegExp("\\,"), QString::SkipEmptyParts);
        chnum = sl.length();
        for (int i=0; i<chnum; i++)
            arrchnvalue[i]=sl[i].toDouble();
        valueFile.close();

        double stdv = getstdvalue(arrchnvalue,chnum);

        channinfo tempchn[chnum];   // getting list of channels
        inputFile.setFileName("E:/EEGdat/Analysis/CSPtopo/subj"+QString::number(j)+"_"+QString::number(blocknum)+".txt");
        inputFile.open(QIODevice::ReadOnly);
        QTextStream fint(&inputFile);
        for (int i=0; i<chnum; i++)
        {
            st = fint.readLine();
            sl = st.split(QRegExp("\\s+"), QString::SkipEmptyParts);
            tempchn[i].chnname=sl[0];
            tempchn[i].value=arrchnvalue[i];
        }
        inputFile.close();

      //  for (int i=0; i<chnum; i++)
      //      cout<<tempchn[i].chnname.toStdString()<<" "<<tempchn[i].value<<endl;
        addchannelvalues(tempchn,chnum,stdv);

    }

  //  for (int i=0; i<31; i++)
  //      cout<<channelsinfo[i].chnname.toStdString()<<" "<<channelsinfo[i].value/(totalsubjs-2)<<endl;

    QFile outputFile("E:/EEGdat/Analysis/CSPtopo/meantopo_postnoise.dat");
    outputFile.open(QIODevice::WriteOnly);
    QTextStream fout(&outputFile);
    for (int i=0; i<31; i++)
    {
        fout << channelsinfo[i].value/totalsubjs;
        if (i<30)
            fout<<",";
    }
    outputFile.close();
}

/* ui processing */

void plotwindow::on_pushButton_clicked()
{
    offlinedata=true;    
    butterord=1;
    lcutoff=8;
    hcutoff=10;
    zerophaseinit(lcutoff, hcutoff, butterord, hnt->srfr);
    if (randomstims)
        setrandomrelations();
    else
        setfixedrelations();
    sequencestim();
}

void plotwindow::on_checkBox_clicked()
{
    if (!ui->checkBox->isChecked())
    {
        hnt->osfr=hnt->defosfr;
        ui->doubleSpinBox->setValue(hnt->defosfr);
    };
}

void plotwindow::on_checkBox_2_clicked()
{
    if (!ui->checkBox_2->isChecked())
    {
        ui->doubleSpinBox_2->setEnabled(false);
        hnt->twooscil = false;
    }
    else
    {
        ui->doubleSpinBox_2->setEnabled(true);
        hnt->twooscil = true;
    }    
}

void plotwindow::on_radioButton_2_clicked()
{    
    randnoisestim=false;
    randomstims=false;
    hnt->inphase=false;
}

void plotwindow::on_radioButton_clicked()
{
    randnoisestim=false;
    randomstims=false;
    hnt->inphase=true;
}

void plotwindow::on_checkBox_3_clicked()
{
    if (ui->checkBox_3->isChecked())
        usefiltering=true;
    else
        usefiltering=false;
}

void plotwindow::on_widget_destroyed()
{
    delete hnt;
}

void plotwindow::on_pushButton_2_clicked()
{
    clearstim();
    clearfiltered();
}

void plotwindow::on_pushButton_3_clicked()
{
    QString fname=QFileDialog::getSaveFileName(this,tr("Save File as"),"C://","Data file (*.dat);;All files (*.*)");
    if (fname!="")
    {
        savedatatofile(fname);
        mw->savelogtofile(fname+"log");
    }
}

void plotwindow::on_pushButton_4_clicked()
{
    slSettings();
}

void plotwindow::on_checkBox_4_clicked()
{
    if (ui->checkBox_4->isChecked())
        addmode=true;
    else
        addmode=false;
}

void plotwindow::on_doubleSpinBox_3_valueChanged(double arg1)
{
    hnt->imstprop=arg1;
}

void plotwindow::on_spinBox_7_valueChanged(int arg1)
{
    if (!tim->isActive())
        hnt->numst=arg1;
}

void plotwindow::on_spinBox_5_valueChanged(int arg1)
{
    hnt->imlength=arg1;
    if (!cleardata)
    {
       // hnt->lfilt=arg1/10;
       // ui->spinBox_2->setValue(hnt->lfilt);
        autoregdegree=arg1/2;
    }
}

void plotwindow::on_spinBox_valueChanged(int arg1)
{
    hnt->posstim=arg1;
}

void plotwindow::on_spinBox_2_valueChanged(int arg1)
{
    hnt->lfilt=arg1;
}

void plotwindow::on_spinBox_3_valueChanged(int arg1)
{
    hnt->phase=arg1;
}

void plotwindow::on_doubleSpinBox_valueChanged(double arg1)
{
    hnt->osfr=arg1;
}

void plotwindow::on_spinBox_6_valueChanged(int arg1)
{
    hnt->stampl=arg1;
}

void plotwindow::on_spinBox_4_valueChanged(int arg1)
{
    hnt->noise=arg1;
}

void plotwindow::on_doubleSpinBox_2_valueChanged(double arg1)
{
    hnt->thfr=arg1;
}

void plotwindow::on_spinBox_8_valueChanged(int arg1)
{
    hnt->srfr=arg1;
}

void plotwindow::on_spinBox_5_editingFinished()
{
    //autoregdegree=ui->spinBox_5->value()-10;
}

void plotwindow::on_pushButton_6_clicked()
{
    this->showMinimized();
}

void plotwindow::on_radioButton_3_clicked()
{
    randnoisestim=true;
    randphasestim=false;
    randomstims=false;
}

void plotwindow::on_radioButton_4_clicked()
{
    randphasestim=true;
    randnoisestim=false;
    randomstims=false;
}

void plotwindow::on_radioButton_5_clicked()
{
    randomstims=true;
}

int plotwindow::findoptardegree()
{
    double minar_err = 1000;
    int optar_degree = 100;   
    int lb = 40; // hnt->imlength / 2 - hnt->imlength / 10;
    int hb = 90; // hnt->imlength / 2 + hnt->imlength / 10;
    for (int k=lb; k<hb; k++)
    {
        totalar_error = 0;
        autoregdegree = k;
        for (int i=0; i<nstimforopt; i++)
        {
            runautoreg(hnt->posstim+(hnt->imlength+hnt->stlength)*i);
            totalar_error+=ar_error;
        }
       // qDebug()<<totalar_error;
        if (totalar_error<minar_err)
        {
            minar_err=totalar_error;
            optar_degree=k;
        }
    }

    autoregdegree = optar_degree;

   // for (int i=0; i<np; i++)
   //     runautoreg(hnt->posstim+(hnt->imlength+hnt->stlength)*i);
   // ui->widget->graph(2)->setData(arrc.xc, arrc.of2);

    return optar_degree;
}

void plotwindow::findoptzcparams()
{   
    int optlfb = 0;
    int optrfb = 0;
    int optord = 0;
    double averageacc = 0;
    double maxacc = 0;
    int pos = hnt->posstim;

    offlinedata = true;
    for (int q=2; q<4; q++)
    for (int t=6; t<9; t++)
    for (int k=12; k<17; k++)
    {
        zerophaseinit(t, k, q, hnt->srfr);
        for (int i=0; i<nstimforopt; i++)
        {
            zerophasefiltzc(pos+(hnt->imlength+hnt->stlength)*i-hnt->imlength,hnt->imlength);
            zerocrossingstim(pos+(hnt->imlength+hnt->stlength)*i,hnt->imlength);

            double diffZC = zcdifference(pos+(hnt->imlength+hnt->stlength)*i);
            double proczc = 0;

            offlinevaluation(pos,i,3);

            if (hnt->inphase)
                proczc=hnt->optphasediff/diffZC;
            else
                proczc=diffZC/hnt->optphasediff;

            averageacc+=proczc;
        }
        averageacc/=nstimforopt;
      //  qDebug()<<averageacc;
        if (averageacc>maxacc)
        {
            maxacc=averageacc;
            optlfb=t; optrfb=k; optord=q;
        }
    }

    lcutoff = optlfb;
    hcutoff = optrfb;
    butterord = optord;

  //  qDebug()<<optord;
  //  qDebug()<<optlfb;
  //  qDebug()<<optrfb;

  /*  zerophaseinit(optlfb, optrfb, optord, hnt->srfr);
    for (int i=0; i<np; i++)
    {
        zerophasefiltzc(pos+(hnt->imlength+hnt->stlength)*i-hnt->imlength,hnt->imlength);
        zerocrossingstim(pos+(hnt->imlength+hnt->stlength)*i,hnt->imlength);

        double diffZC = zcdifference(pos+(hnt->imlength+hnt->stlength)*i);
        double proczc = 0;

        offlinevaluation(pos,i,3);

        if (hnt->inphase)
            proczc=hnt->optphasediff/diffZC;
        else
            proczc=diffZC/hnt->optphasediff;

        averageacc+=proczc;
    }
    averageacc/=np;
   // qDebug()<<averageacc;

    ui->widget->graph(5)->setData(arrc.xc, arrc.of4); */
}

void plotwindow::saveoptparams()
{
    QFile outputFile("C:/optparams.txt");
    outputFile.open(QIODevice::Append);
    QTextStream fout(&outputFile);
    fout << butterord << " " << lcutoff << " " << hcutoff << " " << autoregdegree << endl;
    outputFile.close();
}

void plotwindow::dophasepred()
{
    for (int i=0; i<4; i++)
    {
        hnt->imlength=100+50*i;
        ui->spinBox_5->setValue(hnt->imlength);
        on_radioButton_clicked();

        clearstim();
        clearfiltered();
        findoptzcparams();  // ZC-optimization
        findoptardegree();  // AR-optimization
        saveoptparams();
        clearstim();
        clearfiltered();

        on_pushButton_clicked();
        on_radioButton_2_clicked();
        on_pushButton_clicked();
        on_radioButton_5_clicked();
        on_pushButton_clicked();
        ui->widget->replot();
    }
}

void plotwindow::read2ndexpdata(QString fname)
{
    QFile inputFile(fname);
    inputFile.open(QIODevice::ReadOnly);
    QTextStream fin(&inputFile);
    QStringList sl;

    QString line = fin.readLine();
    sl = line.split(QRegExp("\\s+"), QString::SkipEmptyParts);
    trialnums = sl[1].toInt();
    line = fin.readLine();
    int leng=0;

    if (line.length()>0)
    {
        sl = line.split(QRegExp("\\s+"), QString::SkipEmptyParts);
        leng = sl.size();
        for (int k=0; k < leng; k++)
            arrc.amp1[k]=sl[k].toDouble();
    }
    ui->widget->graph(0)->setData(arrc.xc, arrc.amp1);
    ui->widget->replot();
    inputFile.close();

}

void plotwindow::analyzersspectrapoz(int pos1,  int pos2, int numtr)
{
    lcutoff=5; hcutoff=41; butterord=2;

    for (int i=lcutoff; i<hcutoff; i++)
    {
        spectECpre[i-lcutoff]=0;
        spectECpost[i-lcutoff]=0;
    }
    rsalphapre=0;
    rsalphapost=0;

    zerophaseinit(lcutoff, hcutoff, butterord, hnt->srfr);
    int length = 500;
    int zp = 12;       // 12 - for 1000 ms
    int start = 0;
    int len = 512;
    double* frampl;
    int fr = (int)round(iafrs); // IAF

    for (int j=0; j<numtr; j++)
    {
        start=pos1+j*length;
        Complex t[len];
        zerophasefilt(start,length);
        for (int i=start; i<start+length; i++)
            t[i-start].real(arrc.amp1[i]);

        for (int i=0; i<len; i++)
            t[i].imag(0);
        for (int i=0; i<zp; i++)  // zero-padding
            t[length+i].real(0);

        cdata = CArray(t,len);
        hnt->fft(cdata);
        frampl = new double[len];
        for (int i=0; i<len; i++)
        {
            frampl[i]=2*sqrt(cdata[i].real()*cdata[i].real()+cdata[i].imag()*cdata[i].imag());
            frampl[i]/=len;
        }

        for (int i=lcutoff; i<hcutoff; i++)
            spectECpre[i-lcutoff]+=frampl[i];
        rsalphapre+=frampl[fr-1]+frampl[fr]+frampl[fr+1];

        start=pos2+j*length;
        zerophasefiltzc(start,length);
        for (int i=start; i<start+length; i++)
            t[i-start].real(arrc.amp2[i]);

        for (int i=0; i<len; i++)
            t[i].imag(0);
        for (int i=0; i<zp; i++)  // zero-padding
            t[length+i].real(0);

        cdata = CArray(t,len);
        hnt->fft(cdata);
        frampl = new double[len];
        for (int i=0; i<len; i++)
        {
            frampl[i]=2*sqrt(cdata[i].real()*cdata[i].real()+cdata[i].imag()*cdata[i].imag());
            frampl[i]/=len;
        }

        for (int i=lcutoff; i<hcutoff; i++)
            spectECpost[i-lcutoff]+=frampl[i];
      //  cout<<fixed<<setprecision(2)<<arrpowerpost[j]-arrpowerpre[j]<<" ";
        rsalphapost+=frampl[fr-1]+frampl[fr]+frampl[fr+1];;
    }

    for (int i=lcutoff; i<hcutoff; i++)
    {
        spectECpre[i-lcutoff]/=numtr;
        spectECpost[i-lcutoff]/=numtr;
    }
    rsalphapre/=numtr;
    rsalphapost/=numtr;
}

void plotwindow::analyzespectra(int numtr, int timew)
{
    lcutoff=5; hcutoff=41; butterord=2;

    for (int i=lcutoff; i<hcutoff; i++)
    {
        spectECpre[i-lcutoff]=0;
        spectECpost[i-lcutoff]=0;
    }

    if (!spatial_fitering)
        zerophaseinit(lcutoff, hcutoff, butterord, hnt->srfr);

    int length = 750;
    int zp;
    int start = 0;
    int sh1;
    int sh2 = 60;  // shift from start of post-stim
    int len = 512;
    int lentr;

    if (timew==1000)
    {
        lentr = 500;
        sh1 = 150;
        zp = 12;
    }
    else if (timew==500)
    {
        lentr = 250;
        sh1 = 400;
        zp = 262;
    }

    double* frampl;

    for (int j=0; j<numtr; j++)
    {
        start=j*(length*2);
        Complex t[len];
        if (!spatial_fitering)
            zerophasefilt(start,length);
        for (int i=start+sh1; i<start+sh1+lentr; i++)
            t[i-start-sh1].real(arrc.amp1[i]);

        for (int i=0; i<len; i++)
            t[i].imag(0);
        for (int i=0; i<zp; i++)  // zero-padding
            t[lentr+i].real(0);

        cdata = CArray(t,len);
        hnt->fft(cdata);
        frampl = new double[len];
        for (int i=0; i<len; i++)
        {
            frampl[i]=2*sqrt(cdata[i].real()*cdata[i].real()+cdata[i].imag()*cdata[i].imag());
            frampl[i]/=len;
        }

        for (int i=lcutoff; i<hcutoff; i++)
            spectECpre[i-lcutoff]+=frampl[i];

        start=length+j*(length*2);
        if (!spatial_fitering)
            zerophasefilt(start,length);
        for (int i=start+sh2; i<start+sh2+lentr; i++)
            t[i-start-sh2].real(arrc.amp1[i]);

        for (int i=0; i<len; i++)
            t[i].imag(0);
        for (int i=0; i<zp; i++)  // zero-padding
            t[lentr+i].real(0);

        cdata = CArray(t,len);
        hnt->fft(cdata);
        frampl = new double[len];
        for (int i=0; i<len; i++)
        {
            frampl[i]=2*sqrt(cdata[i].real()*cdata[i].real()+cdata[i].imag()*cdata[i].imag());
            frampl[i]/=len;
        }

        for (int i=lcutoff; i<hcutoff; i++)
            spectECpost[i-lcutoff]+=frampl[i];
      //  cout<<fixed<<setprecision(2)<<arrpowerpost[j]-arrpowerpre[j]<<" ";
    }

    for (int i=lcutoff; i<hcutoff; i++)
    {
        spectECpre[i-lcutoff]/=numtr;
        spectECpost[i-lcutoff]/=numtr;
    }

}

void plotwindow::extractpozrsspectra()
{
    QString fname, foldername;
    int subjnums = 2;
    for (int k=1; k<subjnums; k++)
    {
        if (k<10)
            foldername = "E:/EEGdat/Subj0"+QString::number(k)+"/";
        else
            foldername = "E:/EEGdat/Subj"+QString::number(k)+"/";

        for (int i=1; i<4; i++)
        {
            int t=3;
            if ((k==13))
                t=2;
            for (int j=1; j<t; j++)
            {
                fname=foldername+QString::number(j)+"rest"+QString::number(i)+".dat";
                QFile inpFile(fname);                                        // read data from file
                inpFile.open(QIODevice::ReadOnly);
                QTextStream finp(&inpFile);
                QString line; QStringList sl;
                QFile outputFile(foldername+"all"+QString::number(j)+"rest"+QString::number(i)+".dat");
                outputFile.open(QIODevice::WriteOnly);
                QTextStream fout(&outputFile);
                for (int r=0; r<31; r++)
                {
                    //finp.readLine();
                    line = finp.readLine();
                    sl = line.split(QRegExp("\\s+"), QString::SkipEmptyParts);
                    int length = sl.length();
                    for (int c=1; c<length; c++)
                    if (c<length-1)
                        fout << sl[c] <<",";
                    else
                        fout << sl[c];
                    fout << endl;
                }
                outputFile.close();
            }
        }
    }
}

void plotwindow::analyzersspectra()
{
    QString foldername, fname1, fname2;
    int stpos1, stpos2;
    stpos1 = 100000;
    stpos2 = 3000;

    if (currentsubj<10)
        foldername = "E:/EEGdat/Subj0"+QString::number(currentsubj)+"/";
    else
        foldername = "E:/EEGdat/Subj"+QString::number(currentsubj)+"/";
    switch (currentblock) {
     case 1:
        fname1 = foldername+"POz1rest1.dat";
        fname2 = foldername+"POz1rest2.dat";
        iafrs=alphaarr[currentsubj-1][0];
     break;
     case 2:
        fname1 = foldername+"POz1rest2.dat";
        fname2 = foldername+"POz1rest3.dat";
        iafrs=alphaarr[currentsubj-1][0];
     break;
     case 3:
        fname1 = foldername+"POz2rest1.dat";
        fname2 = foldername+"POz2rest2.dat";
        iafrs=alphaarr[currentsubj-1][1];
     break;
     case 4:
        fname1 = foldername+"POz2rest2.dat";
        fname2 = foldername+"POz2rest3.dat";
        iafrs=alphaarr[currentsubj-1][1];
     break;
    }

    clean();
    QFile inpFile1(fname1);
    QFile inpFile2(fname2);
    inpFile1.open(QIODevice::ReadOnly);
    inpFile2.open(QIODevice::ReadOnly);
    QTextStream finp1(&inpFile1);
    QTextStream finp2(&inpFile2);
    QString line; QStringList sl1,sl2;
    line = finp1.readLine();
    sl1 = line.split(QRegExp("\\s+"), QString::SkipEmptyParts);
    line = finp2.readLine();
    sl2 = line.split(QRegExp("\\s+"), QString::SkipEmptyParts);
    inpFile1.close(); inpFile2.close();
    int length1 = sl1.length();
    int length2 = sl2.length();
    for (int i=0; i<length1; i++)
        arrc.amp1[i]=sl1[i].toDouble();
    for (int i=0; i<length2; i++)
        arrc.amp2[i]=sl2[i].toDouble();
    analyzersspectrapoz(stpos1, stpos2, 200);

}

void plotwindow::mergepowers()
{
    int subjnums = 20;
        QString fres1, fres2, fres3;

    for (int i=1; i<subjnums; i++)
    {
        cout<<"Subject "<<i<<endl;

        if (i<10)
        {
            fres1 = "E:/EEGdat/analysis/results/500plv/subj0"+QString::number(i)+"_plvPOz.dat";
            fres2 = "E:/EEGdat/analysis/results/500plv/subj0"+QString::number(i)+"_plvCSPpre.dat";
            fres3 = "E:/EEGdat/analysis/results/500plv/subj0"+QString::number(i)+"_plvCSPpost.dat";
        }
        else
        {
            fres1 = "E:/EEGdat/analysis/results/500plv/subj"+QString::number(i)+"_plvPOz.dat";
            fres2 = "E:/EEGdat/analysis/results/500plv/subj"+QString::number(i)+"_plvCSPpre.dat";
            fres3 = "E:/EEGdat/analysis/results/500plv/subj"+QString::number(i)+"_plvCSPpost.dat";
        }

        mergeresfiles(fres1,fres2,fres3);
    }
}

void plotwindow::mergepowerswithpoc()
{
    int subjnums = 20;
        QString fres1, fres2;

    for (int i=1; i<subjnums; i++)
    {
        cout<<"Subject "<<i<<endl;

        if (i<10)
        {
            fres1 = "E:/EEGdat/analysis/results/500ms/subj0"+QString::number(i)+".dat";
            fres2 = "E:/EEGdat/analysis/results/500poc/subj0"+QString::number(i)+".dat";
        }
        else
        {
            fres1 = "E:/EEGdat/analysis/results/500ms/subj"+QString::number(i)+".dat";
            fres2 = "E:/EEGdat/analysis/results/500poc/subj"+QString::number(i)+".dat";
        }

        mergeresfiles2(fres1,fres2);
    }
}

void plotwindow::analyze2ndplv(int numbl)
{
    for (int i=0; i<numbl; i++)
    {
        arrplvpre[i]=0;
        arrplvpost[i]=0;
    }
    lcutoff=5; hcutoff=41; butterord=2;
    zerophaseinit(lcutoff, hcutoff, butterord, hnt->srfr);
    int length = 750;
    int start = 0;
    int sh1 = 400; // shift from start of pre-stim
    int sh2 = 60;  // shift from start of post-stim
    int lentr = 250;
    double st1[lentr];
    double st2[lentr];

    for (int j=0; j<numbl; j++)
    {
        hnt->phase=qrand()%50;

        start=j*(length*2);
        if (!spatial_fitering)
            zerophasefilt(start,length);

        fillstim(start,length*2);

        for (int k=0; k<lentr; k++)
            st1[k]=arrc.of1[start+sh1+k];
        for (int k=0; k<lentr; k++)
            st2[k]=arrc.amp1[start+sh1+k];
        hnt->phsdif = phasediff(st1, st2, lentr);
        arrplvpre[j]=plv;

        start=length+j*(length*2);
        if (!spatial_fitering)
            zerophasefilt(start,length);

        for (int k=0; k<lentr; k++)
            st1[k]=arrc.of1[start+sh2+k];
        for (int k=0; k<lentr; k++)
            st2[k]=arrc.amp1[start+sh2+k];
        hnt->phsdif = phasediff(st1, st2, lentr);
        arrplvpost[j]=plv;

       // cout<<fixed<<setprecision(3)<<arrplvpre[j]<<" "<<arrplvpost[j]<<endl;
    }
  //  ui->widget->graph(0)->setData(arrc.xc, arrc.amp1);
  //  ui->widget->graph(3)->setData(arrc.xc, arrc.of1);
  //  ui->widget->replot();
}

void plotwindow::fillblockparts(int trnums)
{
    for (int i=0; i<blockparts; i++)
    {
        prestimtrials[i]=0;
        poststimtrials[i]=0;
    }

    for (int i=0; i<trnums; i++)
    {        
        prestimtrials[i/40]+=arrpowerpre[i]; // arrplvpre
        poststimtrials[i/40]+=arrpowerpost[i];
    }

    int rt;
    if (trnums<200)
        rt = trnums/40;
    else
        rt = 5;

    for (int i=0; i<rt; i++)
    {
        prestimtrials[i]/=40;
        poststimtrials[i]/=40;
    }
    if (rt<5)
    {
        prestimtrials[rt]/=trnums-rt*40;
        poststimtrials[rt]/=trnums-rt*40;
    }
}

void plotwindow::analyze2ndexp()
{
    cleareeg(0,505850);

    int subjnums = 20;
    QString filename;
    QString fresname;

    for (int i=1; i<subjnums; i++)
    {                
        cout<<"Subject "<<i<<endl;

        int t=5;
        if (i==13)
            t=3;

        for (int j=1; j<t; j++)
        {            
            if (i<10)
            {
                filename = "E:/EEGdat/Subj0"+QString::number(i)+"/block"+QString::number(j)+".dat";
                fresname = "E:/EEGdat/analysis/results/subj0"+QString::number(i)+"_plvCSPpost.dat";
            }
            else
            {
                filename = "E:/EEGdat/Subj"+QString::number(i)+"/block"+QString::number(j)+".dat";
                fresname = "E:/EEGdat/analysis/results/subj"+QString::number(i)+"_plvCSPpost.dat";
            }

            fsubjshort = filename.left(17)+"rawPOz"+QString::number(j)+".dat";
         //   fsubjshort = filename.left(17)+"cspPost"+QString::number(j)+".dat";
         //   spatial_fitering=true;
            read2ndexpdata(fsubjshort);
           // trialnums = extractblock(filename);

            cout<<trialnums<<" ";                        
            if (j<3)
                hnt->osfr = alphaarr[i-1][0];
            else
                hnt->osfr = alphaarr[i-1][1];

           // savesingledatatofile(fsubjshort,hnt->osfr);

            currentblock=j; currentsubj=i;
          //  getcsptopo2(filename,trialnums);
          //  applycsp2ndexp(filename,trialnums,cspcomparr2[i-1][j-1],"S2");
          //  QString fcsp = filename.left(17)+"cspPost"+QString::number(j)+".dat";
          //  savecspdatatofile(fcsp,hnt->osfr);

            analyze2ndplv(trialnums);
          //  analyzersspectra();
          //  analyzeblockpower(trialnums,1000);
            meanpowerfromopregion2(filename,trialnums,j-1,i-1);
          //  fillblockparts(trialnums);
            saveresulttofile(j-1,i-1,fresname);
          //  savepartstofiles(j-1,i-1,"CSPpre500");
          //  analyzespectra(trialnums,1000);
          //  meanspectfromopregion2(filename,trialnums,j-1,i-1);
          //  savespecttofiles(j-1,i-1,"sp1POC500");

          //  ui->widget->graph(1)->setData(arrc.xc, arrc.amp2);
        }
        cout<<endl;
    }
}

void plotwindow::savepartstofiles(int blocki, int currentsubj, QString suffix)
{
    QFile outresFile1;

    if (blocksarr[currentsubj][blocki]==1)
        outresFile1.setFileName("E:/EEGdat/analysis/"+suffix+"_pre-in.dat");
    else if (blocksarr[currentsubj][blocki]==-1)
        outresFile1.setFileName("E:/EEGdat/analysis/"+suffix+"_pre-anti.dat");
    else if (blocksarr[currentsubj][blocki]==3)
        outresFile1.setFileName("E:/EEGdat/analysis/"+suffix+"_pre-rand.dat");
    else
        outresFile1.setFileName("E:/EEGdat/analysis/"+suffix+"_pre-noise.dat");

    if (currentsubj==0)
        outresFile1.open(QIODevice::WriteOnly);
    else
        outresFile1.open(QIODevice::Append);

    if (((currentsubj==0) && (blocki==2)) || ((currentsubj==6) && (blocki==1)) || ((currentsubj==15) && (blocki==3)) || ((currentsubj==16) && (blocki==1)))
    {
        for (int i=0; i<blockparts; i++) // trialnums
        {
            prestimtrials[i]=0;
            poststimtrials[i]=0;
        }
        rsalphapre=0; rsalphapost=0;
    }

    QTextStream fresout1(&outresFile1);

    for (int i=0; i<blockparts; i++) // trialnums
    {
        fresout1<<prestimtrials[i];// arrpowerpre[i];
        if (i<blockparts-1)
            fresout1<<",";
    }
   // fresout1<<rsalphapre;

    fresout1<<endl;
    outresFile1.close();

    if (blocksarr[currentsubj][blocki]==1)
        outresFile1.setFileName("E:/EEGdat/analysis/"+suffix+"_post-in.dat");
    else if (blocksarr[currentsubj][blocki]==-1)
        outresFile1.setFileName("E:/EEGdat/analysis/"+suffix+"_post-anti.dat");
    else if (blocksarr[currentsubj][blocki]==3)
        outresFile1.setFileName("E:/EEGdat/analysis/"+suffix+"_post-rand.dat");
    else
        outresFile1.setFileName("E:/EEGdat/analysis/"+suffix+"_post-noise.dat");

    if (currentsubj==0)
        outresFile1.open(QIODevice::WriteOnly);
    else
        outresFile1.open(QIODevice::Append);

    QTextStream fresout2(&outresFile1);

    for (int i=0; i<blockparts; i++)
    {
        fresout2<<poststimtrials[i];
        if (i<blockparts-1)
            fresout2<<",";
    }
   // fresout2<<rsalphapost;

    fresout2<<endl;
    outresFile1.close();
}

void plotwindow::savespecttofiles(int blocki, int currentsbj, QString suffix)
{
    if (((currentsbj==0) && (blocki==2)) || ((currentsbj==6) && (blocki==1)) || ((currentsbj==15) && (blocki==3)) || ((currentsbj==16) && (blocki==1)))
        return;

    QFile outresFile1;

    int lborder = 5;
    int hborder = 41;

    if (blocksarr[currentsbj][blocki]==1)
        outresFile1.setFileName("E:/EEGdat/analysis/"+suffix+"_pre-in.dat");
    else if (blocksarr[currentsbj][blocki]==-1)
        outresFile1.setFileName("E:/EEGdat/analysis/"+suffix+"_pre-anti.dat");
    else if (blocksarr[currentsbj][blocki]==3)
        outresFile1.setFileName("E:/EEGdat/analysis/"+suffix+"_pre-rand.dat");
    else
        outresFile1.setFileName("E:/EEGdat/analysis/"+suffix+"_pre-noise.dat");

    QTextStream fresout1(&outresFile1);

    if ((currentsbj==0) || ((currentsbj==1) && (blocki==0)))
    {  
        outresFile1.open(QIODevice::WriteOnly);
        for (int i=lborder; i<hborder; i++)
        {
            fresout1<<i;
            if (i<hborder-1)
                fresout1<<",";
        }
        fresout1<<endl;
    }
    else
        outresFile1.open(QIODevice::Append);

    for (int i=lborder; i<hborder; i++)
    {
        fresout1<<spectECpre[i-lborder];
        if (i<hborder-1)
            fresout1<<",";
    }

    fresout1<<endl;
    outresFile1.close();

    if (blocksarr[currentsbj][blocki]==1)
        outresFile1.setFileName("E:/EEGdat/analysis/"+suffix+"_post-in.dat");
    else if (blocksarr[currentsbj][blocki]==-1)
        outresFile1.setFileName("E:/EEGdat/analysis/"+suffix+"_post-anti.dat");
    else if (blocksarr[currentsbj][blocki]==3)
        outresFile1.setFileName("E:/EEGdat/analysis/"+suffix+"_post-rand.dat");
    else
        outresFile1.setFileName("E:/EEGdat/analysis/"+suffix+"_post-noise.dat");

    QTextStream fresout2(&outresFile1);

    if ((currentsbj==0) || ((currentsbj==1) && (blocki==0)))
    {
        outresFile1.open(QIODevice::WriteOnly);
        for (int i=lborder; i<hborder; i++)
        {
            fresout2<<i;
            if (i<hborder-1)
                fresout2<<",";
        }
        fresout2<<endl;
    }
    else
        outresFile1.open(QIODevice::Append);

    for (int i=lborder; i<hborder; i++)
    {
        fresout2<<spectECpost[i-lborder];
        if (i<hborder-1)
            fresout2<<",";
    }

    fresout2<<endl;
    outresFile1.close();
}

void plotwindow::mergeresfiles(QString fres1, QString fres2, QString fres3)
{
    QFile inpFile1(fres1);
    QFile inpFile2(fres2);
    QFile inpFile3(fres3);
    inpFile1.open(QIODevice::ReadOnly);
    inpFile2.open(QIODevice::ReadOnly);
    inpFile3.open(QIODevice::ReadOnly);
    QFile outFile(fres1.left(fres1.length()-4)+"_n.dat");
    outFile.open(QIODevice::WriteOnly);

    QTextStream outFilest(&outFile);
    QTextStream inFile1st(&inpFile1);
    QTextStream inFile2st(&inpFile2);
    QTextStream inFile3st(&inpFile3);

    QStringList sl;
    QString line;

    line = inFile1st.readLine();

    while (line.length()>0)
    {
        outFilest<<line;
        line = inFile2st.readLine();
        sl = line.split(QRegExp("\\,"), QString::SkipEmptyParts);
        outFilest<<","<<sl[2]<<","<<sl[3];
        line = inFile3st.readLine();
        sl = line.split(QRegExp("\\,"), QString::SkipEmptyParts);
        outFilest<<","<<sl[2]<<","<<sl[3]<<endl;
        line = inFile1st.readLine();
    }

    outFile.close();
    inpFile1.close();
    inpFile2.close();
    inpFile3.close();
}

void plotwindow::mergeresfiles2(QString fres1, QString fres2)
{
    QFile inpFile1(fres1);
    QFile inpFile2(fres2);
    inpFile1.open(QIODevice::ReadOnly);
    inpFile2.open(QIODevice::ReadOnly);
    QFile outFile(fres1.left(fres1.length()-4)+"_n.dat");
    outFile.open(QIODevice::WriteOnly);

    QTextStream outFilest(&outFile);
    QTextStream inFile1st(&inpFile1);
    QTextStream inFile2st(&inpFile2);

    QStringList sl;
    QString line;

    line = inFile1st.readLine();

    while (line.length()>0)
    {
        outFilest<<line;
        line = inFile2st.readLine();
        sl = line.split(QRegExp("\\,"), QString::SkipEmptyParts);
        outFilest<<","<<sl[2]<<","<<sl[3]<<endl;
        line = inFile1st.readLine();
    }

    outFile.close();
    inpFile1.close();
    inpFile2.close();
}

void plotwindow::saveresulttofile(int blocki, int currentsubj, QString fresname)
{
    QFile outresFile(fresname);
    if (blocki==0)
        outresFile.open(QIODevice::WriteOnly);
    else
        outresFile.open(QIODevice::Append);
    QTextStream fresout(&outresFile);

    for (int i=0; i<trialnums; i++)
    {
        if (blocksarr[currentsubj][blocki]==1)
            fresout<<"in,";
        else if (blocksarr[currentsubj][blocki]==-1)
            fresout<<"anti,";
        else if (blocksarr[currentsubj][blocki]==3)
            fresout<<"rand,";
        else
            fresout<<"noise,";
        fresout<<i+1<<",";
        fresout<<arrplvpre[i]<<","<<arrplvpost[i]<<endl;
    }
    outresFile.close();
}

void plotwindow::on_pushButton_5_clicked() // DO BUTTON
{
  //  extractblock("C:/EEGdat/Subj01/block4.dat");
  //  analyzeblockpower(149);
  //  extractpozrsspectra();
    analyze2ndexp();
  //  parselrtcres("noise");

  //  mergepowerswithpoc();
  //  mergepowers();

  //  dophasepred();

  //  savescreen();

  //  sortreferencelist();

  //  getchanlabels();
  //  geteceopower();

  //  fsubjshort = "C:/EEGdat/Subj02/"; currentsubj = 2;
  //  loadstimfromfile(fsubjshort+"/rawdata.dat");
  //  loadstimfromfile(fsubjshort+"/csppre.dat");
  //  extractstimdata(fsubjshort,10,true);

  //  applycsp(fsubjshort+"10blocks.dat",10,1,"pre/post","S1");
  //  getcsptopo(fsubjshort+"10blocks.dat",10);
  //  analyzepower(1);
  //  analyzegammapow(1);

    //  analyzepower(10);
    // sortblocks("C:/EEGdat/Subj18/",10, true);
   // analyzestimdata(1,1,2,true);
    //  analyzeallsubjs();
   //  analyzeallsubjsmean();
   //  ppiforallsubjs();

  //  sortresultblocks("NCSP800(post)_pre_in");
 //   sortresultblocks("NCSP800(post)_post_in");
  //  sortresultblocks("NCSP800(post)_pre_anti");
 //   sortresultblocks("NCSP800(post)_post_anti");

  //   mergeresults("500_OCP");
  //  getmeantopo();

 //  flipsign();

  //   alignallresultblocks("MeanOP800");
  // combineresultblocks();

  //  for (int i=0; i<10; i++)
  //     phasepreservation(i,3);

  //  readstimsynchs();
    //currentsubj=8;
   // loadrestingstate("C:/EEGdat/Subj10/restingstateEC.dat");
   // analyzerestingstate();
   //   analyzerspower();

    // readstimplvs();
   // readstimsynchs();
   //  analyzeallrestingstate();
   // getpowerspectrumrs();

   //  testsyntdata();

    // fftphasepred();

  //  parseresfromfile();
   // parseoptparams();

    /*  float* inpd = new float[250];
    for (int i=0; i<250; i++)
        inpd[i]=arrc.amp1[hnt->posstim+i];
    morletwavelet(inpd,0,250); */
}
