#ifndef PLOT_H
#define PLOT_H

#include <QWidget>
#include <qcustomplot.h>
#include <hilbert.h>
#include <vector>
#include "mainwindow.h"
#include "omp.h"
#include "stdio.h"
#include <NIDAQmx.h>
#include <filt.h>
#include <Eigen/Dense>

namespace Ui {
class plot;
}

using namespace Eigen;

typedef std::vector<int> vectori;
typedef std::vector<double> vectord;

struct coords {
    QVector<double> xc,amp0,amp1,amp2,of1,of2,of3,of4,of5,of6;
    int lscale,lmax,posstim;   
};

struct freqband {
    int lf,hf;
};

struct channinfo {
    QString chnname;
    double value;
};

static const int shift=100;

class appconnect;
class MainWindow;
class Settings;

class plotwindow : public QWidget
{
    Q_OBJECT       

public:  
    int stepsPerPress, iti1, iti2, stimshift, receivedelay, avgdevhr, trialnums;
    int counter,stims,startpos,tscale,ampval,corrprednumb,recparts, autoregdegree, daqscalevalue, transdelay, flength, chnums, sampleblock, exch, sourcech;
    QString daqport;
    QString fsubjshort;
    bool start, offlinedata, addmode, addmodeon, estimation, adaptampl, addmoderec, hideardata, usefiltering, filterardata, carfilter, zerocr, phasepr;
    bool showprediction, estimpredstim, filtereddata, filterstims, correctdelay, maxentropy, leastsqr, regforstim, hronar, hronarrun, filtfir, fftpredict;
    bool randomstims, extractoptphase, spatial_fitering, use_mean_trial, randnoisestim, randphasestim;
    channinfo channelsinfo[31];
    double *meanpretrial;
    double rsalphapre, rsalphapost, iafrs;
    double *meanposttrial;
    double *prestimtrials;
    double *poststimtrials;
    double **arrres;
    double *allopttimes;
    double alltimesc;
    int opttimecount, nstimforopt;
    double *autoreginput;
    double *autoregcoeffs;    
    double *prestim;
    double *poststim;
    double *synchprestim;
    double *synchpoststim;
    double *meanpre_in;
    double *meanpost_in;
    double *meanpre_anti;
    double *meanpost_anti;
    double plv, averopttime, fsampr, flowpass;
    int* exlchannels;
    int* blockstarts;
    int* clearoptphase;
    double* arrphasediffhr;
    double* arrphasediffar;
    double* arrphasediffzc;
    double* arrphasedifffft;
    double* arrsynchr;
    double* arrsyncar;
    double* arrsynczc;
    double* arrsyncfft;
    double* arrpowerpre;
    double* arrpowerpost;
    double* arrplvpre;
    double* arrplvpost;
    double* meanpreplvs;
    double* meanpostplvs;
    double* meanpresynchs;
    double* meanpostsynchs;
    double* meanprepower;
    double* meanpostpower;
    double ar_error, totalar_error;
    int currentblock;
    bool recurbutter, zerobutter, imagingafterstim, skipdelayedata, readypredict, optstim, compavgdeviation, estimateclear, estimatenoisy, cleardata;
    vectord acoeffs, bcoeffs, ac, bc, b_f, a_f, b_s, a_s; // Filters for SSD
    int ssdscale, ssdshift;
    int butterord, lcutoff, hcutoff; // Butterworth filter
    Complex plvtotal;
    Complex fftarr[512];
    CArray cdata;
    double offsetcomp;
    MatrixXf out;
    MatrixXf inpssd;
    double* perms;
    int* indexesforperms;
    double** rsinitarray;
    double* relativeacc;
    double* tx;
    double* ty;
    double* fftframpl;

  //  int po_list[21] = {10, 11, 11, 11, 11, 12, 12, 11, 9, 10, 11, 11, 11, 11, 11, 12, 11, 9, 12, 10, 10};
  //  int po_list[21] = {10, 10, 10, 10, 10, 10, 10, 10, 9, 10, 10, 10, 10, 10, 10, 10, 10, 9, 10, 10, 10};
    int po_list[21] = {9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9};

    int parietal_occip_list[21][12] =
    {
       {9, 10, 11, 12, 13, 22, 23, 24, 25, 26, 0, 0},
       {9, 10, 11, 12, 13, 14, 22, 23, 24, 25, 26, 0},
       {9, 10, 11, 12, 13, 14, 22, 23, 24, 25, 26, 0},
       {9, 10, 11, 12, 13, 14, 23, 24, 25, 26, 27, 0},
       {8, 9, 10, 11, 12, 13, 20, 21, 22, 23, 24, 0},
       {9, 10, 11, 12, 13, 14, 21, 22, 23, 24, 25, 26},
       {7, 8, 9, 10, 11, 12, 21, 22, 23, 24, 25, 26},
       {8, 9, 10, 11, 12, 21, 22, 23, 24, 25, 26, 0},
       {9, 10, 11, 12, 13, 24, 25, 26, 27, 0, 0, 0},
       {8, 9, 10, 11, 12, 21, 22, 23, 24, 25, 0, 0},
       {9, 10, 11, 12, 13, 14, 22, 23, 24, 25, 26, 0},
       {8, 9, 10, 11, 12, 21, 22, 23, 24, 25, 26, 0},
       {7, 8, 9, 10, 11, 12, 21, 22, 23, 24, 25, 0},
       {7, 8, 9, 10, 11, 20, 21, 22, 23, 24, 25, 0},
       {8, 9, 10, 11, 12, 19, 20, 21, 22, 23, 24, 0},
       {9, 10, 11, 12, 13, 14, 20, 21, 22, 23, 24, 25},
       {8, 9, 10, 11, 12, 19, 20, 21, 22, 23, 24, 0},
       {8, 9, 10, 11, 19, 20, 21, 22, 23, 0, 0, 0},
       {9, 10, 11, 12, 13, 14, 22, 23, 24, 25, 26, 27},
       {8, 9, 10, 11, 12, 20, 21, 22, 23, 24, 0, 0},
       {8, 9, 10, 11, 12, 19, 20, 21, 22, 23, 0, 0}
    };

    double alphaarr[22][2] =
    {  {10.2, 10.8},
       {10.3, 10.5},
       {10.5, 10.3},
       {10.8, 10.9},
       {9.1, 9.1},
       {10.2, 10.1},
       {10.4, 9.3},
       {9.4, 9.6},
       {9.4, 9.8},
       {10, 11},
       {10.4, 10.3},
       {10.7, 10.6},
       {11, 0.0},
       {9.6, 9.4},
       {9.8, 10.2},
       {10, 10.6},
       {10.1, 11},
       {10.8, 11.1},
       {9.8, 10},
       {10.2, 10.6},
       {10.5, 10.8},
       {10.2, 10.10}
    };

    int inphasearr[19] = {1,4,3,3,2,3,2,3,1,3,4,4,4,2,3,2,3,1,4};
    int antiphasearr[19] = {3,1,2,4,4,4,4,4,3,4,3,3,1,3,4,4,2,2,1};
    int randphasearr[19] = {4,2,4,2,1,1,1,1,2,2,2,1,3,4,1,3,1,3,2};
    int noisearr[19] = {2,3,1,1,3,2,3,2,4,1,1,2,2,1,2,1,4,4,3};

    int blocksarr[19][4] =
    {  {1, 0, -1, 3}, // 1 0 1 3
       {-1, 3, 0, 1},
       {0, -1, 1, 3},
       {0, 3, 1, -1},
       {3, 1, 0, -1},
       {3, 0, 1, -1},
       {3, 1, 0, -1}, // 3 3 0 -1
       {3, 0, 1, -1},
       {1, 3, -1, 0},
       {0, 3, 1, -1},
       {0, 3, -1, 1},
       {3, 0, -1, 1},
       {-1, 0, 3, 1},
       {0, 1, -1, 3},
       {3, 0, 1, -1},
       {0, 1, 3, -1}, // 0 1 3 3
       {3, -1, 1, 0},  // 3 3 1 0
       {1, -1, 3, 0},
       {-1, 3, 0, 1}
      // {-1, 3, 0, 1},
      // {0, 3, 1, -1},
      // {0, -1, 1, 3}
    };

    int cspcomparr2[19][4] =
    {  {3, 2, 5, 5},
       {7, 8, 6, 3},
       {7, 7, 4, 5},
       {5, 8, 5, 3},
       {6, 5, 6, 3},
       {5, 6, 5, 7},
       {8, 4, 3, 4},
       {7, 6, 6, 6},
       {8, 8, 7, 6},
       {2, 8, 5, 10},
       {1, 6, 5, 3},
       {6, 7, 2, 8},
       {10, 3, 0, 0},
       {3, 2, 6, 2},
       {8, 7, 6, 10},
       {8, 7, 5, 7},
       {9, 9, 6, 10},
       {2, 8, 6, 5},
       {9, 6, 10, 3}
    };

    int channelsarr[19][4][9] =
    {
        {
        {11,12,13,21,22,23,24,25,26},
        {11,12,13,22,23,24,25,26,27},
        {11,12,13,14,21,22,23,24,25},
        {11,12,13,14,21,22,23,24,25}},
        {
        {10,11,12,13,22,23,24,25,26},
        {10,11,12,13,22,23,24,25,26},
        {10,11,12,13,22,23,24,25,26},
        {10,11,12,13,22,23,24,25,26}},
        {
        {10,11,12,13,22,23,24,25,26},
        {10,11,12,13,22,23,24,25,26},
        {10,11,12,13,22,23,24,25,26},
        {10,11,12,13,22,23,24,25,26}},
        {
        {10,11,12,13,22,23,24,25,26},
        {10,11,12,13,22,23,24,25,26},
        {10,11,12,13,22,23,24,25,26},
        {10,11,12,13,22,23,24,25,26}},
        {
        {10,11,12,13,22,23,24,25,26},
        {10,11,12,13,22,23,24,25,26},
        {10,11,12,13,23,24,25,26,27},
        {9,10,11,12,13,22,23,24,25}},
        {
        {9,10,11,12,22,23,24,25,26},
        {9,10,11,12,22,23,24,25,26},
        {9,10,11,12,21,22,23,24,25},
        {10,11,12,13,23,24,25,26,27}},
        {
        {10,11,12,13,22,23,24,25,26},
        {10,11,12,13,23,24,25,26,27},
        {10,11,12,13,24,25,26,27,28},
        {9,10,11,12,21,22,23,24,25}},
        {
        {10,11,12,13,22,23,24,25,26},
        {10,11,12,13,22,23,24,25,26},
        {10,11,12,13,21,22,23,24,25},
        {10,11,12,13,22,23,24,25,26}},
        {
        {10,11,12,13,22,23,24,25,26},
        {9,10,11,12,22,23,24,25,26},
        {10,11,12,13,22,23,24,25,26},
        {10,11,12,13,22,23,24,25,26}},
        {
        {10,11,12,13,22,23,24,25,26},
        {10,11,12,13,23,24,25,26,27},
        {9,10,11,12,21,22,23,24,25},
        {8,9,10,11,20,21,22,23,24}},
        {
        {10,11,12,13,22,23,24,25,26},
        {10,11,12,13,22,23,24,25,26},
        {10,11,12,13,21,22,23,24,25},
        {10,11,12,13,22,23,24,25,26}},
        {
        {10,11,12,13,22,23,24,25,26},
        {10,11,12,13,22,23,24,25,26},
        {10,11,12,13,14,22,23,24,25},
        {10,11,12,13,14,22,23,24,25}},
        {
        {10,11,12,13,22,23,24,25,26},
        {10,11,12,13,22,23,24,25,26},
        {0,0,0,0,0,0,0,0,0},
        {0,0,0,0,0,0,0,0,0}},
        {
        {10,11,12,13,21,22,23,24,25},
        {10,11,12,13,21,22,23,24,25},
        {9,10,11,12,13,14,22,23,24},
        {9,10,11,12,13,14,24,25,26}},
        {
        {7,8,9,10,20,21,22,23,24},
        {8,9,10,11,12,24,25,26,27},
        {10,11,12,13,21,22,23,24,25},
        {10,11,12,13,20,21,22,23,24}},
        {
        {10,11,12,13,20,21,22,23,24},
        {10,11,12,13,21,22,23,24,25},
        {10,11,12,13,14,15,23,24,25},
        {10,11,12,13,14,15,24,25,26}},
        {
        {10,11,12,13,20,21,22,23,24},
        {10,11,12,13,14,15,24,25,26},
        {9,10,11,12,13,22,23,24,25},
        {9,10,11,12,13,14,22,23,24}},
        {
        {10,11,12,13,21,22,23,24,25},
        {10,11,12,13,22,23,24,25,26},
        {9,10,11,12,20,21,22,23,24},
        {10,11,12,13,21,22,23,24,25}},
        {
        {10,11,12,13,22,23,24,25,26},
        {10,11,12,13,22,23,24,25,26},
        {10,11,12,13,21,22,23,24,25},
        {10,11,12,13,21,22,23,24,25}}
    };

    int eyescondition[22][10] =
    { {1, -1, 1, -1, 1, -1, 1, -1, 1, -1},
      {1, -1, 1, -1, 1, -1, 1, -1, 1, -1},
      {-1, 1, -1, 1, -1, -1, 1, 1, -1, 0},   //  {-1, 1, -1, 1, -1, -1, 1, 1, -1, 0},
      {1, -1, 1, -1, 1, -1, 1, -1, 1, -1},
      {-1, 1, -1, 1, -1, 1, -1, 1, -1, 0},
      {1, -1, 1, -1, 1, -1, 1, -1, 1, -1},
      {-1, 1, -1, 1, -1, 1, -1, 1, -1, 1},
      {-1, 1, -1, 1, -1, 1, -1, 1, -1, 1},
      {-1, 1, -1, 1, -1, 1, -1, 1, -1, 1},
      {1, -1, 1, -1, 1, -1, 1, -1, 1, -1},
      {1, -1, 1, -1, 1, -1, 1, -1, 1, -1},
      {1, -1, 1, -1, 1, -1, 1, -1, 1, -1},
      {-1, 1, -1, 1, -1, 1, -1, 1, -1, 1},
      {1, -1, 1, -1, 1, -1, 1, -1, 1, -1},
      {1, -1, 1, -1, 1, -1, 1, -1, 1, -1},
      {-1, 1, -1, 1, -1, 1, -1, 1, -1, 1},
      {-1, 1, -1, 1, -1, 1, -1, 1, -1, 1},
      {1, -1, 1, -1, 1, -1, 1, -1, 1, -1},
      {-1, 1, -1, 1, -1, 1, -1, 1, -1, 1},
      {1, -1, 1, -1, 1, -1, 1, -1, 1, -1},
      {1, -1, 1, -1, 1, -1, 1, -1, 1, -1},
      {-1, 1, -1, 1, -1, 1, -1, 1, -1, 1}
    };
   // int cspcomparr[22] = {3,3,3,3,6,3,3,3,3,7,3,2,2,5,1,2,2,4,3,1,3}; // bp 5-40 1 sec
  //  int cspcomparr[22] = {8,8,5,7,5,3,4,2,5,2,5,3,2,5,1,1,3,4,6,2,5}; // bp 7-14 1 sec
  //  int cspcomparr[22] = {3,2,5,5,6,2,4,3,4,5,3,4,3,4,4,2,4,2,4,1,5}; // bp 7-14 1/2 sec
  //  int cspcomparr[22] = {4,2,3,3,4,1,4,2,6,5,2,4,1,4,2,2,2,1,1,1,6}; // bp 7-14 1/2 sec used
  //  int cspcomparr[22] = {5,3,2,6,5,5,6,2,2,2,4,4,2,4,1,2,1,4,1,2,4}; // bp 7-14 800 msec

  //  int cspcomparr[22] = {4,6,3,3,7,1,4,2,6,5,2,6,1,4,2,2,2,1,1,1,5}; // bp 7-14 1/2 sec new
  //  int cspcomparr[22] = {4,7,7,5,9,8,5,4,6,6,8,4,4,7,1,2,4,4,1,2,8}; // bp 7-14 800 msec new

    int cspcomparr[22] = {4,4,3,3,8,1,4,5,8,5,2,6,1,4,2,2,4,7,1,1,10}; // bp 7-14 1/2 sec novel
  //  int cspcomparr[22] = {9,8,10,5,7,5,7,4,9,6,8,7,4,7,1,2,2,4,8,2,10}; // bp 7-14 800 msec novel

    int blockparts;
    int currentsubj;
    double** prestimplvs;
    double** poststimplvs;
    double** prestimsynchs;
    double** poststimsynchs;
    double** phaseshiftdiffs;
    double** postinstphasediffs;
    double** preinstphasediffs;
    double** rawdata;
    double** resblock; int blockstart; bool artifactremove;
    double* postppi;
    double* preppi;
    double* spectEOpre;
    double* spectEOpost;
    double* spectECpre;
    double* spectECpost;
    double* specttemp_pre;
    double* specttemp_post;
    double* meantrEOin;
    double* meantrEOanti;
    double* meantrECin;
    double* meantrECanti;
    int* indexes;
    int* relations; int totalstims;
    double* DenC;
    double* NumC;
    int samplenums,samplesize;
    int fftphase;
    double plvhr, plvzc, plvHRonAR, plvfft;
    double ECin, ECanti, EOin, EOanti;
    int cECin, cECanti, cEOin, cEOanti;
    freqband sigb,noibp,noibs;
    MatrixXf Xssd, wW; bool extractingssd,applyssd; int ssdcomp, maxssdcomp, ssdtime, chnum; double ssdt; // SSD
    MainWindow* mw;
    Settings* sw;      
    QElapsedTimer timer;
    QTimer* tim;
    QTimer* ssdtim;
    QByteArray stimarr;  
    hilbert* hnt;
    coords arrc;
    appconnect* appcn;
    TaskHandle taskHandle = 0;
    TaskHandle task0;
    QStringList strLst1, strLst2, strLst3;
    QStringListModel *strLstM1, *strLstM2, *strLstM3;

    explicit plotwindow(QWidget *parent = 0);    
    ~plotwindow();
    void plot(QCustomPlot *customPlot);
    bool eventFilter(QObject *target, QEvent *event);

    void clean();
    void setrandomrelations();
    void delay(int temp);
    double gethindex();
    void saverestofile();
    void parseresfromfile();    
    void testsyntdata();
    void sortblock(QString fname, int blocknum, int blocklength);
    double phasediff(double* stim1, double* stim2, int length);
    void minphasediff(int posstim, int length);
    void maxphasediff(int posstim, int length);
    void fillstimforstim(int posstim, int length);
    void fillstim(int posstim, int length);
    void fillrandstim(int posstim, int length);
    void clearstim(int posstim, int length, int method);
    void doplot();
    void strategy(int pos);
    void singlestim();
    void sequencestim();
    void clearstim();
    void cleareeg(int start, int end);
    void clearfiltered();
    void draw(int start);
    void fordemoshow();
    void bcidata(double d1);
    void makestimoneeg(int start);
    void savedatatofile(QString fname);
    void writerelationstofile(QString fname);
    void loadstimfromfile(QString fname);
    void meanpowerfromopregion2(QString fname, int trnum, int blocki, int currsubj);
    void generatesignals(int pos, int length, double fr1, double fr2, double fr3, double fr4, int phase, double noise);
    void refresh();
    int estimateoptlength(int n, int l1, int l2, int pos);
    double estimateoptprop(int n, double p1, double p2, int pos);
    void setaddmode(bool f);  
    void setfixedrelations();
    int maxabsstim(int pos);
    void stimulationonly(int pos);
    double offlinevaluation(int pos, int i, int method);
    void writedaq(int pos);
    void maxentropyf(double *input, int length, int degree, double **ar, double *per, double *pef, double *h, double *g);
    void autoregress(double *input, int length, int degree, double *coeffs);
    void runautoreg(int pos);
    double zcdifference(int pos);
    void savescreen();
    void printtoresultstring(QString str);
    void printtoresultbox(QString str);
    void printtoresultsync(QString str);

    void fillstimforregr(int posstim, int length);
    void hilbonar(int pos);
    double arphasediff(int posstim, int length);

    complex<double>* morletwavelet(float* inp, int poswavelet, int length);
    void fftphasepred(int posstim, int length);
    void fillstimforfft(int posstim, int length);

    void filloptstim(int posstim, int length, int phase, int shift);
    void getrawdata(int chn, double val);
    void delayt(int temp);

    void filterdata(int posstim, int length);
    void filterar(int posstim, int length);

    void zerocrossingstim(int posstim, int length);
    void fillstimforzc(int posstim, int length, int phaseshift);    

    void recurbutterf(int order, double sfr, double hpf, double* x, int length);
    void recurbutterfilter(int posstim, int length);
    void recurbutterfilterar(int posstim, int length);

    double* ComputeLP(int FilterOrder );
    double* ComputeHP(int FilterOrder );
    double* TrinomialMultiply(int FilterOrder, double *b, double *c );
    double* ComputeNumCoeffs(int FilterOrder);
    double* ComputeNumCoeffs(int FilterOrder, double Lcutoff, double Ucutoff, double *DenC);
    double* ComputeDenCoeffs(int FilterOrder, double Lcutoff, double Ucutoff);
    void filter(int ord, double *a, double *b, int np, double *x, double *y);
    void butterfiltcoefs(int lcut, int hcut, int order, int sampr);

    void filtfilt(vectord B, vectord A, const vectord &X, vectord &Y);
    void filter(vectord B, vectord A, const vectord &X, vectord &Y, vectord &Zi);
    vectord subvector_reverse(const vectord &vec, int idx_end, int idx_start);
    void append_vector(vectord &vec, const vectord &tail);
    void add_index_const(vectori &indices, int value, size_t numel);
    void add_index_range(vectori &indices, int beg, int end);
    void zerophaseinit(int lcut, int hcut, int order, int sampr);
    void zerophasefilt(int posstim, int length);
    void zerophasefiltar(int posstim, int length);
    void zerophasefiltzc(int posstim, int length);
    void zerophasefiltssd(int posstim, int length);
    void testfiltfilt();

    void getcoeffs(vectord& acc, vectord& bcc, int ord);    
    MatrixXf covarian(MatrixXf mat);
    MatrixXf filtfiltmatr(vectord b, vectord a, MatrixXf inp, int snums, int ssize);
    MatrixXf filtfiltmatrbytrials(vectord b, vectord a, MatrixXf inp, int trialnums, int ssize);
    MatrixXf readdata(int k, int length, QString fname);
    MatrixXf readdatafromeeg(int t);
    void savessdtofile(MatrixXf xssd, int dim, QString fname);
    void readxsddfromfile(QString fname);
    void ssd(QString fname, int trialnum);
    void getssd(int t, int length, QString fname, int blnum);
    void applycsp(QString fname, int blnum, int comp, QString condition, QString state);
    void applycsphalf(QString fname, int blnum, int comp, QString condition, QString state);
    MatrixXf getpart(MatrixXf inp, int p);
    MatrixXf cutpart(MatrixXf inp, int p, int trnum, int int_len);
    void fill_result_in_vector(VectorXf pre, VectorXf post, int p);
    void getrandompermutation();
    double getpvalue(int size, double value, bool side);

    void readstimplvs();
    void readstimsynchs();
    void savedataforplv(QString fname);
    void extractandsaveresult(QString fname, int blockstartpos);
    void loadresultblock(QString fname, int blockstartpos);
    void blockinline(QString fname, int blocklength);
    void extractstimdata(QString fname, int numbl, bool draw);
    void analyzestimdata(int predperiods, int postperiods, int numbl, bool draw);
    void analyzeallsubjs();
    void analyzeallsubjsmean();
    void sortblocks(QString fname, int blnum, bool draw);
    double getphaseoninterval(int start, int length);
    void phasepreservation(int blocknum, int nump);
    void ppiforallsubjs();
    void loadrestingstate(QString fname);
    void analyzerestingstate();
    void analyzeallrestingstate();
    void analyzepower(int numbl);
    void analyzerspower();
    void analyzerspowerspect();
    void getmeantrials(VectorXf datapre, VectorXf datapost);
    void fillresfile(QString fresname, int numbl);
    void mergeresults(QString prefix);
    void analyzegammapow(int numbl);
    void getpowerspectrum(int numbl);
    void getpowerspectrumrs();
    void getcsptopo(QString fname, int blnum);
    void fillmeantrials(int numbl);
    void addchannelvalues(channinfo* chnin, int chnnums, double stdv);
    int getchnind(QString str);
    void getmeantopo();
    double getstdvalue(double *arr, int sz);
    void sortresultblocks(QString prefix);
    void alignresultblocks(QString prefix);
    void alignallresultblocks(QString prefix);
    void combineresultblocks();
    void meandatafromopregion(QString fname, int blnum, int currsubj);
    void meandatafromopregionDT(QString fname, int blnum, int currsubj);
    void meanspectfromopregion(QString fname, int blnum, int currsubj);
    void meanplvfromopregionDT(QString fname, int blnum, int currsubj);
    double plvfromchannel(MatrixXf inp, int chnum, int blnum, int trialnum, bool pre);
    double alphafromchannel(MatrixXf inp, int chnum, int sh, int zp);
    void spectfromchannel(MatrixXf inp, int chnum, int sh, int zp, bool fl);
    double getalphapower(double* arr, int sh, int zp);
    VectorXf getmpart(VectorXf inp, int p);
    void detrend(double* y, double* x, int m);
    void flipsign();
    void sortreferencelist();
    void geteceopower();
    void getchanlabels();

    int findoptardegree();
    void findoptzcparams();
    void saveoptparams();
    void updatedata();
    void dophasepred();
    void parseoptparams();

    int extractblock(QString fname);
    void analyzeblockpower(int numtr,int timew);
    void analyzespectra(int numtr, int timew);
    void analyze2ndexp();
    void saveresulttofile(int blocki, int currentsubj, QString fresname);
    void savepartstofiles(int blocki, int currentsubj, QString suffix);
    void savespecttofiles(int blocki, int currentsbj, QString suffix);
    void savesingledatatofile(QString fname, double iaf);
    void read2ndexpdata(QString fname);
    void getcsptopo2(QString fname, int trnum);
    void applycsp2ndexp(QString fname, int trnum, int comp, QString state);
    void extractpozrsspectra();
    void analyzersspectra();
    void analyzersspectrapoz(int pos1, int pos2, int numtr);
    void savecspdatatofile(QString fname, double iaf);
    void spectfromchannel2(MatrixXf inp, int chnum, int sh, bool fl);
    void meanspectfromopregion2(QString fname, int trnum, int blocki, int currsubj);
    void parselrtcres(QString stim);
    void mergeresfiles(QString fres1, QString fres2, QString fres3);
    void fillblockparts(int trnums);
    void mergepowers();
    void mergepowerswithpoc();
    void mergeresfiles2(QString fres1, QString fres2);
    void analyze2ndplv(int numbl);

private slots:

    void timerUpdate();

    void ssdtimUpdate();

    void on_widget_destroyed();

    void on_pushButton_clicked();

    void on_checkBox_clicked();

    void on_checkBox_2_clicked();

    void on_radioButton_2_clicked();

    void on_radioButton_clicked();

    void on_checkBox_3_clicked();

    void on_pushButton_2_clicked();

    void on_pushButton_3_clicked();

    void on_pushButton_4_clicked();

    void on_checkBox_4_clicked();

    void slSettings();

    void on_doubleSpinBox_3_valueChanged(double arg1);

    void on_spinBox_7_valueChanged(int arg1);

    void on_spinBox_5_valueChanged(int arg1);

    void on_spinBox_valueChanged(int arg1);

    void on_spinBox_2_valueChanged(int arg1);

    void on_spinBox_3_valueChanged(int arg1);

    void on_doubleSpinBox_valueChanged(double arg1);

    void on_spinBox_6_valueChanged(int arg1);

    void on_spinBox_4_valueChanged(int arg1);

    void on_doubleSpinBox_2_valueChanged(double arg1);

    void on_spinBox_8_valueChanged(int arg1);

    void on_pushButton_5_clicked();

    void on_spinBox_5_editingFinished();

    void on_pushButton_6_clicked();

    void on_radioButton_3_clicked();

    void on_radioButton_4_clicked();

    void on_radioButton_5_clicked();

private:
    Ui::plot *ui;

};

#endif // PLOT_H
