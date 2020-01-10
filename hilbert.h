#ifndef HILBERT_H
#define HILBERT_H

#include <QMainWindow>
#include <QWidget>
#include <QApplication>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <chrono>
#include <iostream>
#include <fstream>
#include <string>
#include <complex>
#include <valarray>
#include <QFile>
#include <QTextStream>
#include <QIODevice>
#include <QDebug>
#include <QElapsedTimer>
#include <QScrollArea>
#include <qcustomplot.h>

#define NMAX  1250000
#define LMAX  30000
#define LMFILT 1024	/* Finite Impulse Response (FIR) filter length - must be even! */

using namespace std;

const double PI = 3.141592653589793238460;
const double ToDeg = 57.2958;

typedef complex<double> Complex;
typedef valarray<Complex> CArray;

class MainWindow;

class hilbert
{
public:    

    int npt, lfilt, pos, posstim, srfr ,phase, initdelay, noise, stlength, imlength, stampl, sh, numst, phsdiffdeg, optphase, optphasehr, optphasear, optphasezc, optphasefft;
    double osfr, thfr, avg1, avg2, defosfr, phsdif, mintheta, maxtheta, mind, maxd, diff, extrphase, extinterval, arinterval;
    double totalprocpred, totalprocregr, totalproczc, totalprocfft, imstprop, totaldiffHR, totaldiffZC, totaldiffHRonAR, totaldiffFFT, degdevHR, degdevAR, degdevZC, degdevFFT, optphasediff;
    double averagedegdiffHR, averagedegdiffAR, averagedegdiffZC, averagedegdiffFFT, totalsynchr, totalsyncar, totalsynczc, totalsyncfft;
    bool twooscil, inphase;

    Complex plv;
    float *time, *x;
    float xh[LMAX + 1], ampl1[LMAX + 1], phase1[LMAX + 1], ampl2[LMAX + 1], phase2[LMAX + 1], omega1[LMAX + 1], omega2[LMAX + 1], intphs[LMAX+1];
    float hilb[LMFILT + 1], xt, xht;
    float* binborders;
    int * phasebins;
    int binnums;

    MainWindow* mww;
    string channel;
    QString filep;
    Complex t[LMAX];
    CArray data;

    hilbert(int srfr, int pstim, int lft, double osfr, double thfr, int phase, int ampl, int numst);
    hilbert(const hilbert &obj);
    void init();   
    void firstinit(QString filepath, int end);
    void fillphasebins(double phase);
    void readfromfile();

    void convol(float *source, float *target, float *filt, int npt, int lfilt);
    void firhilbert();

    void fft(CArray& x);
    void ifft(CArray& x);
    void matlabhilbert();

	void averagefreq();      

};

#endif // HILBERT_H
