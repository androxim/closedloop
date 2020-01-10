#include "settings.h"
#include "ui_settings.h"
#include "exludech.h"

Settings::Settings(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::Settings)
{
    ui->setupUi(this);   
    exch = new exludech();
    ui->checkBox_16->setDisabled(true);
}

Settings::~Settings()
{
    delete ui;
}

void Settings::on_pushButton_3_clicked() // "Apply" button
{
    if ((!ui->checkBox_13->isChecked()) && (!ui->checkBox_14->isChecked()) && (!ui->checkBox_8->isChecked()) && (!ui->checkBox_17->isChecked()))
    {
        QMessageBox msgBox;
        msgBox.setText("At least one method for stimulation should be enabled!");
        msgBox.exec();
    } else
    {
        pwd->showprediction=ui->checkBox_4->isChecked();
        pwd->estimpredstim=ui->checkBox_5->isChecked();
        pwd->filtereddata=ui->checkBox_6->isChecked();
        pwd->correctdelay=ui->checkBox_7->isChecked();        
        pwd->transdelay=ui->spinBox_10->value()/2; // for s.rate 500 Hz
        pwd->daqscalevalue=ui->spinBox_11->value();
        pwd->hronar=ui->checkBox_8->isChecked();
        pwd->hideardata=ui->checkBox_10->isChecked();
        pwd->flowpass=ui->spinBox_12->value();
        pwd->filterardata=ui->checkBox_11->isChecked();
        pwd->flength=ui->spinBox_13->value();
        pwd->carfilter=ui->checkBox_12->isChecked();
        pwd->offsetcomp=ui->doubleSpinBox_3->value();
        pwd->butterord=ui->spinBox_14->value();
        pwd->zerocr=ui->checkBox_13->isChecked();
        pwd->lcutoff=ui->spinBox_15->value();
        pwd->hcutoff=ui->spinBox_16->value();
        pwd->filtfir=ui->radioButton_7->isChecked();
        pwd->recurbutter=ui->radioButton_9->isChecked();
        pwd->zerobutter=ui->radioButton_8->isChecked();
        pwd->phasepr=ui->checkBox_14->isChecked();
        pwd->filterstims=ui->checkBox_15->isChecked();
        pwd->iti1=ui->spinBox_17->value()*2/3;
        pwd->iti2=pwd->iti1;
        pwd->applyssd=ui->checkBox_16->isChecked();
        pwd->sigb.lf=ui->spinBox_18->value();
        pwd->sigb.hf=ui->spinBox_19->value();
        pwd->noibp.lf=ui->spinBox_20->value();
        pwd->noibp.hf=ui->spinBox_21->value();
        pwd->ssdcomp=ui->spinBox_22->value();
        pwd->ssdtime=ui->spinBox_23->value();
        pwd->stimshift=ui->spinBox_24->value();
        pwd->fftpredict=ui->checkBox_17->isChecked();
       // pwd->autoregdegree=ui->spinBox_9->value();
        if ((ui->spinBox_9->value()>pwd->hnt->imlength) && (pwd->hronar))
        {
            QMessageBox msgBox;
            msgBox.setText("Autoregression order must be < imaging interval!");
            msgBox.exec();
        } else
            pwd->autoregdegree=ui->spinBox_9->value();
        if ((ui->checkBox_2->isChecked()) || (ui->checkBox_3->isChecked()))
            pwd->hnt->numst=ui->spinBox_7->value();
        pwd->hnt->osfr=esavg;
        pwd->hnt->imlength=eslength;
        pwd->hnt->imstprop=esprop;
        pwd->hnt->stlength=pwd->hnt->imlength*pwd->hnt->imstprop;
        if (ui->radioButton->isChecked())
            pwd->adaptampl=true;
        else
            pwd->adaptampl=false;
        if (ui->radioButton_3->isChecked())
            pwd->setaddmode(true);
        else
            pwd->setaddmode(false);
        pwd->refresh();   
        close();
    }
}

void Settings::on_pushButton_2_clicked()
{
    close();
}

void Settings::on_pushButton_clicked() // "Run" button
{  
    if (ui->checkBox->isChecked())
    {
        if (ui->spinBox_2->value()<ui->spinBox->value()+pwd->hnt->imlength)
        {
            QMessageBox msgBox;
            msgBox.setText("Incorrect interval for estimation!");
            msgBox.exec();
        } else
        {
            ui->label_6->setText("");
            pwd->hnt->averagefreq();
            esavg=pwd->hnt->avg1;
            ui->label_6->setText("Result: "+QString::number(esavg,'f',1));
        }
    }
    if (ui->checkBox_2->isChecked())
    {              
        if (ui->spinBox_4->value()<ui->spinBox_3->value()+10)
        {
            QMessageBox msgBox;
            msgBox.setText("Incorrect interval for estimation!");
            msgBox.exec();
        } else
        {
            ui->label_5->setText("");
            eslength=pwd->estimateoptlength(ui->spinBox_5->value(),ui->spinBox_3->value(),ui->spinBox_4->value(),ui->spinBox_6->value());
            ui->label_5->setText("Result: "+QString::number(eslength));
        }
    }
    if (ui->checkBox_3->isChecked())
    {
        ui->label_10->setText("");
        esprop=pwd->estimateoptprop(ui->spinBox_7->value(),ui->doubleSpinBox->value(),ui->doubleSpinBox_2->value(),ui->spinBox_8->value());
        ui->label_10->setText("Result: "+QString::number(esprop,'f',2));
    }
}

void Settings::init()
{
    int l=pwd->arrc.amp1.length();
    esavg=pwd->hnt->osfr;
    eslength=pwd->hnt->imlength;
    esprop=pwd->hnt->imstprop;
    ui->spinBox_10->setValue(pwd->transdelay*2);
    ui->spinBox_11->setValue(pwd->daqscalevalue);
    ui->spinBox_23->setValue(pwd->ssdtime);
    ui->spinBox_8->setEnabled(ui->checkBox_3->isChecked());
    ui->spinBox_7->setEnabled(ui->checkBox_3->isChecked());
    ui->doubleSpinBox->setEnabled(ui->checkBox_3->isChecked());
    ui->doubleSpinBox_2->setEnabled(ui->checkBox_3->isChecked());
    ui->doubleSpinBox_3->setValue(pwd->offsetcomp);
    ui->spinBox_3->setEnabled(ui->checkBox_2->isChecked());
    ui->spinBox_4->setEnabled(ui->checkBox_2->isChecked());
    ui->spinBox_5->setEnabled(ui->checkBox_2->isChecked());
    ui->spinBox_6->setEnabled(ui->checkBox_2->isChecked());
    ui->spinBox->setEnabled(ui->checkBox->isChecked());
    ui->spinBox_2->setEnabled(ui->checkBox->isChecked());  
    ui->spinBox_17->setValue(pwd->iti1*3/2);
    ui->radioButton_7->setChecked(pwd->filtfir);    
    ui->radioButton_9->setChecked(pwd->recurbutter);
    ui->radioButton_8->setChecked(pwd->zerobutter);
    ui->checkBox_4->setChecked(pwd->showprediction);
    ui->checkBox_5->setChecked(pwd->estimpredstim);
    ui->checkBox_6->setChecked(pwd->filtereddata);
    ui->checkBox_7->setChecked(pwd->correctdelay);
    ui->checkBox_8->setChecked(pwd->hronar);  
    ui->checkBox_10->setChecked(pwd->hideardata);
    ui->checkBox_11->setChecked(pwd->filterardata);
    ui->checkBox_12->setChecked(pwd->carfilter);
    ui->checkBox_13->setChecked(pwd->zerocr);
    ui->checkBox_14->setChecked(pwd->phasepr);
    ui->checkBox_15->setChecked(pwd->filterstims);
    ui->checkBox_16->setChecked(pwd->applyssd);
    ui->checkBox_17->setChecked(pwd->fftpredict);
    ui->spinBox->setMaximum(l-pwd->hnt->imlength);
    ui->spinBox_2->setMaximum(l);    
    ui->spinBox_5->setValue(pwd->hnt->numst);
    ui->spinBox_7->setValue(pwd->hnt->numst);
    ui->spinBox->setValue(pwd->hnt->posstim);
    ui->spinBox_2->setValue(pwd->hnt->posstim+pwd->hnt->imlength);
    ui->spinBox_6->setMaximum(l-pwd->hnt->imlength);
    ui->spinBox_8->setMaximum(l-pwd->hnt->imlength);
    ui->spinBox_3->setValue(pwd->hnt->imlength);
    ui->spinBox_6->setValue(pwd->hnt->posstim);
    ui->spinBox_8->setValue(pwd->hnt->posstim);
    ui->spinBox_9->setValue(pwd->autoregdegree);
    ui->spinBox_12->setValue(pwd->flowpass);      
    ui->spinBox_13->setValue(pwd->flength);
    ui->spinBox_14->setValue(pwd->butterord);
    ui->spinBox_15->setValue(pwd->lcutoff);
    ui->spinBox_16->setValue(pwd->hcutoff);
    ui->spinBox_18->setValue(pwd->sigb.lf);
    ui->spinBox_19->setValue(pwd->sigb.hf);
    ui->spinBox_20->setValue(pwd->noibp.lf);
    ui->spinBox_21->setValue(pwd->noibp.hf);
    ui->spinBox_22->setValue(pwd->ssdcomp);
    ui->spinBox_22->setMaximum(pwd->maxssdcomp);
    ui->spinBox_24->setValue(pwd->stimshift);
    ui->label_5->setText("");
    ui->label_6->setText("");
    ui->label_10->setText("");
    if (pwd->addmoderec)
        ui->radioButton_3->setChecked(true);
    else
        ui->radioButton_4->setChecked(true);
    if (pwd->adaptampl)
        ui->radioButton->setChecked(true);
    else
        ui->radioButton_2->setChecked(true);
}

void Settings::on_spinBox_valueChanged(int arg1)
{
    if (arg1>ui->spinBox_2->value())
        ui->spinBox_2->setValue(arg1+pwd->hnt->imlength);
    ui->spinBox_2->setMaximum(arg1+LMAX);
}

void Settings::on_spinBox_3_valueChanged(int arg1)
{
    if (arg1>ui->spinBox_4->value())
        ui->spinBox_4->setValue(arg1+10);
}

void Settings::on_checkBox_2_clicked()
{
    if (ui->checkBox_3->isChecked())
        ui->checkBox_3->setChecked(false);
    ui->spinBox_3->setEnabled(ui->checkBox_2->isChecked());
    ui->spinBox_4->setEnabled(ui->checkBox_2->isChecked());
    ui->spinBox_5->setEnabled(ui->checkBox_2->isChecked());
    ui->spinBox_6->setEnabled(ui->checkBox_2->isChecked());
}

void Settings::on_checkBox_3_clicked()
{    
    if (ui->checkBox_2->isChecked())
        ui->checkBox_2->setChecked(false);
    ui->spinBox_8->setEnabled(ui->checkBox_3->isChecked());
    ui->spinBox_7->setEnabled(ui->checkBox_3->isChecked());
    ui->doubleSpinBox->setEnabled(ui->checkBox_3->isChecked());
    ui->doubleSpinBox_2->setEnabled(ui->checkBox_3->isChecked());
}

void Settings::on_checkBox_clicked()
{
    ui->spinBox->setEnabled(ui->checkBox->isChecked());
    ui->spinBox_2->setEnabled(ui->checkBox->isChecked());
}

void Settings::on_toolButton_clicked()
{
    exch->totalchnumb=0;
    exch->plw=pwd;
    exch->show();
}

void Settings::on_checkBox_8_clicked()
{   
    ui->spinBox_9->setDisabled(!ui->checkBox_8->isChecked());
}

void Settings::on_pushButton_4_clicked()
{
    if (pwd->appcn->ready)
    {
        pwd->sigb.lf=ui->spinBox_18->value();
        pwd->sigb.hf=ui->spinBox_19->value();
        pwd->noibp.lf=ui->spinBox_20->value();
        pwd->noibp.hf=ui->spinBox_21->value();
        pwd->ssdcomp=ui->spinBox_22->value();

        siml=pwd->hnt->imlength;
        pwd->hnt->imlength=pwd->hnt->srfr*pwd->ssdtime;
      //  for (int i=0; i<pwd->chnums; i++)
      //      pwd->indexes[i]=0;
        pwd->extractingssd=true;
    }
}

void Settings::ssdready()
{
    ui->checkBox_16->setDisabled(false);
}


void Settings::on_spinBox_10_editingFinished()
{

}

void Settings::on_spinBox_9_editingFinished()
{

}

void Settings::on_checkBox_5_clicked()
{

}
