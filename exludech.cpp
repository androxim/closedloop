#include "exludech.h"
#include "ui_exludech.h"
#include "plotwindow.h"

exludech::exludech(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::exludech)
{
    ui->setupUi(this);
    totalchnumb=0;
}

exludech::~exludech()
{
    delete ui;
}

void exludech::process()
{
    if (ui->checkBox->isChecked())
        plw->exlchannels[0]=0; else plw->exlchannels[0]=1;
    if (ui->checkBox_2->isChecked())
        plw->exlchannels[1]=0; else plw->exlchannels[1]=1;
    if (ui->checkBox_3->isChecked())
        plw->exlchannels[2]=0; else plw->exlchannels[2]=1;
    if (ui->checkBox_4->isChecked())
        plw->exlchannels[3]=0; else plw->exlchannels[3]=1;
    if (ui->checkBox_5->isChecked())
        plw->exlchannels[4]=0; else plw->exlchannels[4]=1;
    if (ui->checkBox_6->isChecked())
        plw->exlchannels[5]=0; else plw->exlchannels[5]=1;
    if (ui->checkBox_7->isChecked())
        plw->exlchannels[6]=0; else plw->exlchannels[6]=1;
    if (ui->checkBox_8->isChecked())
        plw->exlchannels[7]=0; else plw->exlchannels[7]=1;
    if (ui->checkBox_9->isChecked())
        plw->exlchannels[8]=0; else plw->exlchannels[8]=1;
    if (ui->checkBox_10->isChecked())
        plw->exlchannels[9]=0; else plw->exlchannels[9]=1;
    if (ui->checkBox_11->isChecked())
        plw->exlchannels[10]=0; else plw->exlchannels[10]=1;
    if (ui->checkBox_12->isChecked())
        plw->exlchannels[11]=0; else plw->exlchannels[11]=1;
    if (ui->checkBox_13->isChecked())
        plw->exlchannels[12]=0; else plw->exlchannels[12]=1;
    if (ui->checkBox_14->isChecked())
        plw->exlchannels[13]=0; else plw->exlchannels[13]=1;
    if (ui->checkBox_15->isChecked())
        plw->exlchannels[14]=0; else plw->exlchannels[14]=1;
    if (ui->checkBox_16->isChecked())
        plw->exlchannels[15]=0; else plw->exlchannels[15]=1;
    if (ui->checkBox_17->isChecked())
        plw->exlchannels[16]=0; else plw->exlchannels[16]=1;
    if (ui->checkBox_18->isChecked())
        plw->exlchannels[17]=0; else plw->exlchannels[17]=1;
    if (ui->checkBox_19->isChecked())
        plw->exlchannels[18]=0; else plw->exlchannels[18]=1;
    if (ui->checkBox_20->isChecked())
        plw->exlchannels[19]=0; else plw->exlchannels[19]=1;
    if (ui->checkBox_21->isChecked())
        plw->exlchannels[20]=0; else plw->exlchannels[20]=1;
    if (ui->checkBox_22->isChecked())
        plw->exlchannels[21]=0; else plw->exlchannels[21]=1;
    if (ui->checkBox_23->isChecked())
        plw->exlchannels[22]=0; else plw->exlchannels[22]=1;
    if (ui->checkBox_24->isChecked())
        plw->exlchannels[23]=0; else plw->exlchannels[23]=1;
    if (ui->checkBox_25->isChecked())
        plw->exlchannels[24]=0; else plw->exlchannels[24]=1;
    if (ui->checkBox_26->isChecked())
        plw->exlchannels[25]=0; else plw->exlchannels[25]=1;
    if (ui->checkBox_27->isChecked())
        plw->exlchannels[26]=0; else plw->exlchannels[26]=1;
    if (ui->checkBox_28->isChecked())
        plw->exlchannels[27]=0; else plw->exlchannels[27]=1;
    if (ui->checkBox_29->isChecked())
        plw->exlchannels[28]=0; else plw->exlchannels[28]=1;
    if (ui->checkBox_30->isChecked())
        plw->exlchannels[29]=0; else plw->exlchannels[29]=1;
    if (ui->checkBox_31->isChecked())
        plw->exlchannels[30]=0; else plw->exlchannels[30]=1;
    if (ui->checkBox_32->isChecked())
        plw->exlchannels[31]=0; else plw->exlchannels[31]=1;

    if (ui->checkBox_33->isChecked())
        plw->exlchannels[32]=0; else plw->exlchannels[32]=1;
    if (ui->checkBox_34->isChecked())
        plw->exlchannels[33]=0; else plw->exlchannels[33]=1;
    if (ui->checkBox_35->isChecked())
        plw->exlchannels[34]=0; else plw->exlchannels[34]=1;
    if (ui->checkBox_36->isChecked())
        plw->exlchannels[35]=0; else plw->exlchannels[35]=1;
    if (ui->checkBox_37->isChecked())
        plw->exlchannels[36]=0; else plw->exlchannels[36]=1;
    if (ui->checkBox_38->isChecked())
        plw->exlchannels[37]=0; else plw->exlchannels[37]=1;
    if (ui->checkBox_39->isChecked())
        plw->exlchannels[38]=0; else plw->exlchannels[38]=1;
    if (ui->checkBox_40->isChecked())
        plw->exlchannels[39]=0; else plw->exlchannels[39]=1;
    if (ui->checkBox_41->isChecked())
        plw->exlchannels[40]=0; else plw->exlchannels[40]=1;
    if (ui->checkBox_42->isChecked())
        plw->exlchannels[41]=0; else plw->exlchannels[41]=1;
    if (ui->checkBox_43->isChecked())
        plw->exlchannels[42]=0; else plw->exlchannels[42]=1;
    if (ui->checkBox_44->isChecked())
        plw->exlchannels[43]=0; else plw->exlchannels[43]=1;
    if (ui->checkBox_45->isChecked())
        plw->exlchannels[44]=0; else plw->exlchannels[44]=1;
    if (ui->checkBox_46->isChecked())
        plw->exlchannels[45]=0; else plw->exlchannels[45]=1;
    if (ui->checkBox_47->isChecked())
        plw->exlchannels[46]=0; else plw->exlchannels[46]=1;
    if (ui->checkBox_48->isChecked())
        plw->exlchannels[47]=0; else plw->exlchannels[47]=1;
    if (ui->checkBox_49->isChecked())
        plw->exlchannels[48]=0; else plw->exlchannels[48]=1;
    if (ui->checkBox_50->isChecked())
        plw->exlchannels[49]=0; else plw->exlchannels[49]=1;
    if (ui->checkBox_51->isChecked())
        plw->exlchannels[50]=0; else plw->exlchannels[50]=1;
    if (ui->checkBox_52->isChecked())
        plw->exlchannels[51]=0; else plw->exlchannels[51]=1;
    if (ui->checkBox_53->isChecked())
        plw->exlchannels[52]=0; else plw->exlchannels[52]=1;
    if (ui->checkBox_54->isChecked())
        plw->exlchannels[53]=0; else plw->exlchannels[53]=1;
    if (ui->checkBox_55->isChecked())
        plw->exlchannels[54]=0; else plw->exlchannels[54]=1;
    if (ui->checkBox_56->isChecked())
        plw->exlchannels[55]=0; else plw->exlchannels[55]=1;
    if (ui->checkBox_57->isChecked())
        plw->exlchannels[56]=0; else plw->exlchannels[56]=1;
    if (ui->checkBox_58->isChecked())
        plw->exlchannels[57]=0; else plw->exlchannels[57]=1;
    if (ui->checkBox_59->isChecked())
        plw->exlchannels[58]=0; else plw->exlchannels[58]=1;
    if (ui->checkBox_60->isChecked())
        plw->exlchannels[59]=0; else plw->exlchannels[59]=1;
    if (ui->checkBox_61->isChecked())
        plw->exlchannels[60]=0; else plw->exlchannels[60]=1;
    if (ui->checkBox_62->isChecked())
        plw->exlchannels[61]=0; else plw->exlchannels[61]=1;
    if (ui->checkBox_63->isChecked())
        plw->exlchannels[62]=0; else plw->exlchannels[62]=1;
    if (ui->checkBox_64->isChecked())
        plw->exlchannels[63]=0; else plw->exlchannels[63]=1;
}

void exludech::on_buttonBox_accepted()
{
    process();
    for (int i=0; i<plw->chnums; i++)
        if ((plw->exlchannels[i]==0) && (i!=plw->sourcech))
            totalchnumb++;
    plw->exch=plw->chnums-totalchnumb;
  //  qDebug()<<plw->exch;
}

void exludech::setcheckboxs(bool fl)
{
    ui->checkBox->setChecked(fl);
    ui->checkBox_2->setChecked(fl);
    ui->checkBox_3->setChecked(fl);
    ui->checkBox_4->setChecked(fl);
    ui->checkBox_5->setChecked(fl);
    ui->checkBox_6->setChecked(fl);
    ui->checkBox_7->setChecked(fl);
    ui->checkBox_8->setChecked(fl);
    ui->checkBox_9->setChecked(fl);
    ui->checkBox_10->setChecked(fl);
    ui->checkBox_11->setChecked(fl);
    ui->checkBox_12->setChecked(fl);
    ui->checkBox_13->setChecked(fl);
    ui->checkBox_14->setChecked(fl);
    ui->checkBox_15->setChecked(fl);
    ui->checkBox_16->setChecked(fl);
    ui->checkBox_17->setChecked(fl);
    ui->checkBox_18->setChecked(fl);
    ui->checkBox_19->setChecked(fl);
    ui->checkBox_20->setChecked(fl);
    ui->checkBox_21->setChecked(fl);
    ui->checkBox_22->setChecked(fl);
    ui->checkBox_23->setChecked(fl);
    ui->checkBox_24->setChecked(fl);
    ui->checkBox_25->setChecked(fl);
    ui->checkBox_26->setChecked(fl);
    ui->checkBox_27->setChecked(fl);
    ui->checkBox_28->setChecked(fl);
    ui->checkBox_29->setChecked(fl);
    ui->checkBox_30->setChecked(fl);
    ui->checkBox_31->setChecked(fl);
    ui->checkBox_32->setChecked(fl);
    ui->checkBox_33->setChecked(fl);
    ui->checkBox_34->setChecked(fl);
    ui->checkBox_35->setChecked(fl);
    ui->checkBox_36->setChecked(fl);
    ui->checkBox_37->setChecked(fl);
    ui->checkBox_38->setChecked(fl);
    ui->checkBox_39->setChecked(fl);
    ui->checkBox_40->setChecked(fl);
    ui->checkBox_41->setChecked(fl);
    ui->checkBox_42->setChecked(fl);
    ui->checkBox_43->setChecked(fl);
    ui->checkBox_44->setChecked(fl);
    ui->checkBox_45->setChecked(fl);
    ui->checkBox_46->setChecked(fl);
    ui->checkBox_47->setChecked(fl);
    ui->checkBox_48->setChecked(fl);
    ui->checkBox_49->setChecked(fl);
    ui->checkBox_50->setChecked(fl);
    ui->checkBox_51->setChecked(fl);
    ui->checkBox_52->setChecked(fl);
    ui->checkBox_53->setChecked(fl);
    ui->checkBox_54->setChecked(fl);
    ui->checkBox_55->setChecked(fl);
    ui->checkBox_56->setChecked(fl);
    ui->checkBox_57->setChecked(fl);
    ui->checkBox_58->setChecked(fl);
    ui->checkBox_59->setChecked(fl);
    ui->checkBox_60->setChecked(fl);
    ui->checkBox_61->setChecked(fl);
    ui->checkBox_62->setChecked(fl);
    ui->checkBox_63->setChecked(fl);
    ui->checkBox_64->setChecked(fl);
}

void exludech::on_pushButton_clicked()
{
     setcheckboxs(true);
}

void exludech::on_pushButton_2_clicked()
{
    setcheckboxs(false);
}
