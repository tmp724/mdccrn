#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "QFileDialog"
#include "QMessageBox"
#include "QInputDialog"
#include "math.hpp"
#include "crn.hpp"
#include "monotonedependenciescalculator.hpp"

static bool filename_set = 0;
static bool output_dir_name_set = 0;
static bool model_name_set = 0;
QString qoutput_dir_name;
static std::string filename = "example_input.txt";
static std::string output_dir_name = "output";
static std::string model_name = "model";

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
}

MainWindow::~MainWindow()
{
    delete ui;
}

// input file button
void MainWindow::on_pushButton_clicked()
{
    QString qfilename = QFileDialog::getOpenFileName(this,"Please choose an input file", ".");
    QMessageBox::information(this,"",qfilename);
    filename = qfilename.toStdString();
    filename_set = 1;
}

// output dir button
void MainWindow::on_pushButton_2_clicked()
{
  qoutput_dir_name = QFileDialog::getExistingDirectory(this,"Please choose an output directory that all of the output files are going to be written into. It is recommended to create a dedicated directory for this.", ".");
  QMessageBox::information(this,"",qoutput_dir_name);
  output_dir_name = qoutput_dir_name.toStdString();
  output_dir_name_set = 1;
}

// run button
void MainWindow::on_pushButton_3_clicked()
{
  if(!filename_set){
    QMessageBox::information(this,"","Please choose an input file!");
  }else if(!output_dir_name_set){
    QMessageBox::information(this,"","Please choose an output directory!");
  }else if(!model_name_set){
    QMessageBox::information(this,"","Please choose a model name!");
  }else{
      QMessageBox::information(this,"","Calculation started! Please note that, depending on the size of the chosen CRN, this may take some time. "
                                         "Please look into the provided README file for more information on program speed.");
      CRN crn1(filename, model_name);
      MonotoneDependenciesCalculator mdc(crn1, 0);
      mdc.run();
      mdc.log(output_dir_name);
      QString tmp = "Done! Your output files are in " + qoutput_dir_name;
      QMessageBox::information(this,"", tmp);

  }
}

// exit button
void MainWindow::on_pushButton_4_clicked()
{
    MainWindow::close();
}

// model name button
void MainWindow::on_pushButton_5_clicked()
{
    QString qmodel_name = QInputDialog::getText(this,"","Please choose a model name!");
    model_name = qmodel_name.toStdString();
    model_name_set = 1;
}

void MainWindow::on_progressBar_valueChanged(int value)
{

}
