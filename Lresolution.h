#ifndef _LRESOLUTION_H_
#define _LRESOLUTION_H_

#include <iostream>
#include <fstream>
#include <cmath>

#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TRandom.h"
#include "TMath.h"

class Lresolution{
 public:
  Lresolution(const char configuration[]);/* Initialization. The configuration can be PVDIS, SIDIS_He3, JPsi */
  int GetStatus();/* Print current configuration */
  int CheckConfiguration();/* Check the validity of current configuration. Return 0:normal, 1:Invalid pid, 2:out of range*/
  int CheckKinematic(const double * kin, const int pid);/* Check Kinematic in domain. Return 0:in,  1:out */
  int GetDomain(double * dom, const int pid);
  int GetDomain(double * dom, const char pname[]);/* Get kinematic domain */
  int PrintDomain();/* Output kinematic domain of particles on the screen */
  int GetResolution(const double * kin, double * res, const int pid);
  int GetResolution(const double * kin, double * res, const char pname[]);
/* Get resolution to res: delta_p/p(no unit), theta(deg), phi(deg), vertexz(cm)
  multiply by res_factor=1.5 to for safety
  Kinematic variables kin: p(GeV), theta(deg), phi(deg)
  Particle: e+(1), e-(-1), mu+(2), mu-(-2), pi+(3), pi-(-3), K+(4), K-(-4), p(5). Default with e-(-1)
  Return O:normal, 1:error, 2: out of range */
  int PrintResolution(const double * kin, const int pid);
  int PrintResolution(const double * kin, const char pname[]);/* Output the resolution on the screen */
 protected:
  Int_t _status;// configuration status 1:PVDIS, 2:SIDIS, 3:JPsi
  TFile * _file_electron;// resolution files
  TFile * _file_positron;
  TFile * _file_muonp;
  TFile * _file_muonm;
  TFile * _file_pionp;
  TFile * _file_pionm;
  TFile * _file_kaonp;
  TFile * _file_kaonm;
  TFile * _file_proton;
  Double_t _domain_electron[4];// domain of definition: pl, pu, thetal, thetau
  Double_t _domain_positron[4];
  Double_t _domain_muonp[4];
  Double_t _domain_muonm[4];
  Double_t _domain_pionp[4];
  Double_t _domain_pionm[4];
  Double_t _domain_kaonp[4];
  Double_t _domain_kaonm[4];
  Double_t _domain_proton[4];
  TH2D * electron_p;// resolution 2D histograms
  TH2D * electron_theta;
  TH2D * electron_phi;
  TH2D * electron_z;
  TH2D * positron_p;
  TH2D * positron_theta;
  TH2D * positron_phi;
  TH2D * positron_z;
  TH2D * muonp_p;
  TH2D * muonp_theta;
  TH2D * muonp_phi;
  TH2D * muonp_z;
  TH2D * muonm_p;
  TH2D * muonm_theta;
  TH2D * muonm_phi;
  TH2D * muonm_z;
  TH2D * pionp_p;
  TH2D * pionp_theta;
  TH2D * pionp_phi;
  TH2D * pionp_z;
  TH2D * pionm_p;
  TH2D * pionm_theta;
  TH2D * pionm_phi;
  TH2D * pionm_z;
  TH2D * kaonp_p;
  TH2D * kaonp_theta;
  TH2D * kaonp_phi;
  TH2D * kaonp_z;
  TH2D * kaonm_p;
  TH2D * kaonm_theta;
  TH2D * kaonm_phi;
  TH2D * kaonm_z;
  TH2D * proton_p;
  TH2D * proton_theta;
  TH2D * proton_phi;
  TH2D * proton_z;
};

Lresolution::Lresolution(const char configuration[]){
  _status = 0;
  bool control = true;
  //file names
  TString felectron, fpositron, fmuonp, fmuonm, fpionp, fpionm, fkaonp, fkaonm, fproton;
  //PVDIS configuration
  if (strcmp(configuration, "PVDIS") == 0){
    _status = 1;
    felectron = "res_file/PVDIS_2016_01_21/PVDIS_electron_resolution_2d.root";
    fpositron = "res_file/PVDIS_2016_01_21/PVDIS_electron_resolution_2d.root";
    fmuonp = "res_file/PVDIS_2016_01_21/PVDIS_electron_resolution_2d.root";
    fmuonm = "res_file/PVDIS_2016_01_21/PVDIS_electron_resolution_2d.root";
    fpionp = "res_file/PVDIS_2016_01_21/PVDIS_electron_resolution_2d.root";
    fpionm = "res_file/PVDIS_2016_01_21/PVDIS_electron_resolution_2d.root";
    fkaonp = "res_file/PVDIS_2016_01_21/PVDIS_electron_resolution_2d.root";
    fkaonm = "res_file/PVDIS_2016_01_21/PVDIS_electron_resolution_2d.root";
    fproton = "res_file/PVDIS_2016_01_21/PVDIS_electron_resolution_2d.root";
    std::cout << "SoLID PVDIS configuration" << std::endl;
    std::cout << "Notice: missing files for e+, mu+, mu-, pi+, pi-, K+, K-, p" << std::endl;
    std::cout << "        using e- file instead!" << std::endl;
  }
  //SIDIS configuration (Helium-3)
  else if (strcmp(configuration, "SIDIS_He3") == 0){
    _status = 2;
    felectron = "res_file/SIDIS_He3_2016_01_21/SIDIS_electron_resolution_2d.root";
    fpositron = "res_file/SIDIS_He3_2016_01_21/SIDIS_electron_resolution_2d.root";
    fmuonp = "res_file/SIDIS_He3_2016_01_21/SIDIS_pim_resolution_2d.root";
    fmuonm = "res_file/SIDIS_He3_2016_01_21/SIDIS_pim_resolution_2d.root";
    fpionp = "res_file/SIDIS_He3_2016_01_21/SIDIS_pim_resolution_2d.root";
    fpionm = "res_file/SIDIS_He3_2016_01_21/SIDIS_pim_resolution_2d.root";
    fkaonp = "res_file/SIDIS_He3_2016_01_21/SIDIS_pim_resolution_2d.root";
    fkaonm = "res_file/SIDIS_He3_2016_01_21/SIDIS_pim_resolution_2d.root";
    fproton = "res_file/SIDIS_He3_2016_01_21/SIDIS_pim_resolution_2d.root";
    std::cout << "SoLID SIDIS_He3 configuration" << std::endl;
    std::cout << "Notice: missing files for e+, mu+, mu-, pi+, K+, K-, p" << std::endl;
    std::cout << "        using e-, pi- files instead!" << std::endl;
  }
  else if (strcmp(configuration, "JPsi") == 0){
    _status = 3;
    felectron = "res_file/JPsi_2016_01_21/JPsi_electron_resolution_2d.root";
    fpositron = "res_file/JPsi_2016_01_21/JPsi_electron_resolution_2d.root";
    fmuonp = "res_file/JPsi_2016_01_21/JPsi_proton_resolution_2d.root";
    fmuonm = "res_file/JPsi_2016_01_21/JPsi_proton_resolution_2d.root";
    fpionp = "res_file/JPsi_2016_01_21/JPsi_proton_resolution_2d.root";
    fpionm = "res_file/JPsi_2016_01_21/JPsi_proton_resolution_2d.root";
    fkaonp = "res_file/JPsi_2016_01_21/JPsi_proton_resolution_2d.root";
    fkaonm = "res_file/JPsi_2016_01_21/JPsi_proton_resolution_2d.root";
    fproton = "res_file/JPsi_2016_01_21/JPsi_proton_resolution_2d.root";
    std::cout << "SoLID JPsi configuration" << std::endl;
    std::cout << "Notice: missing files for e+, mu+, mu-, pi+, pi-, K+, K-" << std::endl;
    std::cout << "        using e-, p files instead!" << std::endl;
  }
  else {
    std::cout << "Error in Lresolution::Lresolution!" << std::endl;
    std::cout << "  No configuration is matched! Choices: PVDIS, SIDIS_He3, JPsi." << std::endl;
    control = false;
  }
  if (control){
    _file_electron = new TFile(felectron, "r");
    _file_positron = new TFile(fpositron, "r");
    _file_muonp = new TFile(fmuonp, "r");
    _file_muonm = new TFile(fmuonm, "r");
    _file_pionp = new TFile(fpionp, "r");
    _file_pionm = new TFile(fpionm, "r");
    _file_kaonp = new TFile(fkaonp, "r");
    _file_kaonm = new TFile(fkaonm, "r");
    _file_proton = new TFile(fproton, "r");
    electron_p = (TH2D *) _file_electron->GetObjectChecked("p_resolution","TH2D");
    electron_theta = (TH2D *) _file_electron->GetObjectChecked("theta_resolution","TH2D");
    electron_phi = (TH2D *) _file_electron->GetObjectChecked("phi_resolution","TH2D");
    electron_z = (TH2D *) _file_electron->GetObjectChecked("vertexz_resolution","TH2D");
    positron_p = (TH2D *) _file_positron->GetObjectChecked("p_resolution","TH2D");
    positron_theta = (TH2D *) _file_positron->GetObjectChecked("theta_resolution","TH2D");
    positron_phi = (TH2D *) _file_positron->GetObjectChecked("phi_resolution","TH2D");
    positron_z = (TH2D *) _file_positron->GetObjectChecked("vertexz_resolution","TH2D");
    muonp_p = (TH2D *) _file_muonp->GetObjectChecked("p_resolution","TH2D");
    muonp_theta = (TH2D *) _file_muonp->GetObjectChecked("theta_resolution","TH2D");
    muonp_phi = (TH2D *) _file_muonp->GetObjectChecked("phi_resolution","TH2D");
    muonp_z = (TH2D *) _file_muonp->GetObjectChecked("vertexz_resolution","TH2D");
    muonm_p = (TH2D *) _file_muonm->GetObjectChecked("p_resolution","TH2D");
    muonm_theta = (TH2D *) _file_muonm->GetObjectChecked("theta_resolution","TH2D");
    muonm_phi = (TH2D *) _file_muonm->GetObjectChecked("phi_resolution","TH2D");
    muonm_z = (TH2D *) _file_muonm->GetObjectChecked("vertexz_resolution","TH2D");
    pionp_p = (TH2D *) _file_pionp->GetObjectChecked("p_resolution","TH2D");
    pionp_theta = (TH2D *) _file_pionp->GetObjectChecked("theta_resolution","TH2D");
    pionp_phi = (TH2D *) _file_pionp->GetObjectChecked("phi_resolution","TH2D");
    pionp_z = (TH2D *) _file_pionp->GetObjectChecked("vertexz_resolution","TH2D");
    pionm_p = (TH2D *) _file_pionm->GetObjectChecked("p_resolution","TH2D");
    pionm_theta = (TH2D *) _file_pionm->GetObjectChecked("theta_resolution","TH2D");
    pionm_phi = (TH2D *) _file_pionm->GetObjectChecked("phi_resolution","TH2D");
    pionm_z = (TH2D *) _file_pionm->GetObjectChecked("vertexz_resolution","TH2D");
    kaonp_p = (TH2D *) _file_kaonp->GetObjectChecked("p_resolution","TH2D");
    kaonp_theta = (TH2D *) _file_kaonp->GetObjectChecked("theta_resolution","TH2D");
    kaonp_phi = (TH2D *) _file_kaonp->GetObjectChecked("phi_resolution","TH2D");
    kaonp_z = (TH2D *) _file_kaonp->GetObjectChecked("vertexz_resolution","TH2D");
    kaonm_p = (TH2D *) _file_kaonm->GetObjectChecked("p_resolution","TH2D");
    kaonm_theta = (TH2D *) _file_kaonm->GetObjectChecked("theta_resolution","TH2D");
    kaonm_phi = (TH2D *) _file_kaonm->GetObjectChecked("phi_resolution","TH2D");
    kaonm_z = (TH2D *) _file_kaonm->GetObjectChecked("vertexz_resolution","TH2D");
    proton_p = (TH2D *) _file_proton->GetObjectChecked("p_resolution","TH2D");
    proton_theta = (TH2D *) _file_proton->GetObjectChecked("theta_resolution","TH2D");
    proton_phi = (TH2D *) _file_proton->GetObjectChecked("phi_resolution","TH2D");
    proton_z = (TH2D *) _file_proton->GetObjectChecked("vertexz_resolution","TH2D");
    _domain_electron[0] = electron_p->GetXaxis()->GetXmin();
    _domain_electron[1] = electron_p->GetXaxis()->GetXmax();
    _domain_electron[2] = electron_p->GetYaxis()->GetXmin();
    _domain_electron[3] = electron_p->GetYaxis()->GetXmax();
    _domain_positron[0] = positron_p->GetXaxis()->GetXmin();
    _domain_positron[1] = positron_p->GetXaxis()->GetXmax();
    _domain_positron[2] = positron_p->GetYaxis()->GetXmin();
    _domain_positron[3] = positron_p->GetYaxis()->GetXmax();
    _domain_muonp[0] = muonp_p->GetXaxis()->GetXmin();
    _domain_muonp[1] = muonp_p->GetXaxis()->GetXmax();
    _domain_muonp[2] = muonp_p->GetYaxis()->GetXmin();
    _domain_muonp[3] = muonp_p->GetYaxis()->GetXmax();
    _domain_muonm[0] = muonm_p->GetXaxis()->GetXmin();
    _domain_muonm[1] = muonm_p->GetXaxis()->GetXmax();
    _domain_muonm[2] = muonm_p->GetYaxis()->GetXmin();
    _domain_muonm[3] = muonm_p->GetYaxis()->GetXmax();
    _domain_pionp[0] = pionp_p->GetXaxis()->GetXmin();
    _domain_pionp[1] = pionp_p->GetXaxis()->GetXmax();
    _domain_pionp[2] = pionp_p->GetYaxis()->GetXmin();
    _domain_pionp[3] = pionp_p->GetYaxis()->GetXmax();
    _domain_pionm[0] = pionm_p->GetXaxis()->GetXmin();
    _domain_pionm[1] = pionm_p->GetXaxis()->GetXmax();
    _domain_pionm[2] = pionm_p->GetYaxis()->GetXmin();
    _domain_pionm[3] = pionm_p->GetYaxis()->GetXmax();
    _domain_kaonp[0] = kaonp_p->GetXaxis()->GetXmin();
    _domain_kaonp[1] = kaonp_p->GetXaxis()->GetXmax();
    _domain_kaonp[2] = kaonp_p->GetYaxis()->GetXmin();
    _domain_kaonp[3] = kaonp_p->GetYaxis()->GetXmax();
    _domain_kaonm[0] = kaonm_p->GetXaxis()->GetXmin();
    _domain_kaonm[1] = kaonm_p->GetXaxis()->GetXmax();
    _domain_kaonm[2] = kaonm_p->GetYaxis()->GetXmin();
    _domain_kaonm[3] = kaonm_p->GetYaxis()->GetXmax();
    _domain_proton[0] = proton_p->GetXaxis()->GetXmin();
    _domain_proton[1] = proton_p->GetXaxis()->GetXmax();
    _domain_proton[2] = proton_p->GetYaxis()->GetXmin();
    _domain_proton[3] = proton_p->GetYaxis()->GetXmax();
  }
}

int Lresolution::GetStatus(){
  std::cout << "Current status: ";
  if (_status == 1){
    std::cout << "PVDIS" << std::endl;
  }
  else if (_status == 2){
    std::cout << "SIDIS" << std::endl;
  }
  else if (_status == 3){
    std::cout << "JPsi" << std::endl;
  }
  else {
    std::cout << "None!" << std::endl;
  }
  return _status;
}

int Lresolution::CheckConfiguration(){
  if ( _status != 1 && _status != 2 && _status != 3){
    std::cout << "Invalid configuration!" << std::endl;
    return 1;
  }
  else 
    return 0;
}

int Lresolution::CheckKinematic(const double * kin, const int pid){
  if (CheckConfiguration()) {return 1;}
  double * domain;
  if (pid == 1) domain = _domain_positron;
  else if (pid == -1) domain = _domain_electron;
  else if (pid == 2) domain = _domain_muonp;
  else if (pid == -2) domain = _domain_muonm;
  else if (pid == 3) domain = _domain_pionp;
  else if (pid == -3) domain = _domain_pionm;
  else if (pid == 4) domain = _domain_kaonp;
  else if (pid == -4) domain = _domain_kaonm;
  else if (pid == 5) domain = _domain_proton;
  else {
    std::cout << "Error in Lresolution::GetDomain!" << std::endl;
    std::cout << "  Invalid pid." << std::endl;
    return 1;
  }
  if (kin[0] >= domain[0] && kin[0] < domain[1] && kin[1] >= domain[2] && kin[1] < domain[3]){
    return 0;
  }
  else {
    return 2;
  }
}

int Lresolution::GetDomain(double * dom, const int pid){
  if (CheckConfiguration()) {return 1;}
  //dom: pl, pu, thetal, thetau in unit of GeV and Deg
  double * domain;
  if (pid == 1) domain = _domain_positron;
  else if (pid == -1) domain = _domain_electron;
  else if (pid == 2) domain = _domain_muonp;
  else if (pid == -2) domain = _domain_muonm;
  else if (pid == 3) domain = _domain_pionp;
  else if (pid == -3) domain = _domain_pionm;
  else if (pid == 4) domain = _domain_kaonp;
  else if (pid == -4) domain = _domain_kaonm;
  else if (pid == 5) domain = _domain_proton;
  else {
    std::cout << "Error in Lresolution::GetDomain!" << std:: endl;
    std::cout << "  Invalid pid." << std::endl;
    return 1;
  }
  for (int i = 0; i < 4; i++){
    dom[i] = domain[i];
  }
  return 0;
}

int Lresolution::GetDomain(double * dom, const char pname[]){
  if (CheckConfiguration()) {return 1;}
  int pid;
  if (strcmp(pname, "e+") == 0) pid = 1;
  else if (strcmp(pname, "e-") == 0) pid = -1;
  else if (strcmp(pname, "mu+") == 0) pid = 2;
  else if (strcmp(pname, "mu-") == 0) pid = -2;
  else if (strcmp(pname, "pi+") == 0) pid = 3;
  else if (strcmp(pname, "pi-") == 0) pid = -3;
  else if (strcmp(pname, "K+") == 0) pid = 4;
  else if (strcmp(pname, "K-") == 0) pid = -4;
  else if (strcmp(pname, "p") == 0) pid = 5;
  else pid = -999;
  return GetDomain(dom, pid);
}

int Lresolution::PrintDomain(){
  if (CheckConfiguration()) {return 1;}
  printf("particle:  plow/GeV  pup/GeV  thetalow/Deg  thetaup/Deg\n");
  printf("e-      :  %.1f     %.1f       %.1f     %.1f\n", _domain_electron[0], _domain_electron[1], _domain_electron[2], _domain_electron[3]);
  printf("e+      :  %.1f     %.1f       %.1f     %.1f\n", _domain_positron[0], _domain_positron[1], _domain_positron[2], _domain_positron[3]);
  printf("mu+     :  %.1f     %.1f       %.1f     %.1f\n", _domain_muonp[0], _domain_muonp[1], _domain_muonp[2], _domain_muonp[3]);
  printf("mu-     :  %.1f     %.1f       %.1f     %.1f\n", _domain_muonm[0], _domain_muonm[1], _domain_muonm[2], _domain_muonm[3]);
  printf("pi+     :  %.1f     %.1f       %.1f     %.1f\n", _domain_pionp[0], _domain_pionp[1], _domain_pionp[2], _domain_pionp[3]);
  printf("pi-     :  %.1f     %.1f       %.1f     %.1f\n", _domain_pionm[0], _domain_pionm[1], _domain_pionm[2], _domain_pionm[3]);
  printf("K+      :  %.1f     %.1f       %.1f     %.1f\n", _domain_kaonp[0], _domain_kaonp[1], _domain_kaonp[2], _domain_kaonp[3]);
  printf("K-      :  %.1f     %.1f       %.1f     %.1f\n", _domain_kaonm[0], _domain_kaonm[1], _domain_kaonm[2], _domain_kaonm[3]);
  printf("p       :  %.1f     %.1f       %.1f     %.1f\n", _domain_proton[0], _domain_proton[1], _domain_proton[2], _domain_proton[3]);
  return 0;
} 

int Lresolution::GetResolution(const double * kin, double * res, const int pid){
  if (CheckConfiguration()) {return 1;}
  //kin: p/GeV, theta/Deg, phi/Deg
  //res: dp/GeV, dtheta/Deg, dphi/Deg, dz/cm
  /* pid: e+: 1      e-: -1
	 mu+: 2     mu-: -2
	 pi+: 3     pi-: -3
	  K+: 4      K-: -4
	   p: 5
  */
  TH2D * Rp;
  TH2D * Rtheta;
  TH2D * Rphi;
  TH2D * Rz;
  if (pid == 1){
    Rp = positron_p;
    Rtheta = positron_theta;
    Rphi = positron_phi;
    Rz = positron_z;
  }
  else if (pid == -1){
    Rp = electron_p;
    Rtheta = electron_theta;
    Rphi = electron_phi;
    Rz = electron_z;
  }
  else if (pid == 2){
    Rp = muonp_p;
    Rtheta = muonp_theta;
    Rphi = muonp_phi;
    Rz = muonp_z;
  }
  else if (pid == -2){
    Rp = muonm_p;
    Rtheta = muonm_theta;
    Rphi = muonm_phi;
    Rz = muonm_z;
  }
  else if (pid == 3){
    Rp = pionp_p;
    Rtheta = pionp_theta;
    Rphi = pionp_phi;
    Rz = pionp_z;
  }
  else if (pid == -3){
    Rp = pionm_p;
    Rtheta = pionm_theta;
    Rphi = pionm_phi;
    Rz = pionm_z;
  }
  else if (pid == 4){
    Rp = kaonp_p;
    Rtheta = kaonp_theta;
    Rphi = kaonp_phi;
    Rz = kaonp_z;
  }
  else if (pid == -4){
    Rp = kaonm_p;
    Rtheta = kaonm_theta;
    Rphi = kaonm_phi;
    Rz = kaonm_z;
  }
  else if (pid == 5){
    Rp = proton_p;
    Rtheta = proton_theta;
    Rphi = proton_phi;
    Rz = proton_z;
  }
  else {
    std::cout << "Error in Lresolution::GetResolution!" << std::endl;
    std::cout << "  No pid is matched." << std::endl;
    return 1;
  }
  if (kin[0] <= 0){
    std::cout << "Error in Lresolution::GetResolution!" << std::endl;
    std::cout << "  Nonphysical momentum." << std::endl;
    return 1;
  }
  if (kin[1] > 180.0 || kin[1] < 0.0){
    std::cout << "Error in Lresolution::GetResolution!" << std::endl;
    std::cout << "  Polar angle out of range [0, 180]." << std::endl;
    return 1;
  }
  if (CheckKinematic(kin, pid)) {return 2;}  
  const double radtodeg = 180.0 / TMath::Pi();
  double res_factor=1.5;
  res[0] = res_factor*Rp->GetBinContent(Rp->GetXaxis()->FindBin(kin[0]), Rp->GetYaxis()->FindBin(kin[1])) / 100.0;
  res[1] = res_factor*Rtheta->GetBinContent(Rtheta->GetXaxis()->FindBin(kin[0]), Rtheta->GetYaxis()->FindBin(kin[1])) / 1000.0 * radtodeg;
  res[2] = res_factor*Rphi->GetBinContent(Rphi->GetXaxis()->FindBin(kin[0]), Rphi->GetYaxis()->FindBin(kin[1])) / 1000.0 * radtodeg;
  res[3] = res_factor*Rz->GetBinContent(Rz->GetXaxis()->FindBin(kin[0]), Rz->GetYaxis()->FindBin(kin[1]));
  return 0;
}
  
int Lresolution::GetResolution(const double * kin, double * res, const char pname[]){
  if (CheckConfiguration()) {return 1;}
  int pid;
  if (strcmp(pname, "e+") == 0) pid = 1;
  else if (strcmp(pname, "e-") == 0) pid = -1;
  else if (strcmp(pname, "mu+") == 0) pid = 2;
  else if (strcmp(pname, "mu-") == 0) pid = -2;
  else if (strcmp(pname, "pi+") == 0) pid = 3;
  else if (strcmp(pname, "pi-") == 0) pid = -3;
  else if (strcmp(pname, "K+") == 0) pid = 4;
  else if (strcmp(pname, "K-") == 0) pid = -4;
  else if (strcmp(pname, "p") == 0) pid = 5;
  else pid = -999;
  return GetResolution(kin, res, pid);
}

int Lresolution::PrintResolution(const double * kin, const int pid){
  if (CheckConfiguration()) {return 1;}
  double res[4];
  int rt = GetResolution(kin, res, pid);
  std::cout << "dp/p:     " << res[0] << " " << std::endl;
  std::cout << "dtheta: " << res[1] << " Deg" << std::endl;
  std::cout << "dphi:   " << res[2] << " Deg" << std::endl;
  std::cout << "dz:     " << res[3] << " cm"  << std::endl;
  return rt;
}

int Lresolution::PrintResolution(const double * kin, const char pname[]){
  if (CheckConfiguration()) {return 1;}
  double res[4];
  int rt = GetResolution(kin, res, pname);
  std::cout << "dp/p:     " << res[0] << " " << std::endl;
  std::cout << "dtheta: " << res[1] << " Deg" << std::endl;
  std::cout << "dphi:   " << res[2] << " Deg" << std::endl;
  std::cout << "dz:     " << res[3] << " cm"  << std::endl;
  return rt;
}


#endif
