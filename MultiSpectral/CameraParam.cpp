/*
//**********************************************************************
This program is developed by Mozhdeh Shahbazi
at Centre de geomatique du Quebec.
This program is distributed under Creative Commons
Attribution 4.0 International (CC BY 4.0) license.
Under this license, You are free to:
Share — copy and redistribute the material in any medium or format
Adapt — remix, transform, and build upon the material
for any purpose, even commercially.
//**********************************************************************
*/
#include "Coregistration.h"

CameraParam::CameraParam(int nsensorwithmodelid0, int nsensorwithmodelid1,
	int nsensorwithmodelid2, int nsensorwithmodelid3, int nsensorwithmodelid4, int nsensorwithmodelid6) {

	if (nsensorwithmodelid6 > 0) {
		//for traditional perspective pinhole model of Zed camera
		f_l_00.resize(nsensorwithmodelid6, 1);
		xpp_00.resize(nsensorwithmodelid6, 1);
		ypp_00.resize(nsensorwithmodelid6, 1);
		K1_00.resize(nsensorwithmodelid6, 1);
		K2_00.resize(nsensorwithmodelid6, 1);
		K3_00.resize(nsensorwithmodelid6, 1);
		K4_00.resize(nsensorwithmodelid6, 1);
		K5_00.resize(nsensorwithmodelid6, 1);
		P1_00.resize(nsensorwithmodelid6, 1);
		P2_00.resize(nsensorwithmodelid6, 1);
		SS5_00.resize(nsensorwithmodelid6, 1);
		SS4_00.resize(nsensorwithmodelid6, 1);
		SS3_00.resize(nsensorwithmodelid6, 1);
		SS2_00.resize(nsensorwithmodelid6, 1);
		SS1_00.resize(nsensorwithmodelid6, 1);
		PS_00.resize(nsensorwithmodelid6, 1);
		Cn_00.resize(nsensorwithmodelid6, 1);
		Rn_00.resize(nsensorwithmodelid6, 1);
		Kcoefs_00.resize(nsensorwithmodelid6, 1);
		Calibrated_6.resize(nsensorwithmodelid6, 1);
	}

	if (nsensorwithmodelid0 > 0) {
		//for traditional perspective pinhole model
		f_l_0.resize(nsensorwithmodelid0, 1);
		xpp_0.resize(nsensorwithmodelid0, 1);
		ypp_0.resize(nsensorwithmodelid0, 1);
		K1_0.resize(nsensorwithmodelid0, 1);
		K2_0.resize(nsensorwithmodelid0, 1);
		K3_0.resize(nsensorwithmodelid0, 1);
		K4_0.resize(nsensorwithmodelid0, 1);
		K5_0.resize(nsensorwithmodelid0, 1);
		P1_0.resize(nsensorwithmodelid0, 1);
		P2_0.resize(nsensorwithmodelid0, 1);
		S2_0.resize(nsensorwithmodelid0, 1);
		S1_0.resize(nsensorwithmodelid0, 1);
		PS_0.resize(nsensorwithmodelid0, 1);
		Cn_0.resize(nsensorwithmodelid0, 1);
		Rn_0.resize(nsensorwithmodelid0, 1);
		Kcoefs_0.resize(nsensorwithmodelid0, 1);
		Calibrated_0.resize(nsensorwithmodelid0, 1);
	}
	if (nsensorwithmodelid1 > 0) {
		//for opencv pinhole model
		Cn_1.resize(nsensorwithmodelid1, 1);
		Rn_1.resize(nsensorwithmodelid1, 1);
		fx_1.resize(nsensorwithmodelid1, 1);
		fy_1.resize(nsensorwithmodelid1, 1);
		cx_1.resize(nsensorwithmodelid1, 1);
		cy_1.resize(nsensorwithmodelid1, 1);
		K1_1.resize(nsensorwithmodelid1, 1);
		K2_1.resize(nsensorwithmodelid1, 1);
		K3_1.resize(nsensorwithmodelid1, 1);
		P1_1.resize(nsensorwithmodelid1, 1);
		P2_1.resize(nsensorwithmodelid1, 1);
		Calibrated_1.resize(nsensorwithmodelid1, 1);
	}

	if (nsensorwithmodelid2 > 0) {
		//forPix4D fisheye model
		Cn_2.resize(nsensorwithmodelid2, 1);
		Rn_2.resize(nsensorwithmodelid2, 1);
		C_2.resize(nsensorwithmodelid2, 1);
		D_2.resize(nsensorwithmodelid2, 1);
		E_2.resize(nsensorwithmodelid2, 1);
		F_2.resize(nsensorwithmodelid2, 1);
		cx_2.resize(nsensorwithmodelid2, 1);
		cy_2.resize(nsensorwithmodelid2, 1);
		P2_2.resize(nsensorwithmodelid2, 1);
		P3_2.resize(nsensorwithmodelid2, 1);
		P4_2.resize(nsensorwithmodelid2, 1);
		Calibrated_2.resize(nsensorwithmodelid2, 1);
	}

	if (nsensorwithmodelid3 > 0) {
		//for traditional perspective pinhole model
		f_l_3.resize(nsensorwithmodelid3, 1);
		xpp_3.resize(nsensorwithmodelid3, 1);
		ypp_3.resize(nsensorwithmodelid3, 1);
		K1_3.resize(nsensorwithmodelid3, 1);
		K2_3.resize(nsensorwithmodelid3, 1);
		K3_3.resize(nsensorwithmodelid3, 1);
		K4_3.resize(nsensorwithmodelid3, 1);
		K5_3.resize(nsensorwithmodelid3, 1);
		P1_3.resize(nsensorwithmodelid3, 1);
		P2_3.resize(nsensorwithmodelid3, 1);
		S2_3.resize(nsensorwithmodelid3, 1);
		S1_3.resize(nsensorwithmodelid3, 1);
		PS_3.resize(nsensorwithmodelid3, 1);
		Cn_3.resize(nsensorwithmodelid3, 1);
		Rn_3.resize(nsensorwithmodelid3, 1);
		Kcoefs_3.resize(nsensorwithmodelid3, 1);
		Calibrated_3.resize(nsensorwithmodelid3, 1);
	}

	if (nsensorwithmodelid4 > 0) {
		//for pix4d perspective pinhole model
		f_l_4.resize(nsensorwithmodelid4, 1);
		Cx_4.resize(nsensorwithmodelid4, 1);
		Cy_4.resize(nsensorwithmodelid4, 1);
		K1_4.resize(nsensorwithmodelid4, 1);
		K2_4.resize(nsensorwithmodelid4, 1);
		K3_4.resize(nsensorwithmodelid4, 1);
		P1_4.resize(nsensorwithmodelid4, 1);
		P2_4.resize(nsensorwithmodelid4, 1);
		PS_4.resize(nsensorwithmodelid4, 1);
		Cn_4.resize(nsensorwithmodelid4, 1);
		Rn_4.resize(nsensorwithmodelid4, 1);
		Calibrated_4.resize(nsensorwithmodelid4, 1);
	}
	
	num_sensors_model0 = nsensorwithmodelid0;
	num_sensors_model1 = nsensorwithmodelid1;
	num_sensors_model2 = nsensorwithmodelid2;
	num_sensors_model3 = nsensorwithmodelid3;
	num_sensors_model4 = nsensorwithmodelid4;
	num_sensors_model6 = nsensorwithmodelid6;
};

CameraParam::~CameraParam() {
	//this destructs the object
};


void CameraParam::Re_init_Estimori(int i, double Omega, double Phi, double Kappa, double Xo, double Yo, double Zo) {

	int nrow = 5 * i;

	Matrix3b3 R;
	Rotation_g2i(Omega, Phi, Kappa, R);

	Bundestim_ori.row(nrow + 0) << R.row(0);
	Bundestim_ori.row(nrow + 1) << R.row(1);
	Bundestim_ori.row(nrow + 2) << R.row(2);

	Bundestim_ori.row(nrow + 3) << Xo, Yo, Zo;

	Bundestim_ori.row(nrow + 4) << Omega, Phi, Kappa;
};

void CameraParam::ResetEOPs(Matrixdb3& Estimori) {

	int n_cams = (this->Bundestim_ori.rows()) / 5;
	for (int i = 0; i<n_cams; i++) {

		Matrix3b3 R;
		int nrow = i * 5;

		R = Estimori.block(nrow, 0, 3, 3);

		this->Bundestim_ori.block(nrow, 0, 3, 3) = R;

		this->Bundestim_ori.row(nrow + 3) << Estimori(nrow + 3, 0), Estimori(nrow + 3, 1), Estimori(nrow + 3, 2);

		this->Bundestim_ori.row(nrow + 4) << Estimori(nrow + 4, 0), Estimori(nrow + 4, 1), Estimori(nrow + 4, 2);
	}

};


void CameraParam::ResetIOPs0(int i, bool calibrated, double PSv, double Cnv, double Rnv, double K1v, double K2v, double K3v, double K4v, double K5v, double P1v, double P2v, double S1v, double S2v, double xppv, double yppv, double f_lv, double kcoeffsv) {

	f_l_0(i, 0) = f_lv;
	xpp_0(i, 0) = xppv;
	ypp_0(i, 0) = yppv;
	PS_0(i, 0) = PSv;
	Cn_0(i, 0) = Cnv;
	Rn_0(i, 0) = Rnv;
	K1_0(i, 0) = K1v;
	K2_0(i, 0) = K2v;
	K3_0(i, 0) = K3v;
	K4_0(i, 0) = K4v;
	K5_0(i, 0) = K5v;
	P1_0(i, 0) = P1v;
	P2_0(i, 0) = P2v;
	S1_0(i, 0) = S1v;
	S2_0(i, 0) = S2v;
	Kcoefs_0(i, 0) = kcoeffsv;
	Calibrated_0(i, 0) = calibrated;
};


void CameraParam::ResetIOPs1(int i, bool calibrated, double Cnv, double Rnv, double K1v, double K2v, double K3v, double P1v, double P2v, double cxv, double cyv, double fxv, double fyv) {

	Cn_1(i, 0) = Cnv;
	Rn_1(i, 0) = Rnv;
	fx_1(i, 0) = fxv;
	fy_1(i, 0) = fyv;
	cx_1(i, 0) = cxv;
	cy_1(i, 0) = cyv;
	K1_1(i, 0) = K1v;
	K2_1(i, 0) = K2v;
	K3_1(i, 0) = K3v;
	P1_1(i, 0) = P1v;
	P2_1(i, 0) = P2v;
	Calibrated_1(i, 0) = calibrated;
};


void CameraParam::ResetIOPs2(int i, bool calibrated, double Cnv, double Rnv, double P2v, double P3v, double P4v, double cxv, double cyv, double Cv, double Dv, double Fv) {
	//i is the sensor index; must be between 0 to num_sensors_model2-1

	Cn_2(i, 0) = Cnv;
	Rn_2(i, 0) = Rnv;
	cx_2(i, 0) = cxv;
	cy_2(i, 0) = cyv;
	P2_2(i, 0) = P2v;
	P3_2(i, 0) = P3v;
	P4_2(i, 0) = P4v;
	C_2(i, 0) = Cv;
	D_2(i, 0) = Dv;
	F_2(i, 0) = Fv;
	E_2(i, 0) = 0.0;
	Calibrated_2(i, 0) = calibrated;
};

void CameraParam::ResetIOPs3(int i, bool calibrated, double PSv, double Cnv, double Rnv, double K1v, double K2v, double K3v, double K4v, double K5v, double P1v, double P2v, double S1v, double S2v, double xppv, double yppv, double f_lv, double kcoeffsv) {
	//i is the sensor index; must be between 0 to num_sensors_model3-1
	f_l_3(i, 0) = f_lv;
	xpp_3(i, 0) = xppv;
	ypp_3(i, 0) = yppv;
	PS_3(i, 0) = PSv;
	Cn_3(i, 0) = Cnv;
	Rn_3(i, 0) = Rnv;
	K1_3(i, 0) = K1v;
	K2_3(i, 0) = K2v;
	K3_3(i, 0) = K3v;
	K4_3(i, 0) = K4v;
	K5_3(i, 0) = K5v;
	P1_3(i, 0) = P1v;
	P2_3(i, 0) = P2v;
	S1_3(i, 0) = S1v;
	S2_3(i, 0) = S2v;
	Kcoefs_3(i, 0) = kcoeffsv;
	Calibrated_3(i, 0) = calibrated;
};

void CameraParam::ResetIOPs4(int i, bool calibrated, double PSv, double Cnv, double Rnv, double K1v, double K2v, double K3v, double P1v, double P2v, double Cxv, double Cyv, double f_lv) {

	f_l_4(i, 0) = f_lv;
	Cx_4(i, 0) = Cxv;
	Cy_4(i, 0) = Cyv;
	PS_4(i, 0) = PSv;
	Cn_4(i, 0) = Cnv;
	Rn_4(i, 0) = Rnv;
	K1_4(i, 0) = K1v;
	K2_4(i, 0) = K2v;
	K3_4(i, 0) = K3v;
	P1_4(i, 0) = P1v;
	P2_4(i, 0) = P2v;
	Calibrated_4(i, 0) = calibrated;
};

void CameraParam::ResetIOPs6(int i, bool calibrated, double PSv, double Cnv, double Rnv, double K1v, double K2v, double K3v, double K4v, double K5v, double P1v, double P2v, double SS1v, double SS2v, double SS3v, double SS4v, double SS5v, double xppv, double yppv, double f_lv, double kcoeffsv) {
	//i is the sensor index; must be between 0 to num_sensors_model6-1
	f_l_00(i, 0) = f_lv;
	xpp_00(i, 0) = xppv;
	ypp_00(i, 0) = yppv;
	PS_00(i, 0) = PSv;
	Cn_00(i, 0) = Cnv;
	Rn_00(i, 0) = Rnv;
	K1_00(i, 0) = K1v;
	K2_00(i, 0) = K2v;
	K3_00(i, 0) = K3v;
	K4_00(i, 0) = K4v;
	K5_00(i, 0) = K5v;
	P1_00(i, 0) = P1v;
	P2_00(i, 0) = P2v;
	SS1_00(i, 0) = SS1v;
	SS2_00(i, 0) = SS2v;
	SS3_00(i, 0) = SS3v;
	SS4_00(i, 0) = SS4v;
	SS5_00(i, 0) = SS5v;
	Kcoefs_00(i, 0) = kcoeffsv;
	Calibrated_6(i, 0) = calibrated;

};


void CameraParam::WriteCameraToFile(char *FileName, char *FileName2) {

	ofstream matfile;
	matfile.open(FileName, ios::out); 
	if (matfile.fail()) 
	{
		cout << "There was a problem opening the follwoing file: " << endl << FileName << endl;
		system("pause");
		exit(EXIT_FAILURE);


	}
	matfile.flags(ios::fixed);


	if (num_sensors_model0 > 0) {
		matfile << "IOPs of Sensors with Traditional Perspective Pinhole Model " << endl;
		for (int i = 0; i < f_l_0.rows(); i++) {
			matfile << "IOPs of Sensor " << setprecision(0) << i << " : " << endl;
			matfile << setprecision(25) << "PS=" << this->PS_0(i, 0) << endl;
			matfile << setprecision(25) << "Cn=" << this->Cn_0(i, 0) << endl;
			matfile << setprecision(25) << "Rn=" << this->Rn_0(i, 0) << endl;
			matfile << setprecision(25) << "xpp=" << this->xpp_0(i, 0) << endl;
			matfile << setprecision(25) << "ypp=" << this->ypp_0(i, 0) << endl;
			matfile << setprecision(25) << "f=" << this->f_l_0(i, 0) << endl;
			matfile << setprecision(15) << "K1,K2,K3=" << this->K1_0(i, 0) << "," << this->K2_0(i, 0) << "," << this->K3_0(i, 0) << endl;
			matfile << setprecision(15) << "K4,K5=" << this->K4_0(i, 0) << "," << this->K5_0(i, 0) << endl;
			matfile << setprecision(15) << "P1,P2=" << this->P1_0(i, 0) << "," << this->P2_0(i, 0) << endl;
			matfile << setprecision(15) << "S1,S2=" << this->S1_0(i, 0) << "," << this->S2_0(i, 0) << endl;
		}
	}


	if (num_sensors_model1 > 0) {
		matfile << "IOPs of Sensors with OpenCV Perspective Pinhole Model " << endl;
		for (int i = 0; i < fx_1.rows(); i++) {
			matfile << "IOPs of Sensor " << setprecision(0) << i << " : " << endl;
			matfile << setprecision(25) << "cx=" << this->cx_1(i, 0) << endl;
			matfile << setprecision(25) << "cy=" << this->cy_1(i, 0) << endl;
			matfile << setprecision(25) << "fx=" << this->fx_1(i, 0) << endl;
			matfile << setprecision(25) << "fy=" << this->fy_1(i, 0) << endl;
			matfile << setprecision(15) << "K1,K2,K3=" << this->K1_1(i, 0) << "," << this->K2_1(i, 0) << "," << this->K3_1(i, 0) << endl;
			matfile << setprecision(15) << "P1,P2=" << this->P1_1(i, 0) << "," << this->P2_1(i, 0) << endl;
		}
	}

	if (num_sensors_model2 > 0) {
		matfile << "IOPs of Sensors with Pix4D Fisheye Model " << endl;
		for (int i = 0; i < cx_2.rows(); i++) {
			matfile << "IOPs of Sensor " << setprecision(0) << i << " : " << endl;
			matfile << setprecision(25) << "Cn=" << this->Cn_2(i, 0) << endl;
			matfile << setprecision(25) << "Rn=" << this->Rn_2(i, 0) << endl;
			matfile << setprecision(25) << "cx=" << this->cx_2(i, 0) << endl;
			matfile << setprecision(25) << "cy=" << this->cy_2(i, 0) << endl;
			matfile << setprecision(25) << "C,D,E,F=" << this->C_2(i, 0) << "," << this->D_2(i, 0) << "," << this->E_2(i, 0) << "," << this->F_2(i, 0) << endl;
			matfile << setprecision(15) << "P2,P3,P4=" << this->P2_2(i, 0) << "," << this->P3_2(i, 0) << "," << this->P4_2(i, 0) << endl;
		}
	}

	if (num_sensors_model3 > 0) {
		matfile << "IOPs of Sensors with Traditional Fisheye Model " << endl;
		for (int i = 0; i < f_l_3.rows(); i++) {
			matfile << "IOPs of Sensor " << setprecision(0) << i << " : " << endl;
			matfile << setprecision(25) << "PS=" << this->PS_3(i, 0) << endl;
			matfile << setprecision(25) << "Cn=" << this->Cn_3(i, 0) << endl;
			matfile << setprecision(25) << "Rn=" << this->Rn_3(i, 0) << endl;
			matfile << setprecision(25) << "xpp=" << this->xpp_3(i, 0) << endl;
			matfile << setprecision(25) << "ypp=" << this->ypp_3(i, 0) << endl;
			matfile << setprecision(25) << "f=" << this->f_l_3(i, 0) << endl;
			matfile << setprecision(25) << "K1,K2,K3=" << this->K1_3(i, 0) << "," << this->K2_3(i, 0) << "," << this->K3_3(i, 0) << endl;
			matfile << setprecision(25) << "K4,K5=" << this->K4_3(i, 0) << "," << this->K5_3(i, 0) << endl;
			matfile << setprecision(25) << "P1,P2=" << this->P1_3(i, 0) << "," << this->P2_3(i, 0) << endl;
			matfile << setprecision(25) << "S1,S2=" << this->S1_3(i, 0) << "," << this->S2_3(i, 0) << endl;
		}
	}
	
	if (num_sensors_model4 > 0) {
		matfile << "IOPs of Sensors with Pix4D pinhole Model " << endl;
		for (int i = 0; i < f_l_4.rows(); i++) {
			matfile << "IOPs of Sensor " << setprecision(0) << i << " : " << endl;
			matfile << setprecision(25) << "PS=" << this->PS_4(i, 0) << endl;
			matfile << setprecision(25) << "Cn=" << this->Cn_4(i, 0) << endl;
			matfile << setprecision(25) << "Rn=" << this->Rn_4(i, 0) << endl;
			matfile << setprecision(25) << "Cx=" << this->Cx_4(i, 0) << endl;
			matfile << setprecision(25) << "Cy=" << this->Cy_4(i, 0) << endl;
			matfile << setprecision(25) << "f=" << this->f_l_4(i, 0) << endl;
			matfile << setprecision(25) << "K1,K2,K3=" << this->K1_4(i, 0) << "," << this->K2_4(i, 0) << "," << this->K3_4(i, 0) << endl;
			matfile << setprecision(25) << "P1,P2=" << this->P1_4(i, 0) << "," << this->P2_4(i, 0) << endl;
		}
	}

	if (num_sensors_model6 > 0) {
		matfile << "IOPs of ZED Sensors with Traditional Perspective Pinhole Model " << endl;
		for (int i = 0; i < f_l_00.rows(); i++) {
			matfile << "IOPs of Sensor " << setprecision(0) << i << " : " << endl;
			matfile << setprecision(25) << "PS=" << this->PS_00(i, 0) << endl;
			matfile << setprecision(25) << "Cn=" << this->Cn_00(i, 0) << endl;
			matfile << setprecision(25) << "Rn=" << this->Rn_00(i, 0) << endl;
			matfile << setprecision(25) << "xpp=" << this->xpp_00(i, 0) << endl;
			matfile << setprecision(25) << "ypp=" << this->ypp_00(i, 0) << endl;
			matfile << setprecision(25) << "f=" << this->f_l_00(i, 0) << endl;
			matfile << setprecision(15) << "K1,K2,K3=" << this->K1_00(i, 0) << "," << this->K2_00(i, 0) << "," << this->K3_00(i, 0) << endl;
			matfile << setprecision(15) << "K4,K5=" << this->K4_00(i, 0) << "," << this->K5_00(i, 0) << endl;
			matfile << setprecision(15) << "P1,P2=" << this->P1_00(i, 0) << "," << this->P2_00(i, 0) << endl;
			matfile << setprecision(15) << "S1,S2,S3=" << this->SS1_00(i, 0) << "," << this->SS2_00(i, 0) << "," << this->SS3_00(i, 0) << endl;
			matfile << setprecision(15) << "S4,S5=" << this->SS4_00(i, 0) << "," << this->SS5_00(i, 0) << endl;
		}
	}

	int n_cams = (this->Bundestim_ori.rows()) / 5;
	for (int i = 0; i<n_cams; i++) {

		matfile << "EOPs of Image " << setprecision(0) << this->Camera_labels(i) << " :" << endl;
		int nrow = i * 5;

		matfile << setprecision(15) << this->Bundestim_ori.block(nrow, 0, 3, 3) << endl; 

		matfile << setprecision(15) << this->Bundestim_ori.row(nrow + 3) << endl; 

		matfile << setprecision(20) << this->Bundestim_ori.row(nrow + 4) << endl;
	}


	matfile.close();

	
	matfile.open(FileName2, ios::out); 
	if (matfile.fail()) 
	{
		cout << "There was a problem opening the follwoing file: " << endl << FileName2 << endl;
		system("pause");
		exit(EXIT_FAILURE);
	}
	matfile.flags(ios::fixed);


	for (int i = 0; i<n_cams; i++) {

		matfile << setprecision(0) << this->Camera_labels(i) << " ";
		int nrow = i * 5;

		matfile << setprecision(15) << this->Bundestim_ori(nrow + 3, 0) << " " << this->Bundestim_ori(nrow + 3, 1) << " " << this->Bundestim_ori(nrow + 3, 2) << " "; //Xo, Yo, Zo

		matfile << setprecision(30) << this->Bundestim_ori(nrow + 4, 0) * 180 / pi << " " << this->Bundestim_ori(nrow + 4, 1) * 180 / pi << " " << this->Bundestim_ori(nrow + 4, 2) * 180 / pi << endl;//Omega, Phi, Kappa
	}

	matfile.close();
	return;

};