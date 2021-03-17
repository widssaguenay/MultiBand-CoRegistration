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

ImageGeometry::ImageGeometry(int Cnv, int Rnv)
{

	Cn = Cnv;
	Rn = Rnv;
	UndistortionX.setOnes(1, 1);
	UndistortionY.setOnes(1, 1);

};

ImageGeometry::~ImageGeometry() {
	//this destructs the object
};


void ImageGeometry::Undistortion_mat(CameraParam& im_cam, double model_id, double sensor_id, char* FeatSavePP)
{
	cout << "Creating undistortion look-up tables ..." << endl;

	MatrixXd UndistortionXd; 
	MatrixXd UndistortionYd;

	int retvalx=1;
	char currentfilex[512] = "";
	currentfilex[0] = '\0';
	strcat(currentfilex, FeatSavePP);
	string samplename = "UndistortionX_M" + to_string((int)model_id)+ "_"+ to_string((int)sensor_id) + ".result";
	char* samplenamef = new char[samplename.length() + 1];
	strcpy(samplenamef, samplename.c_str());
	strcat(currentfilex, samplenamef);
	

	int retvaly=1;
	char currentfiley[512] = "";
	currentfiley[0] = '\0';
	strcat(currentfiley, FeatSavePP);
	string samplename3 = "UndistortionY_M" + to_string((int)model_id) + "_" + to_string((int)sensor_id) + ".result";
	char* samplenamef3 = new char[samplename3.length() + 1];
	strcpy(samplenamef3, samplename3.c_str());
	strcat(currentfiley, samplenamef3);
	
	ifstream matfile;
	matfile.open(currentfilex, ios::in);
	if (matfile.fail()) { retvalx = 0; }
	matfile.close();
	matfile.open(currentfiley, ios::in); 
	if (matfile.fail()) { retvaly = 0; }
	matfile.close();



	if (retvalx == 1 && retvaly==1)
	{
		cout << "Success in search for the file path of : " << currentfilex << endl;
		cout << "Success in search for the file path of : " << currentfiley << endl;
		cout << "We will use these files for undistorting images" << endl;
		cout << "Wait for reading these files ... " << endl;
		Read_Mat(currentfilex, UndistortionXd);
		UndistortionX = UndistortionXd.cast<float>();
		UndistortionXd.resize(1, 1);
		Read_Mat(currentfiley, UndistortionYd);
		UndistortionY = UndistortionYd.cast<float>();
		UndistortionYd.resize(1, 1);
		cout << "Ended reading the files! " << endl;
		return;
	}
	if (model_id==0) {

		UndistortionXd.setOnes(Rn, Cn); UndistortionXd = UndistortionXd * MYNAN;
		UndistortionYd.setOnes(Rn, Cn); UndistortionYd = UndistortionYd * MYNAN;

		double K1 = im_cam.K1_0(sensor_id, 0);
		double K2 = im_cam.K2_0(sensor_id, 0);
		double K3 = im_cam.K3_0(sensor_id, 0);
		double K4 = im_cam.K4_0(sensor_id, 0);
		double K5 = im_cam.K5_0(sensor_id, 0);
		double P1 = im_cam.P1_0(sensor_id, 0);
		double P2 = im_cam.P2_0(sensor_id, 0);
		double S1 = im_cam.S1_0(sensor_id, 0);
		double S2 = im_cam.S2_0(sensor_id, 0);

		double PS = im_cam.PS_0(sensor_id, 0);
		double Cx = (-Cn / 2 + 0.5)*PS - im_cam.xpp_0(sensor_id, 0);
		double Cy = -(-Rn / 2 + 0.5)*PS - im_cam.ypp_0(sensor_id, 0);
		double f = im_cam.f_l_0(sensor_id, 0);
		double PSx = PS; 
		double PSy = PS;
		MatrixXd A(3, 3);
		A<< PSx ,0.0 ,Cx,
			0.0, - PSy ,Cy,
			0.0, 0.0, - f;
		MatrixXd Kmat1 = (A.inverse());
		

		int lastcount = 3;
		if (abs(K1) >0.001) {
			lastcount = 20;
		}

		int ystart = 1; int yend = Rn + 1; 
		Concurrency::parallel_for(ystart, yend, [&](int y)
		{
			for (int x = 1; x <= Cn; x++) {
				MatrixXd pt(3, 1);
				pt << (double)x, (double)y, 1.0;

				MatrixXd tmp = A * pt; 
				double xdc = (tmp(0,0) / tmp(2,0))*(A(2, 2));
				double ydc = (tmp(1,0) / tmp(2,0))*(A(2, 2));

				MatrixXd L(2, 1);
				L<<xdc, ydc;

				double xd = xdc; 
				double yd = ydc;

				for (int count = 1; count <= lastcount; count++) {

					double xdd = xd;
					double ydd = yd;
					double r2 = xdd*xdd + ydd*ydd;

					MatrixXd F0(2, 1);
					F0<<(xdd + xdd * (K1*r2 + K2 * r2 * r2 + K3 * r2 * r2 * r2 + K4 * r2 * r2 * r2 * r2 + K5 * r2 * r2 * r2 * r2 * r2) + P1 * (r2 + 2 * xdd * xdd) + 2 * P2*xdd*ydd + S1 * xdd + S2 * ydd),
						(ydd + ydd * (K1*r2 + K2 * r2 * r2 + K3 * r2 * r2 * r2 + K4 * r2 * r2 * r2 * r2 + K5 * r2 * r2 * r2 * r2 * r2) + P2 * (r2 + 2 * ydd * ydd) + 2 * P1*xdd*ydd);

					MatrixXd deltaL = L - F0;

					double r2prime_x = 2 * (xdd);
					double r2prime_y = 2 * (ydd);

					MatrixXd A_mat(2, 2);
					A_mat(0, 0) = S1 + K2 * r2 * r2 + K3 * r2 * r2 * r2 + K4 * r2 * r2 * r2 * r2 + K5 * r2 * r2 * r2 * r2 * r2 + 6 * P1*xdd + 2 * P2*ydd + xdd * (2 * K1*xdd + 4 * K2*xdd*r2 + 6 * K3*xdd*r2 * r2 + 8 * K4*xdd*r2 * r2 * r2 + 10 * K5*xdd*r2 * r2 * r2 * r2) + K1 * r2 + 1;
					A_mat(0, 1) = S2 + 2 * P2*xdd + 2 * P1*ydd + xdd * (2 * K1*ydd + 4 * K2*ydd*r2 + 6 * K3*ydd*r2 * r2 + 8 * K4*ydd*r2 * r2 * r2 + 10 * K5*ydd*r2 * r2 * r2 * r2);
					A_mat(1, 0) = 2 * P2*xdd + 2 * P1*ydd + ydd * (2 * K1*xdd + 4 * K2*xdd*r2 + 6 * K3*xdd*r2 * r2 + 8 * K4*xdd*r2 * r2 * r2 + 10 * K5*xdd*r2 * r2 * r2 * r2);
					A_mat(1, 1) = K2 * r2 * r2 + K3 * r2 * r2 * r2 + K4 * r2 * r2 * r2 * r2 + K5 * r2 * r2 * r2 * r2 * r2 + 2 * P1*xdd + 6 * P2*ydd + ydd * (2 * K1*ydd + 4 * K2*ydd*r2 + 6 * K3*ydd*r2 * r2 + 8 * K4*ydd*r2 * r2 * r2 + 10 * K5*ydd*r2 * r2 * r2 * r2) + K1 * r2 + 1;
					
					MatrixXd some = A_mat.transpose();
					MatrixXd deltaxcap = (some*A_mat).lu().solve(some*deltaL);

					xd = xd + deltaxcap(0,0);
					yd = yd + deltaxcap(1,0);

				}

				MatrixXd pt_c(3,1);
				pt_c << xd, yd, A(2, 2);
				MatrixXd pixels = Kmat1 * pt_c;

				UndistortionXd(y-1, x-1) = pixels(0,0) / pixels(2,0);
				UndistortionYd(y-1, x-1) = pixels(1,0) / pixels(2,0);
			}
		});
		
	}


	if (model_id == 1) {

		UndistortionXd.setOnes(Rn, Cn); UndistortionXd = UndistortionXd * MYNAN;
		UndistortionYd.setOnes(Rn, Cn); UndistortionYd = UndistortionYd * MYNAN;

		double fx = im_cam.fx_1(sensor_id, 0);
		double fy = im_cam.fy_1(sensor_id, 0);
		double cx = im_cam.cx_1(sensor_id, 0);
		double cy = im_cam.cy_1(sensor_id, 0);
		double K1 = im_cam.K1_1(sensor_id, 0);
		double K2 = im_cam.K2_1(sensor_id, 0);
		double K3 = im_cam.K3_1(sensor_id, 0);
		double P1 = im_cam.P1_1(sensor_id, 0);
		double P2 = im_cam.P2_1(sensor_id, 0);

		MatrixXd Kmat1(3, 3);
		Kmat1 <<fx, 0.0, cx,
			0.0, fy, cy,
			0.0, 0.0, 1.0;
		MatrixXd A = (Kmat1.inverse());
		
		int ystart = 1; int yend = Rn + 1;
		Concurrency::parallel_for(ystart, yend, [&](int y)
		{
			for (int x = 1; x <= Cn; x++) {
				MatrixXd pt(3, 1);
				pt << (double)x, (double)y, 1.0;

				MatrixXd tmp = A * pt;
				double r = tmp(0, 0);
				double s = tmp(1, 0);
				double q = tmp(2, 0);

				double xprime = r / q;
				double yprime = s / q;

				double r2 = xprime * xprime + yprime * yprime;

				double xzegond = xprime * (1 + K1 * r2 + K2 * r2 * r2 + K3 * r2 * r2 * r2) + 2 * P1 * xprime * yprime + P2 * (r2 + 2 * xprime * xprime);
				double yzegond = yprime * (1 + K1 * r2 + K2 * r2 * r2 + K3 * r2 * r2 * r2) + 2 * P2 * xprime * yprime + P1 * (r2 + 2 * yprime * yprime);

				double xd = (-fx * xzegond + cx);
				double yd = (fy * yzegond + cy);

				UndistortionXd(y - 1, x - 1) = xd;
				UndistortionYd(y - 1, x - 1) = yd;
			}
		});

	}


	if (model_id == 2) {

		UndistortionXd.setOnes(Rn, Cn); UndistortionXd = UndistortionXd * MYNAN;
		UndistortionYd.setOnes(Rn, Cn); UndistortionYd = UndistortionYd * MYNAN;

		double cx_2 = im_cam.cx_2(sensor_id, 0);
		double cy_2 = im_cam.cy_2(sensor_id, 0);
		double C_2 = im_cam.C_2(sensor_id, 0);
		double D_2 = im_cam.D_2(sensor_id, 0);
		double E_2 = im_cam.E_2(sensor_id, 0);
		double F_2 = im_cam.F_2(sensor_id, 0);
		double P2_2 = im_cam.P2_2(sensor_id, 0);
		double P3_2 = im_cam.P3_2(sensor_id, 0);
		double P4_2 = im_cam.P4_2(sensor_id, 0);

		MatrixXd Kmat1(3, 3);
		double f = 2 * C_2 / pi;
		Kmat1 << -f, 0.0, cx_2,
			0.0, f, cy_2,
			0.0, 0.0, 1.0;
		MatrixXd A = (Kmat1.inverse());

		int ystart = 1; int yend = Rn + 1;
		Concurrency::parallel_for(ystart, yend, [&](int y)
		{
			for (int x = 1; x <= Cn; x++) {
				MatrixXd pt(3, 1);
				pt << (double)x, (double)y, 1.0;

				MatrixXd tmp = A * pt;
				double r = tmp(0, 0);
				double s = tmp(1, 0);
				double q = tmp(2, 0);

				double theta = 2 / pi * atan(sqrt((r * r + s * s)) / q);

				double rho = theta + P2_2 * theta * theta + P3_2 * pow(theta, 3.0) + P4_2 * pow(theta, 4.0);

				double xh = rho * r / sqrt(r * r + s * s);
				double yh = rho * s / sqrt(r * r + s * s);

				double xd = (-C_2 * xh + D_2 * yh + cx_2);
				double yd = (E_2 * xh + F_2 * yh + cy_2);

				UndistortionXd(y - 1, x - 1) = xd;
				UndistortionYd(y - 1, x - 1) = yd;
			}
		});

	}


	if (model_id == 3) {

		UndistortionXd.setOnes(Rn, Cn); UndistortionXd = UndistortionXd * MYNAN;
		UndistortionYd.setOnes(Rn, Cn); UndistortionYd = UndistortionYd * MYNAN;


		double PS = im_cam.PS_3(sensor_id, 0);
		double f = im_cam.f_l_3(sensor_id, 0);
		double K1 = im_cam.K1_3(sensor_id, 0);
		double K2 = im_cam.K2_3(sensor_id, 0);
		double K3 = im_cam.K3_3(sensor_id, 0);
		double K4 = im_cam.K4_3(sensor_id, 0);
		double K5 = im_cam.K5_3(sensor_id, 0);
		double P1 = im_cam.P1_3(sensor_id, 0);
		double P2 = im_cam.P2_3(sensor_id, 0);
		double S1 = im_cam.S1_3(sensor_id, 0);
		double S2 = im_cam.S2_3(sensor_id, 0);

		double Cx = (-Cn / 2 + 0.5) * PS - im_cam.xpp_3(sensor_id, 0);
		double Cy = -(-Rn / 2 + 0.5) * PS - im_cam.ypp_3(sensor_id, 0);
		double PSx = PS;
		double PSy = PS;
		MatrixXd A(3, 3);
		A << PSx, 0.0, Cx,
			0.0, -PSy, Cy,
			0.0, 0.0, -f;
		MatrixXd Kmat1 = (A.inverse());
		Kmat1 = Kmat1 / Kmat1(2, 2);
		A = (Kmat1.inverse());

		int lastcount = 3;
		if (abs(K1) > 0.001) {
			lastcount = 20;
		}

		int ystart = 1; int yend = Rn + 1;
		Concurrency::parallel_for(ystart, yend, [&](int y)
		{
			for (int x = 1; x <= Cn; x++) {
				MatrixXd pt(3, 1);
				pt << (double)x, (double)y, 1.0;

				MatrixXd tmp = A * pt;
				double r = tmp(0, 0);
				double s = tmp(1, 0);
				double q = tmp(2, 0);

				double xdc = f * r * atan(sqrt((r * r + s * s)) / q) / sqrt(r * r + s * s);
				double ydc = f * s * atan(sqrt((r * r + s * s)) / q) / sqrt(r * r + s * s);

				MatrixXd L(2, 1);
				L << xdc, ydc;

				double xd = xdc;
				double yd = ydc;

				for (int count = 1; count <= lastcount; count++) {

					double xdd = xd;
					double ydd = yd;
					double r2 = xdd * xdd + ydd * ydd;

					MatrixXd F0(2, 1);
					F0 << (xdd + xdd * (K1 * r2 + K2 * r2 * r2 + K3 * r2 * r2 * r2 + K4 * r2 * r2 * r2 * r2 + K5 * r2 * r2 * r2 * r2 * r2) + P1 * (r2 + 2 * xdd * xdd) + 2 * P2 * xdd * ydd + S1 * xdd + S2 * ydd),
						(ydd + ydd * (K1 * r2 + K2 * r2 * r2 + K3 * r2 * r2 * r2 + K4 * r2 * r2 * r2 * r2 + K5 * r2 * r2 * r2 * r2 * r2) + P2 * (r2 + 2 * ydd * ydd) + 2 * P1 * xdd * ydd);

					MatrixXd deltaL = L - F0;

					double r2prime_x = 2 * (xdd);
					double r2prime_y = 2 * (ydd);

					MatrixXd A_mat(2, 2);
					A_mat(0, 0) = S1 + K2 * r2 * r2 + K3 * r2 * r2 * r2 + K4 * r2 * r2 * r2 * r2 + K5 * r2 * r2 * r2 * r2 * r2 + 6 * P1 * xdd + 2 * P2 * ydd + xdd * (2 * K1 * xdd + 4 * K2 * xdd * r2 + 6 * K3 * xdd * r2 * r2 + 8 * K4 * xdd * r2 * r2 * r2 + 10 * K5 * xdd * r2 * r2 * r2 * r2) + K1 * r2 + 1;
					A_mat(0, 1) = S2 + 2 * P2 * xdd + 2 * P1 * ydd + xdd * (2 * K1 * ydd + 4 * K2 * ydd * r2 + 6 * K3 * ydd * r2 * r2 + 8 * K4 * ydd * r2 * r2 * r2 + 10 * K5 * ydd * r2 * r2 * r2 * r2);
					A_mat(1, 0) = 2 * P2 * xdd + 2 * P1 * ydd + ydd * (2 * K1 * xdd + 4 * K2 * xdd * r2 + 6 * K3 * xdd * r2 * r2 + 8 * K4 * xdd * r2 * r2 * r2 + 10 * K5 * xdd * r2 * r2 * r2 * r2);
					A_mat(1, 1) = K2 * r2 * r2 + K3 * r2 * r2 * r2 + K4 * r2 * r2 * r2 * r2 + K5 * r2 * r2 * r2 * r2 * r2 + 2 * P1 * xdd + 6 * P2 * ydd + ydd * (2 * K1 * ydd + 4 * K2 * ydd * r2 + 6 * K3 * ydd * r2 * r2 + 8 * K4 * ydd * r2 * r2 * r2 + 10 * K5 * ydd * r2 * r2 * r2 * r2) + K1 * r2 + 1;

					MatrixXd some = A_mat.transpose();
					MatrixXd deltaxcap = (some * A_mat).lu().solve(some * deltaL);

					xd = xd + deltaxcap(0, 0);
					yd = yd + deltaxcap(1, 0);

				}

				MatrixXd pt_c(3, 1);
				pt_c << xd, yd, 1.0;
				MatrixXd pixels = Kmat1 * pt_c;

				UndistortionXd(y - 1, x - 1) = pixels(0, 0) / pixels(2, 0);
				UndistortionYd(y - 1, x - 1) = pixels(1, 0) / pixels(2, 0);
			}
		});

	}


	if (model_id == 4) {

		UndistortionXd.setOnes(Rn, Cn); UndistortionXd = UndistortionXd * MYNAN;
		UndistortionYd.setOnes(Rn, Cn); UndistortionYd = UndistortionYd * MYNAN;

		double K1 = im_cam.K1_4(sensor_id, 0);
		double K2 = im_cam.K2_4(sensor_id, 0);
		double K3 = im_cam.K3_4(sensor_id, 0);
		double P1 = im_cam.P1_4(sensor_id, 0);
		double P2 = im_cam.P2_4(sensor_id, 0);

		double PS = im_cam.PS_4(sensor_id, 0);
		double Cx = im_cam.Cx_4(sensor_id, 0);
		double Cy = im_cam.Cy_4(sensor_id, 0);
		double f = im_cam.f_l_4(sensor_id, 0);
		double PSx = PS;
		double PSy = PS;
		MatrixXd A(3, 3);
		A << PSx, 0.0, -Cx,
			0.0, -PSy, Cy,
			0.0, 0.0, -f;
		MatrixXd Kmat1 = (A.inverse());

		int ystart = 1; int yend = Rn + 1;
		Concurrency::parallel_for(ystart, yend, [&](int y)
		{
			for (int x = 1; x <= Cn; x++) {
				MatrixXd pt(3, 1);
				pt << (double)x, (double)y, 1.0; 

				MatrixXd tmp = A * pt; 
				double r = tmp(0, 0);
				double s = tmp(1, 0);
				double q = tmp(2, 0);

				double xprime = r / q;
				double yprime = s / q;

				double r2 = xprime * xprime + yprime * yprime;

				double xzegond = xprime * (1 + K1 * r2 + K2 * r2*r2 + K3* r2*r2*r2) + 2 * P1*xprime*yprime + P2* (r2 + 2 * xprime*xprime);
				double yzegond = yprime * (1 + K1* r2 + K2 * r2*r2 + K3* r2*r2*r2) + 2 * P2*xprime*yprime + P1* (r2 + 2 * yprime*yprime);

				double xd = (-f * xzegond + Cx) / PS;
				double yd = (f * yzegond + Cy) / PS;

				UndistortionXd(y - 1, x - 1) = xd;
				UndistortionYd(y - 1, x - 1) = yd;
			}
		});

	}


	if (model_id == 6) {

		UndistortionXd.setOnes(Rn, Cn); UndistortionXd = UndistortionXd * MYNAN;
		UndistortionYd.setOnes(Rn, Cn); UndistortionYd = UndistortionYd * MYNAN;

		double PS = im_cam.PS_00(sensor_id, 0);
		double f = im_cam.f_l_00(sensor_id, 0);
		double xpp = im_cam.xpp_00(sensor_id, 0);
		double ypp = im_cam.ypp_00(sensor_id, 0);
		double K1 = im_cam.K1_00(sensor_id, 0);
		double K2 = im_cam.K2_00(sensor_id, 0);
		double K3 = im_cam.K3_00(sensor_id, 0);
		double K4 = im_cam.K4_00(sensor_id, 0);
		double K5 = im_cam.K5_00(sensor_id, 0);
		double P1 = im_cam.P1_00(sensor_id, 0);
		double P2 = im_cam.P2_00(sensor_id, 0);
		double S1 = im_cam.SS1_00(sensor_id, 0);
		double S2 = im_cam.SS2_00(sensor_id, 0);
		double S3 = im_cam.SS3_00(sensor_id, 0);
		double S4 = im_cam.SS4_00(sensor_id, 0);
		double S5 = im_cam.SS5_00(sensor_id, 0);

		double Cx = (-Cn / 2 + 0.5) * PS - xpp;
		double Cy = -(-Rn / 2 + 0.5) * PS - ypp;
		double PSx = PS;
		double PSy = PS;
		MatrixXd A(3, 3);
		A << PSx, 0.0, Cx,
			0.0, -PSy, Cy,
			0.0, 0.0, -f;
		MatrixXd Kmat1 = (A.inverse());


		int lastcount = 3;
		if (abs(K1) > 0.001) {
			lastcount = 20;
		}

		int ystart = 1; int yend = Rn + 1;
		Concurrency::parallel_for(ystart, yend, [&](int y)
		{
			for (int x = 1; x <= Cn; x++) {
				MatrixXd pt(3, 1);
				pt << (double)x, (double)y, 1.0;

				MatrixXd tmp = A * pt;
				double xdc = (tmp(0, 0) / tmp(2, 0)) * (A(2, 2));
				double ydc = (tmp(1, 0) / tmp(2, 0)) * (A(2, 2));

				MatrixXd L(2, 1);
				L << xdc, ydc;

				double xd = xdc;
				double yd = ydc;

				for (int count = 1; count <= lastcount; count++) {

					double xdd = xd;
					double ydd = yd;
					double r2 = xdd * xdd + ydd * ydd;

					MatrixXd F0(2, 1);
					F0 << (xdd + xdd * (K1 * r2 + K2 * r2 * r2 + K3 * r2 * r2 * r2 + K4 * r2 * r2 * r2 * r2 + K5 * r2 * r2 * r2 * r2 * r2) + P1 * (r2 + 2 * xdd * xdd) + 2 * P2 * xdd * ydd + S1 * xdd + S4 * ydd + S2 * (xdd * xdd) + S3 * (xdd * xdd * xdd) + S5 * (xdd * xdd * xdd * xdd)),
						(ydd + ydd * (K1 * r2 + K2 * r2 * r2 + K3 * r2 * r2 * r2 + K4 * r2 * r2 * r2 * r2 + K5 * r2 * r2 * r2 * r2 * r2) + P2 * (r2 + 2 * ydd * ydd) + 2 * P1 * xdd * ydd);

					MatrixXd deltaL = L - F0;

					double r2prime_x = 2 * (xdd);
					double r2prime_y = 2 * (ydd);

					MatrixXd A_mat(2, 2);
					A_mat(0, 0) = S1 + 2 * xdd * S2 + 3 * (xdd * xdd) * S3 + 4 * (xdd * xdd * xdd) * S5 + K2 * r2 * r2 + K3 * r2 * r2 * r2 + K4 * r2 * r2 * r2 * r2 + K5 * r2 * r2 * r2 * r2 * r2 + 6 * P1 * xdd + 2 * P2 * ydd + xdd * (2 * K1 * xdd + 4 * K2 * xdd * r2 + 6 * K3 * xdd * r2 * r2 + 8 * K4 * xdd * r2 * r2 * r2 + 10 * K5 * xdd * r2 * r2 * r2 * r2) + K1 * r2 + 1;
					A_mat(0, 1) = S4 + 2 * P2 * xdd + 2 * P1 * ydd + xdd * (2 * K1 * ydd + 4 * K2 * ydd * r2 + 6 * K3 * ydd * r2 * r2 + 8 * K4 * ydd * r2 * r2 * r2 + 10 * K5 * ydd * r2 * r2 * r2 * r2);
					A_mat(1, 0) = 2 * P2 * xdd + 2 * P1 * ydd + ydd * (2 * K1 * xdd + 4 * K2 * xdd * r2 + 6 * K3 * xdd * r2 * r2 + 8 * K4 * xdd * r2 * r2 * r2 + 10 * K5 * xdd * r2 * r2 * r2 * r2);
					A_mat(1, 1) = K2 * r2 * r2 + K3 * r2 * r2 * r2 + K4 * r2 * r2 * r2 * r2 + K5 * r2 * r2 * r2 * r2 * r2 + 2 * P1 * xdd + 6 * P2 * ydd + ydd * (2 * K1 * ydd + 4 * K2 * ydd * r2 + 6 * K3 * ydd * r2 * r2 + 8 * K4 * ydd * r2 * r2 * r2 + 10 * K5 * ydd * r2 * r2 * r2 * r2) + K1 * r2 + 1;

					MatrixXd some = A_mat.transpose();
					MatrixXd deltaxcap = (some * A_mat).lu().solve(some * deltaL);

					xd = xd + deltaxcap(0, 0);
					yd = yd + deltaxcap(1, 0);

				}

				MatrixXd pt_c(3, 1);
				pt_c << xd, yd, A(2, 2);
				MatrixXd pixels = Kmat1 * pt_c;

				UndistortionXd(y - 1, x - 1) = pixels(0, 0) / pixels(2, 0);
				UndistortionYd(y - 1, x - 1) = pixels(1, 0) / pixels(2, 0);
			}
		});

	}


	cout << "Ended creating undistortion look-up tables!" << endl;

	cout << "Wait for writing these tables out ..." << endl;
	Write_Mat(currentfilex, UndistortionXd, 5);
	UndistortionX = UndistortionXd.cast<float>();
	UndistortionXd.resize(1, 1);
	Write_Mat(currentfiley, UndistortionYd, 5);
	UndistortionY = UndistortionYd.cast<float>();
	UndistortionYd.resize(1, 1);

	cout << "Ended writing these tables!" << endl;

};

void ImageGeometry::Undistort(char* FeatSavePP, string imname, CameraParam& im_cam, int im, string imformat) {

	int sensor_id = (int)im_cam.sensor_ID(im, 0);
	int model_id = (int)im_cam.model_ID(im, 0);

	string foldername = string(FeatSavePP);
	string imfullname = foldername + imname;
	Mat Image = imread(imfullname, IMREAD_COLOR | IMREAD_ANYDEPTH | IMREAD_IGNORE_ORIENTATION);
	if (Image.empty())
	{
		cout << "Cannot read image: " << imfullname<< std::endl;
		return;
	}
	//imshow("Display Window", Image);
	//waitKey(0);

	vector<int> imwrite_jpg_params;
	imwrite_jpg_params.push_back(IMWRITE_JPEG_QUALITY);
	imwrite_jpg_params.push_back(100);

	int found = imname.find(imformat);
	string newname = foldername + "Undist_" + (imname.substr(0, found)); 
	newname += ".tiff";

	int retval = 1;
	ifstream matfile;
	matfile.open(newname, ios::in);
	if (matfile.fail()) { retval = 0; }
	matfile.close();
	
	if (retval == 1)
	{
		cout << "Success in search for the file path of : " << newname << endl;
		return;
	}

	
	if (UndistortionX.rows() == 1 && UndistortionY.cols() == 1) {
		Undistortion_mat(im_cam, model_id, sensor_id, FeatSavePP);
	}

	cout << "Wait for undistorting image " << imname << " ..." << endl;
	Mat Imcorrect(Image.size(), Image.type());

	Mat X(Image.size(), CV_32FC1);
	Mat Y(Image.size(), CV_32FC1);
	eigen2cv(UndistortionX, X);
	eigen2cv(UndistortionY, Y);


	
	remap(Image, Imcorrect, X, Y, INTER_CUBIC, 0, Scalar(255, 255, 255));
	cout << "Ended undistorting image " << imname << " !" << endl;

	X.release();
	Y.release();

	imwrite(newname, Imcorrect);//, imwrite_jpg_params

};

/////////////////////////////////////
/////////////////////////////////////
