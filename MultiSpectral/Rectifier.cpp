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

bool Rectifier(char* FeatSavePP, string ImageNameL, string ImageNameR, int im1, int im2,
	CameraParam& im_cam,
	int process_sizex, int process_sizey, Matrix3b3& R_common,
	MatrixXd& K_common_new1, MatrixXd& K_common_new2, vector<double> scale_factors, bool do_resize, bool process_in_batch,
	Matrix<int, Dynamic, Dynamic>& Mask,
	string imformat,
	MatrixXd& H1i_transform, MatrixXd& H2i_transform,
	Matrix3b4& Pmat1, Matrix3b4& Pmat2,
	bool is_coregistration) {

	string foldername = string(FeatSavePP);

	vector<int> imwrite_pgm_params;
	imwrite_pgm_params.push_back(IMWRITE_PXM_BINARY);
	imwrite_pgm_params.push_back(1);


	int found = ImageNameL.find(imformat);
	string newnameL = foldername + "Undist_" + (ImageNameL.substr(0, found));
	newnameL += ".tiff";
	Mat ImageL = imread(newnameL, IMREAD_GRAYSCALE);
	if (ImageL.empty())
	{
		std::cout << "Cannot read image: " << newnameL << std::endl;
		return false;
	}


	found = ImageNameR.find(imformat);
	string newnameR = foldername + "Undist_" + (ImageNameR.substr(0, found));
	newnameR += ".tiff";
	Mat ImageR = imread(newnameR, IMREAD_GRAYSCALE);
	if (ImageR.empty())
	{
		std::cout << "Cannot read image: " << newnameR << std::endl;
		return false;
	}

	found = ImageNameL.find('.');
	string leftname = (ImageNameL.substr(0, found));
	found = ImageNameR.find('.');
	string rightname = (ImageNameR.substr(0, found));

	string lname = foldername + leftname + "_" + rightname + "_first";
	string rname = foldername + leftname + "_" + rightname + "_second";


	Matrix3b3 r;
	Matrix3b4 RT;
	
	int im;
	MatrixXd Te2i_1;
	MatrixXd Te2i_2;
	Matrix3b3	Rxe2i_1;
	Matrix3b3	Rxe2i_2;


	im = im1;
	int sensor_id = (int)im_cam.sensor_ID(im, 0);
	int model_id = (int)im_cam.model_ID(im, 0);
	double PS=1, Cx=1, Cy=1, f=1, PSx=1, PSy=1, Cn1=1, Rn1=1;

	if (model_id == 0) {
		PS = im_cam.PS_0(sensor_id, 0);
		Cx = (-im_cam.Cn_0(sensor_id, 0) / 2 + 0.5)*PS - im_cam.xpp_0(sensor_id, 0);
		Cy = -(-im_cam.Rn_0(sensor_id, 0) / 2 + 0.5)*PS - im_cam.ypp_0(sensor_id, 0);
		f = im_cam.f_l_0(sensor_id, 0);
		PSx = PS;
		PSy = PS;
		Cn1 = im_cam.Cn_0(sensor_id, 0);
		Rn1 = im_cam.Rn_0(sensor_id, 0);
	}
	if (model_id == 4) {
		PS = im_cam.PS_4(sensor_id, 0);
		Cx = -im_cam.Cx_4(sensor_id, 0);
		Cy = im_cam.Cy_4(sensor_id, 0);
		f = im_cam.f_l_4(sensor_id, 0);
		PSx = PS;
		PSy = PS;
		Cn1 = im_cam.Cn_4(sensor_id, 0);
		Rn1 = im_cam.Rn_4(sensor_id, 0);
	}
	if (model_id == 6) {
		PS = im_cam.PS_00(sensor_id, 0);
		Cx = (-im_cam.Cn_00(sensor_id, 0) / 2 + 0.5) * PS - im_cam.xpp_00(sensor_id, 0);
		Cy = -(-im_cam.Rn_00(sensor_id, 0) / 2 + 0.5) * PS - im_cam.ypp_00(sensor_id, 0);
		f = im_cam.f_l_00(sensor_id, 0);
		PSx = PS;
		PSy = PS;
		Cn1 = im_cam.Cn_00(sensor_id, 0);
		Rn1 = im_cam.Rn_00(sensor_id, 0);
	}


	if (do_resize == true) {
		double Cn1_new = floor(Cn1*scale_factors.at(im));
		double sscale_factor = Cn1_new / Cn1;
		Cn1 = Cn1_new;
		double Rn1_new = round(Rn1*sscale_factor);
		Rn1 = Rn1_new;
		f *= sscale_factor;
		Cx *= sscale_factor;
		Cy *= sscale_factor;
		cout << Cn1 << endl;
		cout << Rn1 << endl;
		resize(ImageL, ImageL, cv::Size(Cn1, Rn1), 0, 0, INTER_CUBIC);
	}

	MatrixXd A(3, 3);
	A << PSx, 0.0, Cx,
		0.0, -PSy, Cy,
		0.0, 0.0, -f;
	MatrixXd Kmat1 = (A.inverse());

	if (model_id == 2) {
		Cn1= im_cam.Cn_2(sensor_id, 0);
		Rn1 = im_cam.Rn_2(sensor_id, 0);
		double cx_2 = im_cam.cx_2(sensor_id, 0);
		double cy_2 = im_cam.cy_2(sensor_id, 0);
		double C_2 = im_cam.C_2(sensor_id, 0);
		double f = 2 * C_2 / pi;

		if (do_resize == true) {
			double Cn1_new = floor(Cn1 * scale_factors.at(im));
			double sscale_factor = Cn1_new / Cn1;
			Cn1 = Cn1_new;
			double Rn1_new = round(Rn1 * sscale_factor);
			Rn1 = Rn1_new;
			f *= sscale_factor;
			cx_2 *= sscale_factor;
			cy_2 *= sscale_factor;
			resize(ImageL, ImageL, cv::Size(Cn1, Rn1), 0, 0, INTER_CUBIC);
		}
		Kmat1.resize(3, 3);
		
		Kmat1 << -f, 0.0, cx_2,
			0.0, f, cy_2,
			0.0, 0.0, 1.0;
	}

	if (model_id == 1) {
		Cn1=im_cam.Cn_1(sensor_id, 0);
		Rn1 = im_cam.Rn_1(sensor_id, 0);
		double fx = im_cam.fx_1(sensor_id, 0);
		double fy = im_cam.fy_1(sensor_id, 0);
		double cx = im_cam.cx_1(sensor_id, 0);
		double cy = im_cam.cy_1(sensor_id, 0);

		if (do_resize == true) {
			double Cn1_new = floor(Cn1 * scale_factors.at(im));
			double sscale_factor = Cn1_new / Cn1;
			Cn1 = Cn1_new;
			double Rn1_new = round(Rn1 * sscale_factor);
			Rn1 = Rn1_new;
			fx *= sscale_factor;
			fy *= sscale_factor;
			cx *= sscale_factor;
			cy *= sscale_factor;
			resize(ImageL, ImageL, cv::Size(Cn1, Rn1), 0, 0, INTER_CUBIC);
		}

		Kmat1.resize(3, 3);
		Kmat1 << fx, 0.0, cx,
			0.0, fy, cy,
			0.0, 0.0, 1.0;
	}

	MatrixXd Rot_x_pi(3, 3);
	Rot_x_pi << 1.0, 0, 0,
		0, -1.0, 0,
		0, 0, -1.0;


	Kmat1 = (Kmat1 / abs(Kmat1(2, 2)))*Rot_x_pi;
	r = Rot_x_pi * (im_cam.Bundestim_ori.block(im * 5 + 0, 0, 3, 3));
	Rxe2i_1 = r;
	Te2i_1 = ((im_cam.Bundestim_ori.row(im * 5 + 3)).transpose());
	RT.block(0, 0, 3, 3) = r;
	RT.block(0, 3, 3, 1) = -r * Te2i_1;
	Pmat1 = Kmat1 * RT;


	im = im2;
	sensor_id = (int)im_cam.sensor_ID(im, 0);
	model_id = (int)im_cam.model_ID(im, 0);
	double Cn2=1, Rn2=1;
	if (model_id == 0) {
		PS = im_cam.PS_0(sensor_id, 0);
		Cx = (-im_cam.Cn_0(sensor_id, 0) / 2 + 0.5)*PS - im_cam.xpp_0(sensor_id, 0);
		Cy = -(-im_cam.Rn_0(sensor_id, 0) / 2 + 0.5)*PS - im_cam.ypp_0(sensor_id, 0);
		f = im_cam.f_l_0(sensor_id, 0);
		PSx = PS;
		PSy = PS;
		Cn2 = im_cam.Cn_0(sensor_id, 0);
		Rn2 = im_cam.Rn_0(sensor_id, 0);
	}
	if (model_id == 6) {
		PS = im_cam.PS_00(sensor_id, 0);
		Cx = (-im_cam.Cn_00(sensor_id, 0) / 2 + 0.5) * PS - im_cam.xpp_00(sensor_id, 0);
		Cy = -(-im_cam.Rn_00(sensor_id, 0) / 2 + 0.5) * PS - im_cam.ypp_00(sensor_id, 0);
		f = im_cam.f_l_00(sensor_id, 0);
		PSx = PS;
		PSy = PS;
		Cn2 = im_cam.Cn_00(sensor_id, 0);
		Rn2 = im_cam.Rn_00(sensor_id, 0);
	}
	if (model_id == 4) {
		PS = im_cam.PS_4(sensor_id, 0);
		Cx = -im_cam.Cx_4(sensor_id, 0);
		Cy = im_cam.Cy_4(sensor_id, 0);
		f = im_cam.f_l_4(sensor_id, 0);
		PSx = PS;
		PSy = PS;
		Cn2 = im_cam.Cn_4(sensor_id, 0);
		Rn2 = im_cam.Rn_4(sensor_id, 0);
	}
	if (do_resize == true) {
		double Cn2_new = floor(Cn2*scale_factors.at(im));
		double sscale_factor = Cn2_new / Cn2;
		Cn2 = Cn2_new;
		double Rn2_new = round(Rn2*sscale_factor);
		Rn2 = Rn2_new;
		f *= sscale_factor;
		Cx *= sscale_factor;
		Cy *= sscale_factor;
		resize(ImageR, ImageR, cv::Size(Cn2, Rn2), 0, 0, INTER_CUBIC);
	}

	A.setZero(3, 3);
	A << PSx, 0.0, Cx,
		0.0, -PSy, Cy,
		0.0, 0.0, -f;
	MatrixXd Kmat2 = (A.inverse());

	if (model_id == 2) {
		Cn2 = im_cam.Cn_2(sensor_id, 0);
		Rn2 = im_cam.Rn_2(sensor_id, 0);
		double cx_2 = im_cam.cx_2(sensor_id, 0);
		double cy_2 = im_cam.cy_2(sensor_id, 0);
		double C_2 = im_cam.C_2(sensor_id, 0);
		double f = 2 * C_2 / pi;

		if (do_resize == true) {
			double Cn2_new = floor(Cn2 * scale_factors.at(im));
			double sscale_factor = Cn2_new / Cn2;
			Cn2 = Cn2_new;
			double Rn2_new = round(Rn2 * sscale_factor);
			Rn2 = Rn2_new;
			f *= sscale_factor;
			cx_2 *= sscale_factor;
			cy_2 *= sscale_factor;
			resize(ImageR, ImageR, cv::Size(Cn2, Rn2), 0, 0, INTER_CUBIC);
		}
		Kmat2.resize(3, 3);

		Kmat2 << -f, 0.0, cx_2,
			0.0, f, cy_2,
			0.0, 0.0, 1.0;
	}

	if (model_id == 1) {
		Cn2 = im_cam.Cn_1(sensor_id, 0);
		Rn2 = im_cam.Rn_1(sensor_id, 0);
		double fx = im_cam.fx_1(sensor_id, 0);
		double fy = im_cam.fy_1(sensor_id, 0);
		double cx = im_cam.cx_1(sensor_id, 0);
		double cy = im_cam.cy_1(sensor_id, 0);

		if (do_resize == true) {
			double Cn2_new = floor(Cn2 * scale_factors.at(im));
			double sscale_factor = Cn2_new / Cn2;
			Cn2 = Cn2_new;
			double Rn2_new = round(Rn2 * sscale_factor);
			Rn2 = Rn2_new;
			fx *= sscale_factor;
			fy *= sscale_factor;
			cx *= sscale_factor;
			cy *= sscale_factor;
			resize(ImageR, ImageR, cv::Size(Cn2, Rn2), 0, 0, INTER_CUBIC);
		}

		Kmat2.resize(3, 3);
		Kmat2 << fx, 0.0, cx,
			0.0, fy, cy,
			0.0, 0.0, 1.0;
	}

	Kmat2 = (Kmat2 / abs(Kmat2(2, 2)))*Rot_x_pi;
	r = Rot_x_pi * (im_cam.Bundestim_ori.block(im * 5 + 0, 0, 3, 3));
	Rxe2i_2 = r;
	Te2i_2 = ((im_cam.Bundestim_ori.row(im * 5 + 3)).transpose());
	RT.block(0, 0, 3, 3) = r;
	RT.block(0, 3, 3, 1) = -r * Te2i_2;
	Pmat2 = Kmat2 * RT;

	MatrixXd tmpv = Kmat1 * Rxe2i_1 * (Te2i_2 - Te2i_1);

	MatrixXd Q1 = Kmat1 * Rxe2i_1;
	MatrixXd Q2 = Kmat2 * Rxe2i_2;

	MatrixXd v1d = (Te2i_2 - Te2i_1);
	Vector3d v1(v1d(0, 0), v1d(1, 0), v1d(2, 0));
	Vector3d some = ((Rxe2i_1.row(2)).transpose());
	Vector3d v2 = some.cross(v1);
	Vector3d v3 = v1.cross(v2);

	
	R_common.setZero(3, 3);
	R_common.col(0) = v1/(v1.norm());
	R_common.col(1) = v2/(v2.norm());
	R_common.col(2) = v3/(v3.norm());
	Matrix3b3 Rtemp = (R_common.transpose());
	R_common = Rtemp;


	MatrixXd K_common = Kmat2;
	K_common(0, 1) = 0;

	MatrixXd H1 = K_common * R_common*(Q1.inverse());
	MatrixXd H2 = K_common * R_common*(Q2.inverse());

	MatrixXd bb1, bb2, bb_common;
	BoundingBox(H1, Rn1, Cn1, bb1);
	BoundingBox(H2, Rn2, Cn2, bb2);

	bb_common.setZero(4, 1);
	bb_common << floor(fmin(bb1(0, 0), bb2(0, 0))),
		floor(fmin(bb1(1, 0), bb2(1, 0))),
		floor(fmax(bb1(2, 0), bb2(2, 0))),
		floor(fmax(bb1(3, 0), bb2(3, 0)));
	bb1 = bb_common;
	bb2 = bb_common;

	double shift_x = bb1(0, 0); double shift_y = bb1(1, 0);
	MatrixXd tmpmat(3, 3); 
	tmpmat << 1, 0, -shift_x,
			0, 1, -shift_y,
			0, 0, 1;
	K_common_new1 = tmpmat * K_common;
	shift_x = bb2(0, 0); shift_y = bb2(1, 0);
	tmpmat << 1, 0, -shift_x,
		0, 1, -shift_y,
		0, 0, 1;
	K_common_new2 = tmpmat * K_common;

	int lenx = bb1(2, 0) - bb1(0, 0)+1;
	int leny = bb1(3, 0) - bb1(1, 0)+1;

	if (lenx > (2.1*Cn1) || leny > (2.1*Rn1)) {
		cout << "Unusual change of size!" << endl;
		return false;
	}
	MatrixXd MapX; MapX.setOnes(leny, lenx);
	MatrixXd MapY; MapY.setOnes(leny, lenx); 

	MatrixXd H1i = (H1.inverse());

	int ystart = bb1(1, 0); int yend = bb1(3, 0) + 1;
	int xstart = bb1(0, 0); int xend = bb1(2, 0) + 1;
	
	H1i_transform.resize(3, 4);
	H1i_transform.block(0, 0, 3, 3) = H1i;
	H1i_transform(0, 3) = xstart;
	H1i_transform(1, 3) = ystart;

	Mask.setZero(leny, lenx); Mask.fill(255);
	Concurrency::parallel_for(ystart, yend, [&](int y)
	{	
		for (int x = xstart; x < xend; x++) {
			MatrixXd Pt(3, 1);  Pt << (double)x, (double)y, 1.0;
			MatrixXd tmp = H1i * Pt;
			tmp = tmp / tmp(2,0);
			MapX(y-ystart, x-xstart) = tmp(0,0);
			MapY(y-ystart, x-xstart) = tmp(1,0);
			if (tmp(0, 0) < 50 || tmp(1, 0) < 50 || tmp(0, 0) > (Cn1- 50) || tmp(1, 0) > (Rn1 - 50)) {
				Mask(y - ystart, x - xstart) = 0;
			}
		}
	});
	
	Mat _IML(leny,lenx, ImageL.type());

	Mat X(leny, lenx, CV_32FC1);
	MatrixXf MapXf = (MapX.cast<float>());
	MapX.resize(1, 1);
	eigen2cv(MapXf, X);
	MapXf.resize(1, 1);


	Mat Y(leny, lenx, CV_32FC1);
	MatrixXf MapYf = (MapY.cast<float>());
	MapY.resize(1, 1);
	eigen2cv(MapYf, Y);
	MapYf.resize(1, 1);
	
	
	remap(ImageL, _IML, X, Y, INTER_CUBIC, 0, Scalar(0));

	imwrite(lname+ ".pgm", _IML, imwrite_pgm_params);
	_IML.release();
	ImageL.release();

	ystart = bb2(1, 0); yend = bb2(3, 0) + 1;
	xstart = bb2(0, 0); xend = bb2(2, 0) + 1;
	MapX.setOnes(leny, lenx); 
	MapY.setOnes(leny, lenx); 

	MatrixXd H2i = (H2.inverse());

	H2i_transform.resize(3, 4);
	H2i_transform.block(0, 0, 3, 3) = H2i;
	H2i_transform(0, 3) = xstart;
	H2i_transform(1, 3) = ystart;

	Concurrency::parallel_for(ystart, yend, [&](int y)
	{
		for (int x = xstart; x < xend; x++) {
			MatrixXd Pt(3, 1);  Pt << (double)x, (double)y, 1.0;
			MatrixXd tmp = H2i * Pt;
			tmp = tmp / tmp(2, 0);
			MapX(y - ystart, x - xstart) = tmp(0, 0);
			MapY(y - ystart, x - xstart) = tmp(1, 0);
		}
	});


	MapXf = (MapX.cast<float>());
	MapX.resize(1, 1);
	eigen2cv(MapXf, X);
	MapXf.resize(1, 1);

    MapYf = (MapY.cast<float>());
	MapY.resize(1, 1);
	eigen2cv(MapYf, Y);
	MapYf.resize(1, 1);

	
	Mat _IMR(leny, lenx, ImageR.type());

	remap(ImageR, _IMR, X, Y, INTER_CUBIC, 0, Scalar(0));
	
	imwrite(rname + ".pgm", _IMR, imwrite_pgm_params);

	X.release();
	Y.release();
	_IMR.release();
	ImageR.release();

	return true;
};
