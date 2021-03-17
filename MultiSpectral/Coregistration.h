#ifndef _COREGISTRATION_H_
#define _COREGISTRATION_H_


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
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/Geometry>
#include <fstream> 
#include <stdlib.h>
#include <vector>
#include <math.h>  
#include <algorithm>
#include <stdexcept>
#include <iomanip> 
#include <ppl.h>
#include "lapacke.h"

#include "opencv2/core/core.hpp"
#include "opencv2/imgcodecs/imgcodecs.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/core/eigen.hpp"
#include "opencv2/opencv.hpp" //for pm process
#include "opencv2/photo.hpp"

#include <stdint.h>
#include <climits>
#include <cstring>
#include <stdio.h>
#include <stdint.h>
#include <emmintrin.h>
#include <pmmintrin.h>
#include <cassert>
#include "elas.h"
#include "descriptor.h"
#include "triangle.h"
#include "matrix.h"
#include "image.h"
#include "filter.h"
#include "cxxopts.hpp"
#include <charconv>

#include <filesystem>
namespace fsys = std::filesystem;



using namespace cv;
using namespace std;
using namespace Eigen;

#define MaxMatSize 10000000
#define pi 3.14159265358979323846 
#define MAX_SHORT 1.1E25
#define EIGEN_USE_BLAS
#define EIGEN_USE_LAPACKE 
#define EIGEN_USE_LAPACKE_STRICT	
#define MYNAN -666666
#define PMYNAN 666666

typedef Matrix<double, 2, 1> Matrix2b1;
typedef Matrix<double, Dynamic, 3> Matrixdb3;
typedef Matrix<double, 3, 3> Matrix3b3;
typedef Matrix<double, 3, 4> Matrix3b4;
typedef Matrix<double, Dynamic, 2> Matrixdby2;
typedef Matrix<double, Dynamic, 3> Matrixdby3;
typedef Matrix<int, Dynamic, 2> Matrixdby2i;
typedef Matrix<bool, Dynamic, Dynamic> MatrixBools;



//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`

void Read_Mat(char *FileName, MatrixXd& m);
void Write_Mat(char *FileName, MatrixXd &m, int decimal_precision); 
void Write_Mat(char *FileName, SparseMatrix<double> &m, int decimal_precision); 
void Write_ROPs(char *FileName2, MatrixXd& ROPS);
int sign(float x); 
void BitMat(MatrixXf& comp, MatrixBools& e);
int HamDist(MatrixBools& I1, MatrixBools& I2);

void Rotation_g2i(double Omega, double Phi, double Kappa, Matrix3b3 & Rot_g2i);

void Convert_R_to_Angles(Matrix3b3 R, double& Omega, double& Phi, double& Kappa);

bool does_exist(vector<double> V, double k, int & I);
bool does_exist_eigen(VectorXi V, int k, int & I);


void removeRow(MatrixXd& matrix, unsigned int rowToRemove);
void removeColumn(MatrixXd& matrix, unsigned int colToRemove);
void ImageGeometry_Remove(char* FeatSavePP, string imname, string imformat);

wstring GetWC(const char *c);

double get_maximum(double a, double b);

template <typename MatrixT>
bool compareRows(MatrixT a, MatrixT b) {
	return a(0, 0) < b(0, 0);
}

template <typename Scalar, int rows, int cols, int options, int maxRows, int maxCols>
void  SortRows(Matrix<Scalar, rows, cols, options, maxRows, maxCols> & sorted, Matrix<Scalar, rows, cols, options, maxRows, maxCols>& target, int coltofocus) {

	Matrix<Scalar, rows, cols, options, maxRows, maxCols> target_temp;
	target_temp = target;

	target_temp.col(0) = target.col(coltofocus);
	target_temp.col(coltofocus) = target.col(0);

	vector<Eigen::Matrix<Scalar, 1, cols>> matrixRows;
	for (unsigned int i = 0; i < target_temp.rows(); i++) {
		matrixRows.push_back(target_temp.row(i));
	}


	sort(matrixRows.begin(),
		matrixRows.end(),
		compareRows<Matrix<Scalar, 1, cols>>);


	Matrix<Scalar, rows, cols, options, maxRows, maxCols> sorted_temp;
	sorted_temp = target;

	for (unsigned int i = 0; i < matrixRows.size(); i++) {
		sorted_temp.row(i) << matrixRows[i];
	}

	sorted = sorted_temp;
	sorted.col(0) = sorted_temp.col(coltofocus);
	sorted.col(coltofocus) = sorted_temp.col(0);

	return;
};
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
template < typename T>
bool findInVector(const vector<T>& vecOfElements, const T& element, int& index)
{
	bool result;

	auto it = find(vecOfElements.begin(), vecOfElements.end(), element);

	if (it != vecOfElements.end())
	{
		index = distance(vecOfElements.begin(), it);
		result = true;
	}
	else
	{
		result = false;
		index = -1;
	}

	return result;
};
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
class CameraParam
{
private:
	//none
public:

	// calibration (IO) parameters for model 0 (traditional photogrammetry)
	MatrixXd  xpp_0;
	MatrixXd  ypp_0;
	MatrixXd  f_l_0;
	MatrixXd  K1_0;
	MatrixXd  K2_0;
	MatrixXd  K3_0;
	MatrixXd  K4_0;
	MatrixXd  K5_0;
	MatrixXd  P1_0;
	MatrixXd  P2_0;
	MatrixXd  S1_0;
	MatrixXd  S2_0;
	MatrixXd  PS_0;
	MatrixXd  Cn_0;
	MatrixXd  Rn_0;
	MatrixXd  Kcoefs_0;

	// calibration (IO) parameters for model 6 (specific to Z-camera)
	MatrixXd  xpp_00; //principal point x-offset (mm)
	MatrixXd  ypp_00; //principal point y-offset (mm)
	MatrixXd  f_l_00; //principal distance (focal length) (mm)
	MatrixXd  K1_00; //radial lens distortrion coeff
	MatrixXd  K2_00;
	MatrixXd  K3_00;
	MatrixXd  K4_00;
	MatrixXd  K5_00;
	MatrixXd  P1_00;//tangential lens distortion coeff
	MatrixXd  P2_00;
	MatrixXd  SS1_00;//
	MatrixXd  SS2_00;//
	MatrixXd  SS3_00;//
	MatrixXd  SS4_00;//
	MatrixXd  SS5_00;//
	MatrixXd  PS_00; //pixel size  (mm)
	MatrixXd  Cn_00; //number of columns of the image (pixels)
	MatrixXd  Rn_00; //number of rows of the image (pixels)
	MatrixXd  Kcoefs_00; //number of radial lens distortions to be calibrated or considered


	// calibration (IO) parameters for model 1 (open cv model)
	MatrixXd  Cn_1;
	MatrixXd  Rn_1;
	MatrixXd  cx_1;
	MatrixXd  cy_1;
	MatrixXd  fx_1;
	MatrixXd  fy_1;
	MatrixXd  K1_1;
	MatrixXd  K2_1;
	MatrixXd  K3_1;
	MatrixXd  P1_1;
	MatrixXd  P2_1;


	// calibration (IO) parameters for model 2 (fisheye- pix4d)
	MatrixXd  Cn_2;
	MatrixXd  Rn_2;
	MatrixXd  cx_2;
	MatrixXd  cy_2;
	MatrixXd  C_2;
	MatrixXd  D_2;
	MatrixXd  E_2;
	MatrixXd  F_2;
	MatrixXd  P2_2;
	MatrixXd  P3_2;
	MatrixXd  P4_2;


	// calibration (IO) parameters for model 3 (fisheye- equi-distant - traditional)
	MatrixXd  xpp_3;
	MatrixXd  ypp_3;
	MatrixXd  f_l_3;
	MatrixXd  K1_3;
	MatrixXd  K2_3;
	MatrixXd  K3_3;
	MatrixXd  K4_3;
	MatrixXd  K5_3;
	MatrixXd  P1_3;
	MatrixXd  P2_3;
	MatrixXd  S1_3;
	MatrixXd  S2_3;
	MatrixXd  PS_3;
	MatrixXd  Cn_3;
	MatrixXd  Rn_3;
	MatrixXd  Kcoefs_3;


	// calibration (IO) parameters for model 4 (pix4D model)
	MatrixXd  Cx_4;
	MatrixXd  Cy_4;
	MatrixXd  f_l_4;
	MatrixXd  K1_4;
	MatrixXd  K2_4;
	MatrixXd  K3_4;
	MatrixXd  P1_4;
	MatrixXd  P2_4;
	MatrixXd  PS_4;
	MatrixXd  Cn_4;
	MatrixXd  Rn_4;


	VectorXi Camera_labels;
	MatrixBools Calibrated_0, Calibrated_1, Calibrated_2, Calibrated_3, Calibrated_4, Calibrated_6;

	int num_sensors_model0; //number of sensors that follow modelid 0
	int num_sensors_model1; //number of sensors that follow modelid 1
	int num_sensors_model2; //number of sensors that follow modelid 2
	int num_sensors_model3; //number of sensors that follow modelid 3
	int num_sensors_model4; //number of sensors that follow modelid 4
	int num_sensors_model6; //number of sensors that follow modelid 6


	MatrixXd sensor_ID;
	MatrixXd model_ID;
	Matrixdb3 Bundestim_ori;

	CameraParam(int nsensorwithmodelid0, int nsensorwithmodelid1,
		int nsensorwithmodelid2, int nsensorwithmodelid3, int nsensorwithmodelid4, int nsensorwithmodelid6);

	~CameraParam();

	void Re_init_Estimori(int i, double Omega, double Phi, double Kappa, double Xo, double Yo, double Zo);

	void ResetIOPs0(int i, bool calibrated, double PSv, double Cnv, double Rnv, double K1v, double K2v, double K3v, double K4v, double K5v, double P1v, double P2v, double S1v, double S2v, double xppv, double yppv, double f_lv, double kcoeffsv);
	void ResetIOPs1(int i, bool calibrated, double Cnv, double Rnv, double K1v, double K2v, double K3v, double P1v, double P2v, double cxv, double cyv, double fxv, double fyv);
	void ResetIOPs4(int i, bool calibrated, double PSv, double Cnv, double Rnv, double K1v, double K2v, double K3v, double P1v, double P2v, double cxv, double cvv, double f_lv);
	void ResetIOPs2(int i, bool calibrated, double Cnv, double Rnv, double P2v, double P3v, double P4v, double cxv, double cyv, double Cv, double Dv, double Fv);
	void ResetIOPs3(int i, bool calibrated, double PSv, double Cnv, double Rnv, double K1v, double K2v, double K3v, double K4v, double K5v, double P1v, double P2v, double S1v, double S2v, double xppv, double yppv, double f_lv, double kcoeffsv);
	void ResetIOPs6(int i, bool calibrated, double PSv, double Cnv, double Rnv, double K1v, double K2v, double K3v, double K4v, double K5v, double P1v, double P2v, double SS1v, double SS2v, double SS3v, double SS4v, double SS5v, double xppv, double yppv, double f_lv, double kcoeffsv);


	void ResetEOPs(Matrixdb3& Estimori);
	void WriteCameraToFile(char* FileName, char* FileName2);

};

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void Data_Preparation(char* FeatSavePP, char* DirectEO,
	CameraParam& camera_params, int& n_cams, string& imformat,
	vector<string>& final_imlist, char* matchlist, MatrixXd& MatchList);

void BoundingBox(MatrixXd H, double Rn, double Cn, MatrixXd& bb);
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class ImageGeometry
{
private:
	//nothing
public:
	int Cn; 
	int Rn;

	MatrixXf UndistortionX; 
	MatrixXf UndistortionY;

	ImageGeometry(int Cn=1, int Rn=1); 
	~ImageGeometry(); 

	void Undistortion_mat(CameraParam& im_cam, double model_id, double sensor_id, char* FeatSavePP);
	void Undistort(char* FeatSavePP, string imname, CameraParam& im_cam, int im, string imformat);

};
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
struct PCpoint {
	double x;
	double y;
	double z;
	int r;
	int g;
	int b;
	int im1;
	int im2;
};
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
size_t WriteByBuffer(vector<char>& buffer, size_t offset, double number, char separator);
void WriteCSV(const string& csv_filepath, const vector<PCpoint>& points);
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
bool Rectifier(char* FeatSavePP, string ImageNameL, string ImageNameR, int im1, int im2,
	CameraParam& im_cam,
	int process_sizex, int process_sizey, Matrix3b3& R_common,
	MatrixXd& K_common_new1, MatrixXd& K_common_new2, vector<double> scale_factor, bool do_resize, bool process_in_batch,
	Matrix<int, Dynamic, Dynamic>& Mask,
	string imformat,
	MatrixXd& H1i_transform, MatrixXd& H2i_transform,
	Matrix3b4& Pmat1, Matrix3b4& Pmat2,
	bool is_coregistration);
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void Elas_process(MatrixXd& DLeft, const char* file_1, const char* file_2, Matrix<int, Dynamic, Dynamic>& Mask);

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void get_K_R_T_P(CameraParam& this_camera, int im, bool do_resize, double Cn1, double& Rn1, MatrixXd& K, MatrixXd& R, MatrixXd& T, Matrix3b4& P);
double Det_FormTo4(MatrixXd A, MatrixXd B, MatrixXd C, MatrixXd D);
void TriFocal_Overlay(string prefix, string postfix, string imformat, string foldername, vector<string> ImageList, CameraParam& this_camera,
	MatrixXd K1, MatrixXd K2, MatrixXd T1, MatrixXd T2, MatrixXd R1, MatrixXd R2,
	Matrix3b4 P1, Matrix3b4 P2,
	MatrixXd A1, MatrixXd A2, MatrixXd B1, MatrixXd B2, MatrixXd C1, MatrixXd C2,
	MatrixXd& DLeft, MatrixXi& Mask,
	MatrixXf& Overlay_1x, MatrixXf& Overlay_1y, MatrixXf& Overlay_2x, MatrixXf& Overlay_2y);
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void PointCloud(char* FeatSavePP, string ImageNameL, int im1, int im2, CameraParam& im_cam,
	MatrixXd& DLeft, string SaveName,
	MatrixXd& delta_trans, Matrix3b3& R_common, MatrixXd& K_common1, MatrixXd& K_common2,
	Matrix<int, Dynamic, Dynamic>& Mask, double maxdist_im2cloud,
	vector<PCpoint>& pcvec, int& pc_counter, string imformat);
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

struct path {
	short rowDiff;
	short colDiff;
};

void PatchMatch_process(MatrixXd& DLeft, string file_1, string file_2); //for pm process
#endif