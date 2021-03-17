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

double get_maximum(double a, double b) {
	double c = a;
	if (b > a) {
		c = b;
	}
	return c;
};

void Read_Mat(char *FileName, MatrixXd& m) {

	m.resize(0, 0);

	ifstream matfile;
	matfile.open(FileName, ios::in); 

	if (matfile.fail())
	{
		cout << "There was a problem reading the follwoing file: " << endl << FileName << endl;
		//exit (EXIT_FAILURE);
		return;
	}


	char* readlinechr = new char[MaxMatSize];
	vector<double> v_all;
	int nrow = 0;

	while (matfile.getline(readlinechr, MaxMatSize, '\n')) {
		nrow++;
		int stln = strlen(readlinechr);
		char* readlinestr = new char[stln + 1];
		for (int i = 0; i<stln; i++)
		{
			readlinestr[i] = readlinechr[i];
		}

		readlinestr[stln] = '\0';

		stringstream rowstream(readlinestr);
		double value;
		while (!rowstream.eof()) {
			rowstream >> value;
			v_all.push_back(value);
		}
	}
	matfile.close();


	int ncol = v_all.size() / nrow;
	m.resize(nrow, ncol);

	for (int i = 0; i<nrow; i++) {
		for (int j = 0; j<ncol; j++) {
			m(i, j) = v_all.at(i*ncol + j);
		}
	}

	return;
};

//////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////

void Write_Mat(char *FileName, MatrixXd &m, int decimal_precision) {

	ofstream matfile;
	matfile.open(FileName, ios::out); 
	if (matfile.fail()) 
	{
		cout << "There was a problem opening the follwoing file: " << endl << FileName << endl;
		//exit (EXIT_FAILURE);
		return;
	}
	matfile.flags(ios::fixed);
	matfile.precision(decimal_precision);
	matfile << m;
	matfile.close();
	return;
};

void Write_Mat(char *FileName, SparseMatrix<double> &m, int decimal_precision) {
	ofstream matfile;
	matfile.open(FileName, ios::out); 
	if (matfile.fail()) 
	{
		cout << "There was a problem opening the follwoing file: " << endl << FileName << endl;
		//exit (EXIT_FAILURE);
		return;
	}
	matfile.flags(ios::fixed);
	matfile.precision(decimal_precision);
	matfile << m;
	matfile.close();
	return;
}
//////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////


void Rotation_g2i(double Omega, double Phi, double Kappa, Matrix3b3 & Rot_g2i) {

	Matrix3b3 Mw0;
	Matrix3b3 Mf0;
	Matrix3b3 Mk0;

	Mw0 << 1, 0, 0,
		0, cos(Omega), sin(Omega),
		0, -sin(Omega), cos(Omega);

	Mf0 << cos(Phi), 0, -sin(Phi),
		0, 1, 0,
		sin(Phi), 0, cos(Phi);

	Mk0 << cos(Kappa), sin(Kappa), 0,
		-sin(Kappa), cos(Kappa), 0,
		0, 0, 1;

	Rot_g2i = Mk0 * Mf0*Mw0;
};

/////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////

void removeRow(MatrixXd& matrix, unsigned int rowToRemove) {
	unsigned int numRows = matrix.rows() - 1;
	unsigned int numCols = matrix.cols();

	if (rowToRemove < numRows) {
		matrix.block(rowToRemove, 0, numRows - rowToRemove, numCols) = matrix.block(rowToRemove + 1, 0, numRows - rowToRemove, numCols);
		matrix.conservativeResize(numRows, numCols);
	}
};

//////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////

void removeColumn(MatrixXd& matrix, unsigned int colToRemove)
{
	unsigned int numRows = matrix.rows();
	unsigned int numCols = matrix.cols() - 1;

	if (colToRemove < numCols) {
		matrix.block(0, colToRemove, numRows, numCols - colToRemove) = matrix.block(0, colToRemove + 1, numRows, numCols - colToRemove);
		matrix.conservativeResize(numRows, numCols);
	}
}

//////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////

bool does_exist(vector<double> V, double k, int & I) {
	bool res = 0;
	I = -1;
	for (int i = 0; i<V.size(); i++) {
		if (k == V.at(i)) {
			I = i;
			res = 1;
			break;
		}
	}

	return res;
};

//////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////

bool does_exist_eigen(VectorXi V, int k, int & I) {
	bool res = 0;
	I = -1;
	for (int i = 0; i<V.size(); i++) {
		if (k == V(i)) {
			I = i;
			res = 1;
			break;
		}
	}

	return res;

};
///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
void Convert_R_to_Angles(Matrix3b3 R, double& Omega, double& Phi, double& Kappa) {

	double A11 = R(0, 0);
	double A12 = R(0, 1);
	double A13 = R(0, 2);
	double A21 = R(1, 0);
	double A22 = R(1, 1);
	double A23 = R(1, 2);
	double A33 = R(2, 2);

	if (A13 != 1 && A13 != -1) {

		Phi = asin(A13);
		Omega = atan2(-A23 / cos(Phi), A33 / cos(Phi));
		Kappa = atan2(-A12 / cos(Phi), A11 / cos(Phi));

	}

	else if (A13 == 1) {
		Phi = pi / 2;
		Omega = 0; //arbitrary
		Kappa = -Omega + atan2(A21, A22);

	}

	else if (A13 == -1) {
		Phi = -pi / 2;
		Omega = 0;
		Kappa = Omega + atan2(A21, A22);
	}

	return;
};
////////////////////////////////////////////////////////
///////////////////////////////////////////////////////
void Write_ROPs(char *FileName2, MatrixXd& ROPS) {
	ofstream matfile;
	//Write EOPs only
	matfile.open(FileName2, ios::out);
	if (matfile.fail()) 
	{
		cout << "There was a problem opening the follwoing file: " << endl << FileName2 << endl;
		system("pause");
		exit(EXIT_FAILURE);
	}
	matfile.flags(ios::fixed);

	if (ROPS.cols() == 9) {
		for (int i = 0; i < ROPS.rows(); i++) {

			matfile << setprecision(0) << ROPS(i, 0) << " " << ROPS(i, 1) << " " << ROPS(i, 2) << " ";

			matfile << setprecision(8) << ROPS(i, 3) << " " << ROPS(i, 4) << " " << ROPS(i, 5) << " ";

			matfile << setprecision(15) << ROPS(i, 6) << " " << ROPS(i, 7) << " " << ROPS(i, 8) << endl;
		}
	}
	if (ROPS.cols() == 8) {
		for (int i = 0; i < ROPS.rows(); i++) {

			matfile << setprecision(0) << ROPS(i, 0) << " " << ROPS(i, 1) << " ";

			matfile << setprecision(8) << ROPS(i, 2) << " " << ROPS(i, 3) << " " << ROPS(i, 4) << " "; 

			matfile << setprecision(15) << ROPS(i, 5) << " " << ROPS(i, 6) << " " << ROPS(i, 7) << endl;
		}
	}
	matfile.close();
};

/////////////////////////////////////
/////////////////////////////////////

void Data_Preparation(char* FeatSavePP, char* DirectEO,
	CameraParam& camera_params, int& n_cams, string& imformat, 
	vector<string>& final_imlist, char* matchlist, MatrixXd& MatchList) {

	
	string path = string(FeatSavePP);
	vector<string> init_imlist;
	for (auto& p : fsys::recursive_directory_iterator(path))
	{
		if ((p.path().extension().string()).compare(imformat) == 0)
		{
			init_imlist.push_back(p.path().filename().string());
		}
	}


	char currentfile[512] = "";
	int dummyint = 0;
	currentfile[0] = '\0';
	strcat(currentfile, FeatSavePP);
	strcat(currentfile, DirectEO);
	
	ifstream matfile;
	matfile.open(currentfile, ios::in); 

	if (matfile.fail())
	{
		cout << "There was a problem reading the follwoing file: " << endl << currentfile << endl;
		//exit (EXIT_FAILURE);
		return;
	}
	string imname;
	double Xo, Yo, Zo, Om, Phi, Kap;
	int id_sensor, id_model;
	vector<string> eo_imname;
	MatrixXd eops;
	vector<int> eo_idsensor;
	vector<int> eo_idmodel;
	int thissize = 0;
	while (matfile >> imname >> Xo >> Yo >> Zo>> Om>>Phi>>Kap>>id_model>>id_sensor)
	{
		eo_imname.push_back(imname);
		eo_idsensor.push_back(id_sensor);
		eo_idmodel.push_back(id_model);
		thissize++;
		eops.conservativeResize(thissize, 6);
		eops.row(thissize-1)<< Xo, Yo, Zo, Om, Phi, Kap;
		
	}
	matfile.close();

	int imcounter = -1;
	final_imlist.clear();
	vector<int> final_eo_idsensor, final_eo_idmodel;
	MatrixXd ori;
	for (int i = 0; i < init_imlist.size(); i++) {
		bool res = false;
		int indx;
		res = findInVector(eo_imname, init_imlist.at(i), indx);
		if (res == true) {
			final_imlist.push_back(init_imlist.at(i));
			final_eo_idsensor.push_back(eo_idsensor.at(indx));
			final_eo_idmodel.push_back(eo_idmodel.at(indx));
			imcounter++;
			ori.conservativeResize(imcounter + 1, 7);
			ori.row(imcounter) << imcounter, eops.row(indx);
		}
	}
	eops.resize(1, 1);
	init_imlist.clear();
	eo_imname.clear();


	n_cams = ori.rows();
	//cout << n_cams << endl;
	(camera_params.Camera_labels).setZero(n_cams, 1);
	(camera_params.Bundestim_ori).setZero(n_cams * 5, 3);
	(camera_params.sensor_ID).setZero(n_cams, 1);
	(camera_params.model_ID).setZero(n_cams, 1);
	for (int i = 0; i < n_cams; i++) {
		camera_params.Re_init_Estimori(i, ori(i, 4)*pi / 180, ori(i, 5)*pi / 180, ori(i, 6)*pi / 180, ori(i, 1), ori(i, 2), ori(i, 3));
		camera_params.Camera_labels(i) = (int)ori(i, 0);
		// for now, we only accept traditional photo model (0) and pix4D model (4)
		int which_model= final_eo_idmodel.at(i);
		if (which_model == 0) {
			camera_params.sensor_ID(i,0) = final_eo_idsensor.at(i);
			camera_params.model_ID(i,0) = 0;
		}
		if (which_model == 1) {
			camera_params.sensor_ID(i, 0) = final_eo_idsensor.at(i);
			camera_params.model_ID(i, 0) = 1;
		}
		if (which_model == 2) {
			camera_params.sensor_ID(i, 0) = final_eo_idsensor.at(i);
			camera_params.model_ID(i, 0) = 2;
		}
		if (which_model == 3) {
			camera_params.sensor_ID(i, 0) = final_eo_idsensor.at(i);
			camera_params.model_ID(i, 0) = 3;
		}
		if (which_model == 4) {
			camera_params.sensor_ID(i, 0) = final_eo_idsensor.at(i);
			camera_params.model_ID(i, 0) =4;
		}
		if (which_model == 6) {
			camera_params.sensor_ID(i, 0) = final_eo_idsensor.at(i);
			camera_params.model_ID(i, 0) = 6;
		}
	}
	
	char currentfile2[512] = "";
	currentfile2[0] = '\0';
	strcat(currentfile2, FeatSavePP);
	strcat(currentfile2, matchlist);
	//DirectEO file has this format: Imagename X(m) Y Z Omega(deg) Phi Kappa

	matfile.open(currentfile2, ios::in);

	if (matfile.fail())
	{
		cout << "There was a problem reading the follwoing file: " << endl << currentfile2 << endl;
		return;
	}

	MatchList.setZero(n_cams, n_cams);
	string sm1, sm2;
	while (matfile >> sm1 >> sm2)
	{
		bool res1 = false; bool res2 = false;
		int indx1, indx2;
		res1 = findInVector(final_imlist, sm1, indx1);
		res2 = findInVector(final_imlist, sm2, indx2);
		if (res1 == true && res2 == true) {
			MatchList(indx1, indx2) = 1;
		}
		
	}
	matfile.close();


	ori.resize(1, 1);
	

	return;
};

//////////////////////////////////////////////
///////////////////////////////////////////////
wstring GetWC(const char *c)
{
	const size_t cSize = strlen(c) + 1;
	wstring wc(cSize, L'#');
	mbstowcs(&wc[0], c, cSize);
	
	return wc;
};
////////////////////////////////////////////////
////////////////////////////////////////////////
void BoundingBox(MatrixXd H, double Rn, double Cn, MatrixXd& bb) {

	MatrixXd corners; corners.setZero(4, 3);
	corners << 1, 1, 1,
		Cn, 1, 1,
		1, Rn, 1,
		Cn, Rn, 1;

	MatrixXd corners_x; corners_x.setZero(2, 4);

	MatrixXd tmp;
	for (int i = 0; i < 4; i++) {
		MatrixXd some = ((corners.row(i)).transpose());
		tmp = H * some;
		tmp = tmp / tmp(2, 0);
		corners_x.col(i) << tmp(0), tmp(1);
	}

	double minx = floor(corners_x.row(0).minCoeff()) ;
	double miny = floor(corners_x.row(1).minCoeff()) ;
	double maxx = ceil(corners_x.row(0).maxCoeff()) ;
	double maxy = ceil(corners_x.row(1).maxCoeff()) ;

	bb.setZero(4, 1);
	bb << minx, miny, maxx, maxy;

};
/////////////////////////////////////////////
/////////////////////////////////////////////
int sign(float x)
{
	return (x > 0) - (x < 0);
};
////////////////////////////////////////////
///////////////////////////////////////////
void BitMat(MatrixXf& comp, MatrixBools& e) {
	e.resize(comp.rows(), comp.cols());
	e.fill(false);
	for (int i = 0; i < comp.rows(); i++) {
		for (int j = 0; j < comp.cols(); j++) {
			if (sign(comp(i, j)) == -1)
				e(i, j) = true;
		}
	}
};
//////////////////////////////////////////
//////////////////////////////////////////
int HamDist(MatrixBools& I1, MatrixBools& I2) {
	int dist = 0;
	for (int i = 0; i < I1.rows(); i++) {
		for (int j = 0; j < I1.cols(); j++) {
			bool d = I1(i, j) ^ I2(i, j);
			if (d == true) {
				dist++;
			}
		}
	}
	return dist;
};
/////////////////////////////////////
///////////////////////////////////
void pnm_read(std::ifstream &file, char *buf) {
	char doc[BUF_SIZE];
	char c;

	file >> c;
	while (c == '#') {
		file.getline(doc, BUF_SIZE);
		file >> c;
	}
	file.putback(c);

	file.width(BUF_SIZE);
	file >> buf;
	file.ignore();
};

image<uchar> *loadPGM(const char *name) {
	char buf[BUF_SIZE];

	// read header
	std::ifstream file(name, std::ios::in | std::ios::binary);
	pnm_read(file, buf);
	if (strncmp(buf, "P5", 2)) {
		std::cout << "ERROR: Could not read file " << name << std::endl;
		throw pnm_error();
	}

	pnm_read(file, buf);
	int width = atoi(buf);
	pnm_read(file, buf);
	int height = atoi(buf);

	pnm_read(file, buf);
	if (atoi(buf) > UCHAR_MAX) {
		std::cout << "ERROR: Could not read file " << name << std::endl;
		throw pnm_error();
	}

	// read data
	image<uchar> *im = new image<uchar>(width, height);
	file.read((char *)imPtr(im, 0, 0), width * height * sizeof(uchar));

	return im;
};

void savePGM(image<uchar> *im, const char *name) {
	int width = im->width();
	int height = im->height();
	std::ofstream file(name, std::ios::out | std::ios::binary);

	file << "P5\n" << width << " " << height << "\n" << UCHAR_MAX << "\n";
	file.write((char *)imPtr(im, 0, 0), width * height * sizeof(uchar));
};
/////////////////////////////////////
/////////////////////////////////////
void ImageGeometry_Remove(char* FeatSavePP, string imname, string imformat) {

	string foldername = string(FeatSavePP);
	int found = imname.find(imformat);
	string newname = foldername + "Undist_" + (imname.substr(0, found));
	newname += ".tiff";
	char* namecstr = new char[newname.length() + 1];
	strcpy(namecstr, newname.c_str());

	int b = remove(namecstr);
}
//////////////////////////////////////
//////////////////////////////////////
void get_K_R_T_P(CameraParam& this_camera, int im, bool do_resize, double Cn1_new, double& Rn1, MatrixXd& K, MatrixXd& R, MatrixXd& T, Matrix3b4& P) {

	int sensor_id = (int)this_camera.sensor_ID(im, 0);
	int model_id = (int)this_camera.model_ID(im, 0);
	double PS=1, Cx=1, Cy=1, f=1, PSx=1, PSy=1, Cn1=1;
	Rn1 = 1;

	if (model_id == 0) {
		PS = this_camera.PS_0(sensor_id, 0);
		Cx = (-this_camera.Cn_0(sensor_id, 0) / 2 + 0.5) * PS - this_camera.xpp_0(sensor_id, 0);
		Cy = -(-this_camera.Rn_0(sensor_id, 0) / 2 + 0.5) * PS - this_camera.ypp_0(sensor_id, 0);
		f = this_camera.f_l_0(sensor_id, 0);
		PSx = PS;
		PSy = PS;
		Cn1 = this_camera.Cn_0(sensor_id, 0);
		Rn1 = this_camera.Rn_0(sensor_id, 0);
	}
	if (model_id == 4) {
		PS = this_camera.PS_4(sensor_id, 0);
		Cx = -this_camera.Cx_4(sensor_id, 0);
		Cy = this_camera.Cy_4(sensor_id, 0);
		f = this_camera.f_l_4(sensor_id, 0);
		PSx = PS;
		PSy = PS;
		Cn1 = this_camera.Cn_4(sensor_id, 0);
		Rn1 = this_camera.Rn_4(sensor_id, 0);
	}
	if (model_id == 6) {
		PS = this_camera.PS_00(sensor_id, 0);
		Cx = (-this_camera.Cn_00(sensor_id, 0) / 2 + 0.5) * PS - this_camera.xpp_00(sensor_id, 0);
		Cy = -(-this_camera.Rn_00(sensor_id, 0) / 2 + 0.5) * PS - this_camera.ypp_00(sensor_id, 0);
		f = this_camera.f_l_00(sensor_id, 0);
		PSx = PS;
		PSy = PS;
		Cn1 = this_camera.Cn_00(sensor_id, 0);
		Rn1 = this_camera.Rn_00(sensor_id, 0);
	}

	if (do_resize == true) {
		double sscale_factor = Cn1_new / Cn1;
		Cn1 = Cn1_new;
		double Rn1_new = round(Rn1 * sscale_factor);
		Rn1 = Rn1_new;
		f *= sscale_factor;
		Cx *= sscale_factor;
		Cy *= sscale_factor;
	}

	MatrixXd A(3, 3);
	A << PSx, 0.0, Cx,
		0.0, -PSy, Cy,
		0.0, 0.0, -f;
	K = (A.inverse());
	K = (K / abs(K(2, 2)));

	if (model_id == 2) {
		Cn1 = this_camera.Cn_2(sensor_id, 0);
		Rn1 = this_camera.Rn_2(sensor_id, 0);
		double cx_2 = this_camera.cx_2(sensor_id, 0);
		double cy_2 = this_camera.cy_2(sensor_id, 0);
		double C_2 = this_camera.C_2(sensor_id, 0);
		double f = 2 * C_2 / pi;

		if (do_resize == true) {
			double sscale_factor = Cn1_new / Cn1;
			Cn1 = Cn1_new;
			double Rn1_new = round(Rn1 * sscale_factor);
			Rn1 = Rn1_new;
			f *= sscale_factor;
			cx_2 *= sscale_factor;
			cy_2 *= sscale_factor;
		}
		K.resize(3, 3);

		K << -f, 0.0, cx_2,
			0.0, f, cy_2,
			0.0, 0.0, 1.0;
	}

	if (model_id == 1) {
		Cn1 = this_camera.Cn_1(sensor_id, 0);
		Rn1 = this_camera.Rn_1(sensor_id, 0);
		double fx = this_camera.fx_1(sensor_id, 0);
		double fy = this_camera.fy_1(sensor_id, 0);
		double cx = this_camera.cx_1(sensor_id, 0);
		double cy = this_camera.cy_1(sensor_id, 0);

		if (do_resize == true) {
			double sscale_factor = Cn1_new / Cn1;
			Cn1 = Cn1_new;
			double Rn1_new = round(Rn1 * sscale_factor);
			Rn1 = Rn1_new;
			fx *= sscale_factor;
			fy *= sscale_factor;
			cx *= sscale_factor;
			cy *= sscale_factor;
		}

		K.resize(3, 3);
		K << fx, 0.0, cx,
			0.0, fy, cy,
			0.0, 0.0, 1.0;
	}

	R = (this_camera.Bundestim_ori.block(im * 5 + 0, 0, 3, 3));
	T = ((this_camera.Bundestim_ori.row(im * 5 + 3)).transpose());
	Matrix3b4 RT;
	RT.block(0, 0, 3, 3) = R;
	RT.block(0, 3, 3, 1) = -R * T;
	P = K * RT;
};
//////////////////////////////////
/////////////////////////////////
double Det_FormTo4(MatrixXd A, MatrixXd B, MatrixXd C, MatrixXd D) {
	MatrixXd E(4, 4);
	E.col(0) = A;
	E.col(1) = B;
	E.col(2) = C;
	E.col(3) = D;
	return E.determinant();
};
///////////////////////////////////
//////////////////////////////////
void TriFocal_Overlay(string prefix, string postfix, string imformat, string foldername, vector<string> ImageList, CameraParam& this_camera,
	MatrixXd K1, MatrixXd K2, MatrixXd T1, MatrixXd T2, MatrixXd R1, MatrixXd R2,
	Matrix3b4 P1, Matrix3b4 P2,
	MatrixXd A1, MatrixXd A2, MatrixXd B1, MatrixXd B2, MatrixXd C1, MatrixXd C2,
	MatrixXd& DLeft, MatrixXi& Mask,
	MatrixXf& Overlay_1x, MatrixXf& Overlay_1y, MatrixXf& Overlay_2x, MatrixXf& Overlay_2y)
{
	string im3_name = prefix + postfix + imformat;
	bool res3 = false;
	int indx3;
	res3 = findInVector(ImageList, im3_name, indx3);
	int im = indx3;

	MatrixXd K3, R3, T3, A3, B3, C3;
	Matrix3b4 P3;
	bool do_resize = false;
	double Cn1_new = 1.0;
	double Rn1_new = 1.0;
	/*if (postfix.compare("RGB") == 0) {
		do_resize = true;
	    Cn1_new = 1280.0;
	}*/
	get_K_R_T_P(this_camera, im, do_resize, Cn1_new, Rn1_new, K3, R3, T3, P3);
	A3 = (P3.row(0).transpose()); B3 = (P3.row(1).transpose()); C3 = (P3.row(2).transpose());

	Matrix3b3 Trif_tens1, Trif_tens2, Trif_tens3;
	Trif_tens1 << Det_FormTo4(B1, C1, A2, A3), Det_FormTo4(B1, C1, A2, B3), Det_FormTo4(B1, C1, A2, C3),
		Det_FormTo4(B1, C1, B2, A3), Det_FormTo4(B1, C1, B2, B3), Det_FormTo4(B1, C1, B2, C3),
		Det_FormTo4(B1, C1, C2, A3), Det_FormTo4(B1, C1, C2, B3), Det_FormTo4(B1, C1, C2, C3);
	Trif_tens1 = -Trif_tens1;

	Trif_tens2 << Det_FormTo4(C1, A1, A2, A3), Det_FormTo4(C1, A1, A2, B3), Det_FormTo4(C1, A1, A2, C3),
		Det_FormTo4(C1, A1, B2, A3), Det_FormTo4(C1, A1, B2, B3), Det_FormTo4(C1, A1, B2, C3),
		Det_FormTo4(C1, A1, C2, A3), Det_FormTo4(C1, A1, C2, B3), Det_FormTo4(C1, A1, C2, C3);
	Trif_tens2 = -Trif_tens2;

	Trif_tens3 << Det_FormTo4(A1, B1, A2, A3), Det_FormTo4(A1, B1, A2, B3), Det_FormTo4(A1, B1, A2, C3),
		Det_FormTo4(A1, B1, B2, A3), Det_FormTo4(A1, B1, B2, B3), Det_FormTo4(A1, B1, B2, C3),
		Det_FormTo4(A1, B1, C2, A3), Det_FormTo4(A1, B1, C2, B3), Det_FormTo4(A1, B1, C2, C3);
	Trif_tens3 = -Trif_tens3;

	MatrixXd B32 = T3 - T2;
	MatrixXd T21 = R2 * B32;
	MatrixXd T32x(3, 3);
	T32x << 0, -T21(2, 0), T21(1, 0),
		T21(2, 0), 0, -T21(0, 0),
		-T21(1, 0), T21(0, 0), 0;
	MatrixXd R3t = (R3.transpose());
	MatrixXd Emat = T32x * (R2 * R3t);
	MatrixXd K2inv = (K2.inverse());
	MatrixXd K2inv_t = (K2inv.transpose());
	MatrixXd K3inv = (K3.inverse());
	MatrixXd Fmat_2_3 = K2inv_t * Emat * K3inv;
	double f1, f2, f3, f4, f5, f6, f7, f8, f9;
	f1 = Fmat_2_3(0, 0); f2 = Fmat_2_3(0, 1); f3 = Fmat_2_3(0, 2);
	f4 = Fmat_2_3(1, 0); f5 = Fmat_2_3(1, 1); f6 = Fmat_2_3(1, 2);
	f7 = Fmat_2_3(2, 0); f8 = Fmat_2_3(2, 1); f9 = Fmat_2_3(2, 2);

	MatrixXf Overlay_3x(DLeft.rows(), DLeft.cols());
	MatrixXf Overlay_3y(DLeft.rows(), DLeft.cols());

	for (int i = 0; i < DLeft.rows(); i++) {
		for (int j = 0; j < DLeft.cols(); j++) {
			double disp = DLeft(i, j);
			if (disp != -10 && Mask(i, j) == 255) {
				double x1, x2, y1, y2;
				x1 = Overlay_1x(i, j);
				x2 = Overlay_2x(i, j);
				y1 = Overlay_1y(i, j);
				y2 = Overlay_2y(i, j);

				MatrixXd Tbar = x1 * Trif_tens1 + y1 * Trif_tens2 + Trif_tens3;
				double Tbar_n = 1 / (Tbar.norm());
				Tbar = Tbar_n * Tbar;

				double t1, t2, t3, t4, t5, t6, t7, t8, t9;
				t1 = Tbar(0, 0); t2 = Tbar(0, 1); t3 = Tbar(0, 2);
				t4 = Tbar(1, 0); t5 = Tbar(1, 1); t6 = Tbar(1, 2);
				t7 = Tbar(2, 0); t8 = Tbar(2, 1); t9 = Tbar(2, 2);

				double y3 = y2 + 3;
				double x3 = -(y3 * (f8 + f2 * x2 + f5 * y2) + f3 * x2 + f6 * y2 + f9) / (f7 + f1 * x2 + f4 * y2);

				int max_counter = 4;
				if (postfix.compare("RGB") == 0) {
					max_counter = 10;
				}
				for (int counter = 0; counter <= max_counter; counter++)
				{
					MatrixXd dL(4, 1);
					dL << t8 * y2 - t5 + y3 * (t6 - t9 * y2),
						t4 - t7 * y2 + ((t6 - t9 * y2) * (f9 + y3 * (f8 + f2 * x2 + f5 * y2) + f3 * x2 + f6 * y2)) / (f7 + f1 * x2 + f4 * y2),
						t2 - t8 * x2 - y3 * (t3 - t9 * x2),
						t7* x2 - t1 - ((t3 - t9 * x2) * (f9 + y3 * (f8 + f2 * x2 + f5 * y2) + f3 * x2 + f6 * y2)) / (f7 + f1 * x2 + f4 * y2);

					MatrixXd A(4, 1);
					A << t6 - t9 * y2,
						((t6 - t9 * y2) * (f8 + f2 * x2 + f5 * y2)) / (f7 + f1 * x2 + f4 * y2),
						t9* x2 - t3,
						-((t3 - t9 * x2) * (f8 + f2 * x2 + f5 * y2)) / (f7 + f1 * x2 + f4 * y2);

					MatrixXd some = (A.transpose());
					MatrixXd bb = some * A;
					MatrixXd ll = some * dL;
					double dX = -1.0 / bb(0, 0) * ll(0, 0);

					y3 = y3 + dX;
					x3 = -(y3 * (f8 + f2 * x2 + f5 * y2) + f3 * x2 + f6 * y2 + f9) / (f7 + f1 * x2 + f4 * y2);

				}
				Overlay_3x(i, j) = x3;
				Overlay_3y(i, j) = y3;
			}
		}
	}

	vector<int> imwrite_jpg_params;
	imwrite_jpg_params.push_back(IMWRITE_JPEG_QUALITY);
	imwrite_jpg_params.push_back(100);

	string imfullname = foldername + "Undist_" + prefix + postfix + ".tiff";
	string newname = foldername + "Lay_" + prefix + postfix + ".tiff";
	Mat Image3 = imread(imfullname, IMREAD_COLOR | IMREAD_ANYDEPTH | IMREAD_IGNORE_ORIENTATION);
	/*if (postfix.compare("RGB") == 0) {
		resize(Image3, Image3, cv::Size(Cn1_new, Rn1_new), 0, 0, INTER_CUBIC);
	}*/
	//imshow("Display Window", Image3);
	//waitKey(0);
	Mat Imcorrect3(Image3.size(), Image3.type());
	Mat X3(Image3.size(), CV_32FC1);
	Mat Y3(Image3.size(), CV_32FC1);
	eigen2cv(Overlay_3x, X3);
	eigen2cv(Overlay_3y, Y3);
	remap(Image3, Imcorrect3, X3, Y3, INTER_CUBIC, 0, Scalar(255, 255, 255));
	X3.release();
	Y3.release();
	Image3.release();

	imwrite(newname, Imcorrect3);//, imwrite_jpg_params
	Imcorrect3.release();
};
