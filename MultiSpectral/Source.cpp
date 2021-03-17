#include "Coregistration.h"

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

int main(int argc, char** argv) {
	// Note: The program's current version has no specifi safe gaurds about bad input arguments including bad files/directories!!!
	int sens_model0, sens_model1, sens_model2, sens_model3, sens_model4, sens_model6;
	string imformat = ".tiff";
	string SaveName = "PC1";
	bool Keep_undistored = true;
	int modelid = 0;
	int process_sizex = 150;
	int process_sizey = 300;
	bool do_resize = false;
	double scale_factor = 1;
	double size_threshold = 1500;
	double maxdist_im2cloud = 60;
	bool is_coregistration = true;
	bool process_in_batch = false;
	char DirectEO[512] = "AllEOPs.txt"; 
	char TransFile[512] = "AllTrans.txt";
	char MatchFile[512] = "AllMatches.txt";
	
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	cxxopts::Options options("CLI", "DenseMatching");

	options.add_options()
		("f,folder", "Input the data directory", cxxopts::value<std::string>(), "[filepath]")
		("a,model0", "Calibration model", cxxopts::value<int>()->default_value(std::to_string(0)), "[0]")
		("b,model1", "Calibration model", cxxopts::value<int>()->default_value(std::to_string(0)), "[0]")
		("c,model2", "Calibration model", cxxopts::value<int>()->default_value(std::to_string(0)), "[0]")
		("d,model3", "Calibration model", cxxopts::value<int>()->default_value(std::to_string(0)), "[0]")
		("e,model4", "Calibration model", cxxopts::value<int>()->default_value(std::to_string(0)), "[0]")
		("g,model6", "Calibration model", cxxopts::value<int>()->default_value(std::to_string(0)), "[0]")
		("i,imformat", "Input the format of your original images", cxxopts::value<std::string>()->default_value("tiff"), "[format]")
		("n,name", "Input the name of the output point cloud", cxxopts::value<std::string>()->default_value("PointCloud"), "[name]")
		("s,size", "Input the maximum rectification width (must be below 1500)", cxxopts::value<double>()->default_value(std::to_string(1700)), "[pixels]")
		("k,keep", "Keep the undistorted images", cxxopts::value<int>()->default_value(std::to_string(1)), "[1]")
		("x,maxdist","Maximum distance of cloud from images", cxxopts::value<double>()->default_value(std::to_string(60)), "[meters]")
		("q,coreg", "Coregistration", cxxopts::value<int>()->default_value(std::to_string(1)), "[1]")
		;

	
	auto parsed_options = options.parse(argc, argv);
	Keep_undistored = (bool)parsed_options["keep"].as<int>();
	size_threshold = parsed_options["size"].as<double>();
	maxdist_im2cloud= parsed_options["maxdist"].as<double>();
	SaveName = parsed_options["name"].as<std::string>();
	imformat = "." + parsed_options["imformat"].as<std::string>();
	string folder= parsed_options["folder"].as<std::string>();
	char * ImageFolder = new char[folder.length() + 1];
	strcpy(ImageFolder, folder.c_str());
	
	sens_model0 = parsed_options["model0"].as<int>();
	sens_model1 = parsed_options["model1"].as<int>();
	sens_model2 = parsed_options["model2"].as<int>();
	sens_model3 = parsed_options["model3"].as<int>();
	sens_model4 = parsed_options["model4"].as<int>();
	sens_model6 = parsed_options["model6"].as<int>();

	is_coregistration = (bool) parsed_options["coreg"].as<int>();

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	MatrixXd delta_trans;
	char currentfile[512] = "";
	currentfile[0] = '\0';
	strcat(currentfile, ImageFolder);
	strcat(currentfile, TransFile);
	Read_Mat(currentfile, delta_trans);

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	//these are the IOPs and sensor parameters of the camera with which the images were taken
	CameraParam this_camera(sens_model0, sens_model1, sens_model2, sens_model3, sens_model4, sens_model6);

	char answer;
	bool calibrated;
	for (int i = 0; i <= sens_model0 - 1; i++) {
	
		//Format of this file:
		//pixel pitch (mm)
		//width (pix)
		//height (pix)
		//K1
		//K2
		//K3
		//K4 (optional)
		//K5 (optional)
		//P1
		//P2
		//S1
		//S2
		//xpp
		//ypp
		//f (mm)
		
		string camfile1 = "Cam_M0_"+to_string(i)+".cam";
		char* camerapamfile = new char[camfile1.length() + 1];
		strcpy(camerapamfile, camfile1.c_str());

		MatrixXd campars;
		char currentfile[512] = "";
		currentfile[0] = '\0';
		strcat(currentfile, ImageFolder);
		strcat(currentfile, camerapamfile);
		Read_Mat(currentfile, campars);

		calibrated = true;

		if (campars.rows() == 10) { //No K coefficient
			double kcoefs = 0;
			this_camera.ResetIOPs0(i, calibrated, campars(0, 0), campars(1, 0), campars(2, 0), 0.0, 0.0, 0.0, 0.0, 0.0, campars(3, 0), campars(4, 0), campars(5, 0), campars(6, 0), campars(7, 0), campars(8, 0), campars(9, 0), kcoefs);
		}
		if (campars.rows() == 11) { //K1
			double kcoefs = 1;
			this_camera.ResetIOPs0(i, calibrated, campars(0, 0), campars(1, 0), campars(2, 0), campars(3, 0), 0.0, 0.0, 0.0, 0.0, campars(4, 0), campars(5, 0), campars(6, 0), campars(7, 0), campars(8, 0), campars(9, 0), campars(10, 0), kcoefs);
		}
		if (campars.rows() == 12) { //K1,K2
			double kcoefs = 2;
			this_camera.ResetIOPs0(i, calibrated, campars(0, 0), campars(1, 0), campars(2, 0), campars(3, 0), campars(4, 0), 0.0, 0.0, 0.0, campars(5, 0), campars(6, 0), campars(7, 0), campars(8, 0), campars(9, 0), campars(10, 0), campars(11, 0), kcoefs);
		}
		if (campars.rows() == 13) { //K1,K2,K3
			double kcoefs = 3;
			this_camera.ResetIOPs0(i, calibrated, campars(0, 0), campars(1, 0), campars(2, 0), campars(3, 0), campars(4, 0), campars(5, 0), 0.0, 0.0, campars(6, 0), campars(7, 0), campars(8, 0), campars(9, 0), campars(10, 0), campars(11, 0), campars(12, 0), kcoefs);
		}
		if (campars.rows() == 14) { //K1,K2,K3,K4
			double kcoefs = 4;
			this_camera.ResetIOPs0(i, calibrated, campars(0, 0), campars(1, 0), campars(2, 0), campars(3, 0), campars(4, 0), campars(5, 0), campars(6, 0), 0.0, campars(7, 0), campars(8, 0), campars(9, 0), campars(10, 0), campars(11, 0), campars(12, 0), campars(13, 0), kcoefs);
		}
		if (campars.rows() == 15) { //K1,K2,K3,K4,K5
			double kcoefs = 5;
			this_camera.ResetIOPs0(i, calibrated, campars(0, 0), campars(1, 0), campars(2, 0), campars(3, 0), campars(4, 0), campars(5, 0), campars(6, 0), campars(7, 0), campars(8, 0), campars(9, 0), campars(10, 0), campars(11, 0), campars(12, 0), campars(13, 0), campars(14, 0), kcoefs);
		}
	}


	for (int i = 0; i <= sens_model1 - 1; i++) {

		/*
		Cn
		Rn
		K1
		K2
		K3
		P1
		P2
		Cx (Width/2)
		Cy (Height/2)
		fx (f/pixelpitch)
		fy (f/pixelpitch)
		*/


		string camfile1 = "Cam_M1_" + to_string(i) + ".cam";
		char* camerapamfile = new char[camfile1.length() + 1];
		strcpy(camerapamfile, camfile1.c_str());

		MatrixXd campars;
		char currentfile[512] = "";
		currentfile[0] = '\0';
		strcat(currentfile, ImageFolder);
		strcat(currentfile, camerapamfile);
		Read_Mat(currentfile, campars);

		calibrated = true;
		this_camera.ResetIOPs1(i, calibrated, campars(0, 0), campars(1, 0), campars(2, 0), campars(3, 0), campars(4, 0), campars(5, 0), campars(6, 0), campars(7, 0), campars(8, 0), campars(9, 0), campars(10, 0));
	}


	for (int i = 0; i <= sens_model2 - 1; i++) {

		/*
		Cn
		Rn
		rho2
		rho3
		rho4
		Cx (width/2)
		Cy (height/2)
		C (f/pixelpicth*pi/2)
		D (usually 0)
		F (f/pixelpicth*pi/2)
		*/

		string camfile1 = "Cam_M2_" + to_string(i) + ".cam";
		char* camerapamfile = new char[camfile1.length() + 1];
		strcpy(camerapamfile, camfile1.c_str());

		MatrixXd campars;
		char currentfile[512] = "";
		currentfile[0] = '\0';
		strcat(currentfile, ImageFolder);
		strcat(currentfile, camerapamfile);
		Read_Mat(currentfile, campars);

		
		calibrated = true;
		//E is zero always
		this_camera.ResetIOPs2(i, calibrated, campars(0, 0), campars(1, 0), campars(2, 0), campars(3, 0), campars(4, 0), campars(5, 0), campars(6, 0), campars(7, 0), campars(8, 0), campars(9, 0));
	}


	for (int i = 0; i <= sens_model3 - 1; i++) {
		
		/*pixel pitch (mm)
		width (pix)
		height (pix)
		K1
		K2
		K3
		P1
		P2
		S1
		S2
		xpp
		ypp
		f (mm)*/

		string camfile1 = "Cam_M3_" + to_string(i) + ".cam";
		char* camerapamfile = new char[camfile1.length() + 1];
		strcpy(camerapamfile, camfile1.c_str());

		MatrixXd campars;
		char currentfile[512] = "";
		currentfile[0] = '\0';
		strcat(currentfile, ImageFolder);
		strcat(currentfile, camerapamfile);
		Read_Mat(currentfile, campars);
		

		calibrated = true;
		if (campars.rows() == 10) { //No K coefficient
			double kcoefs = 0;
			this_camera.ResetIOPs3(i, calibrated, campars(0, 0), campars(1, 0), campars(2, 0), 0.0, 0.0, 0.0, 0.0, 0.0, campars(3, 0), campars(4, 0), campars(5, 0), campars(6, 0), campars(7, 0), campars(8, 0), campars(9, 0), kcoefs);
		}
		if (campars.rows() == 11) { //K1
			double kcoefs = 1;
			this_camera.ResetIOPs3(i, calibrated, campars(0, 0), campars(1, 0), campars(2, 0), campars(3, 0), 0.0, 0.0, 0.0, 0.0, campars(4, 0), campars(5, 0), campars(6, 0), campars(7, 0), campars(8, 0), campars(9, 0), campars(10, 0), kcoefs);
		}
		if (campars.rows() == 12) { //K1,K2
			double kcoefs = 2;
			this_camera.ResetIOPs3(i, calibrated, campars(0, 0), campars(1, 0), campars(2, 0), campars(3, 0), campars(4, 0), 0.0, 0.0, 0.0, campars(5, 0), campars(6, 0), campars(7, 0), campars(8, 0), campars(9, 0), campars(10, 0), campars(11, 0), kcoefs);
		}
		if (campars.rows() == 13) { //K1,K2,K3
			double kcoefs = 3;
			this_camera.ResetIOPs3(i, calibrated, campars(0, 0), campars(1, 0), campars(2, 0), campars(3, 0), campars(4, 0), campars(5, 0), 0.0, 0.0, campars(6, 0), campars(7, 0), campars(8, 0), campars(9, 0), campars(10, 0), campars(11, 0), campars(12, 0), kcoefs);
		}
		if (campars.rows() == 14) { //K1,K2,K3,K4
			double kcoefs = 4;
			this_camera.ResetIOPs3(i, calibrated, campars(0, 0), campars(1, 0), campars(2, 0), campars(3, 0), campars(4, 0), campars(5, 0), campars(6, 0), 0.0, campars(7, 0), campars(8, 0), campars(9, 0), campars(10, 0), campars(11, 0), campars(12, 0), campars(13, 0), kcoefs);
		}
		if (campars.rows() == 15) { //K1,K2,K3,K4,K5
			double kcoefs = 5;
			this_camera.ResetIOPs3(i, calibrated, campars(0, 0), campars(1, 0), campars(2, 0), campars(3, 0), campars(4, 0), campars(5, 0), campars(6, 0), campars(7, 0), campars(8, 0), campars(9, 0), campars(10, 0), campars(11, 0), campars(12, 0), campars(13, 0), campars(14, 0), kcoefs);
		}
	}


	for (int i = 0; i <= sens_model4 - 1; i++) {

		/*
		Pixelpitch (mm)
		width (pix)
		height (pix)
		K1
		K2
		K3
		T1
		T2
		Px (mm)
		Py (mm)
		f (mm)
		*/

		string camfile1 = "Cam_M4_" + to_string(i) + ".cam";
		char* camerapamfile1 = new char[camfile1.length() + 1];
		strcpy(camerapamfile1, camfile1.c_str());

		MatrixXd campars;
		char currentfile[512] = "";
		currentfile[0] = '\0';
		strcat(currentfile, ImageFolder);
		strcat(currentfile, camerapamfile1);
		
		
		Read_Mat(currentfile, campars);

		calibrated = true;
		this_camera.ResetIOPs4(i, calibrated, campars(0, 0), campars(1, 0), campars(2, 0), campars(3, 0), campars(4, 0), campars(5, 0), campars(6, 0), campars(7, 0), campars(8, 0), campars(9, 0), campars(10, 0));
	}


	for (int i = 0; i <= sens_model6 - 1; i++) {

		//Format of this file:
		/*pixel pitch (mm)
		width (pix)
		height (pix)
		K1
		K2
		K3
		K4
		K5
		P1
		P2
		SS1
		SS2
		SS3
		SS4
		SS5
		xpp
		ypp
		f (mm)  */

		string camfile1 = "Cam_M6_" + to_string(i) + ".cam";
		char* camerapamfile = new char[camfile1.length() + 1];
		strcpy(camerapamfile, camfile1.c_str());

		MatrixXd campars;
		char currentfile[512] = "";
		currentfile[0] = '\0';
		strcat(currentfile, ImageFolder);
		strcat(currentfile, camerapamfile);
		Read_Mat(currentfile, campars);

		calibrated = true;

		if (campars.rows() == 13) { //No K coefficient
			double kcoefs = 0;
			this_camera.ResetIOPs6(i, calibrated, campars(0, 0), campars(1, 0), campars(2, 0),
				0.0, 0.0, 0.0, 0.0, 0.0,
				campars(3, 0), campars(4, 0),
				campars(5, 0), campars(6, 0), campars(7, 0),
				campars(8, 0), campars(9, 0),
				campars(10, 0), campars(11, 0), campars(12, 0), kcoefs);

		}
		if (campars.rows() == 14) { //K1
			double kcoefs = 1;
			this_camera.ResetIOPs6(i, calibrated, campars(0, 0), campars(1, 0), campars(2, 0),
				campars(3, 0), 0.0, 0.0, 0.0, 0.0,
				campars(4, 0), campars(5, 0),
				campars(6, 0), campars(7, 0), campars(8, 0),
				campars(9, 0), campars(10, 0),
				campars(11, 0), campars(12, 0), campars(13, 0), kcoefs);
		}
		if (campars.rows() == 15) { //K1,K2
			double kcoefs = 2;
			this_camera.ResetIOPs6(i, calibrated, campars(0, 0), campars(1, 0), campars(2, 0),
				campars(3, 0), campars(4, 0), 0.0, 0.0, 0.0,
				campars(5, 0), campars(6, 0),
				campars(7, 0), campars(8, 0), campars(9, 0),
				campars(10, 0), campars(11, 0),
				campars(12, 0), campars(13, 0), campars(14, 0), kcoefs);
		}
		if (campars.rows() == 16) { //K1,K2,K3
			double kcoefs = 3;
			this_camera.ResetIOPs6(i, calibrated, campars(0, 0), campars(1, 0), campars(2, 0),
				campars(3, 0), campars(4, 0), campars(5, 0), 0.0, 0.0,
				campars(6, 0), campars(7, 0),
				campars(8, 0), campars(9, 0), campars(10, 0),
				campars(11, 0), campars(12, 0),
				campars(13, 0), campars(14, 0), campars(15, 0), kcoefs);
		}
		if (campars.rows() == 17) { //K1,K2,K3,K4
			double kcoefs = 4;
			this_camera.ResetIOPs6(i, calibrated, campars(0, 0), campars(1, 0), campars(2, 0),
				campars(3, 0), campars(4, 0), campars(5, 0), campars(6, 0), 0.0,
				campars(7, 0), campars(8, 0),
				campars(9, 0), campars(10, 0), campars(11, 0),
				campars(12, 0), campars(13, 0),
				campars(14, 0), campars(15, 0), campars(16, 0), kcoefs);

		}
		if (campars.rows() == 18) { //K1,K2,K3,K4,K5
			double kcoefs = 5;
			this_camera.ResetIOPs6(i, calibrated, campars(0, 0), campars(1, 0), campars(2, 0),
				campars(3, 0), campars(4, 0), campars(5, 0), campars(6, 0), campars(7, 0),
				campars(8, 0), campars(9, 0),
				campars(10, 0), campars(11, 0), campars(12, 0),
				campars(13, 0), campars(14, 0),
				campars(15, 0), campars(16, 0), campars(17, 0), kcoefs);
		}
	}

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	int n_cams;
	vector<string> ImageList;
	MatrixXd MatchList;
	Data_Preparation(ImageFolder, DirectEO,
		this_camera, n_cams, imformat,
		ImageList, MatchFile, MatchList);
	
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
	for (int i = 0; i < ImageList.size(); i++) {
		cout << ImageList[i] << ", " << this_camera.model_ID(i, 0) << ", " << this_camera.sensor_ID(i, 0) << endl;
	}
	double general_Cn, general_Rn;
	vector<double> scale_factors(ImageList.size());
	for (int sensorid = 0; sensorid <= sens_model0 - 1; sensorid++) {
		general_Cn = this_camera.Cn_0(sensorid, 0);
		general_Rn = this_camera.Rn_0(sensorid, 0);
		ImageGeometry IMTYPE(general_Cn, general_Rn);
		for (int i = 0; i < ImageList.size(); i++) {
			
			if (this_camera.model_ID(i, 0) == 0 && this_camera.sensor_ID(i, 0) == sensorid) {
				IMTYPE.Undistort(ImageFolder, ImageList.at(i), this_camera, i, imformat);
				if (general_Cn > size_threshold && is_coregistration==false) {
					do_resize = true;
					scale_factor = size_threshold / general_Cn;
					scale_factors.at(i) = scale_factor;
				}
			}
		}
	}

	for (int sensorid = 0; sensorid <= sens_model1 - 1; sensorid++) {
		general_Cn = this_camera.Cn_1(sensorid, 0);
		general_Rn = this_camera.Rn_1(sensorid, 0);
		ImageGeometry IMTYPE(general_Cn, general_Rn);
		for (int i = 0; i < ImageList.size(); i++) {
			if (this_camera.model_ID(i, 0) == 1 && this_camera.sensor_ID(i, 0) == sensorid) {
				IMTYPE.Undistort(ImageFolder, ImageList.at(i), this_camera, i, imformat);
				if (general_Cn > size_threshold && is_coregistration == false) {
					do_resize = true;
					scale_factor = size_threshold / general_Cn;
					scale_factors.at(i) = scale_factor;
				}
			}
		}
	}

	for (int sensorid = 0; sensorid <= sens_model2 - 1; sensorid++) {
		general_Cn = this_camera.Cn_2(sensorid, 0);
		general_Rn = this_camera.Rn_2(sensorid, 0);
		ImageGeometry IMTYPE(general_Cn, general_Rn);
		for (int i = 0; i < ImageList.size(); i++) {
			if (this_camera.model_ID(i, 0) == 2 && this_camera.sensor_ID(i, 0) == sensorid) {
				IMTYPE.Undistort(ImageFolder, ImageList.at(i), this_camera, i, imformat);
				if (general_Cn > size_threshold && is_coregistration == false) {
					do_resize = true;
					scale_factor = size_threshold / general_Cn;
					scale_factors.at(i) = scale_factor;
				}
			}
		}
	}

	for (int sensorid = 0; sensorid <= sens_model3 - 1; sensorid++) {
		general_Cn = this_camera.Cn_3(sensorid, 0);
		general_Rn = this_camera.Rn_3(sensorid, 0);
		ImageGeometry IMTYPE(general_Cn, general_Rn);
		for (int i = 0; i < ImageList.size(); i++) {
			if (this_camera.model_ID(i, 0) == 3 && this_camera.sensor_ID(i, 0) == sensorid) {
				IMTYPE.Undistort(ImageFolder, ImageList.at(i), this_camera, i, imformat);
				if (general_Cn > size_threshold && is_coregistration == false) {
					do_resize = true;
					scale_factor = size_threshold / general_Cn;
					scale_factors.at(i) = scale_factor;
				}
			}
		}
	}

	for (int sensorid = 0; sensorid <= sens_model4 - 1; sensorid++) {
		general_Cn = this_camera.Cn_4(sensorid, 0);
		general_Rn = this_camera.Rn_4(sensorid, 0);
		ImageGeometry IMTYPE(general_Cn, general_Rn);
		for (int i = 0; i < ImageList.size(); i++) {
			if (this_camera.model_ID(i, 0) == 4 && this_camera.sensor_ID(i, 0) == sensorid) {
				IMTYPE.Undistort(ImageFolder, ImageList.at(i), this_camera, i, imformat);
				if (general_Cn > size_threshold && is_coregistration == false) {
					do_resize = true;
					scale_factor = size_threshold / general_Cn;
					scale_factors.at(i) = scale_factor;
				}
			}
		}
	}

	for (int sensorid = 0; sensorid <= sens_model6 - 1; sensorid++) {
		general_Cn = this_camera.Cn_00(sensorid, 0);
		general_Rn = this_camera.Rn_00(sensorid, 0);
		ImageGeometry IMTYPE(general_Cn, general_Rn);
		for (int i = 0; i < ImageList.size(); i++) {
			if (this_camera.model_ID(i, 0) == 6 && this_camera.sensor_ID(i, 0) == sensorid) {
				IMTYPE.Undistort(ImageFolder, ImageList.at(i), this_camera, i, imformat);
				if (general_Cn > size_threshold && is_coregistration == false) {
					do_resize = true;
					scale_factor = size_threshold / general_Cn;
					scale_factors.at(i) = scale_factor;
				}
			}
		}
	}


	Mat _IML, _IMR;
	Matrix3b3 R_common;
	MatrixXd K_common1, K_common2;
	Matrix<int, Dynamic, Dynamic> Mask;
	MatrixXd DLeft, H1i_transform, H2i_transform;
	Matrix3b4 Pmat1, Pmat2;

	int pc_counter = 0;
	vector<PCpoint> PCVECT;

	for (int ii = 0; ii < MatchList.rows(); ii++) {
		for (int jj = 0; jj < MatchList.cols(); jj++) {
			if (MatchList(ii, jj) == 1) {

				string imright = ImageList.at(jj); string imleft = ImageList.at(ii);
				int iright = jj; int ileft = ii;

				bool res=Rectifier(ImageFolder, imleft, imright, ileft, iright, this_camera,
					process_sizex, process_sizey, R_common, K_common1, K_common2, scale_factors, do_resize, process_in_batch,
					Mask, imformat,
					H1i_transform, H2i_transform,
					Pmat1, Pmat2,
					is_coregistration);

				if (res == true) {
					string ImageNameL = imleft;
					string ImageNameR = imright;
					string foldername = string(ImageFolder);
					int found = ImageNameL.find('.');
					string leftname = (ImageNameL.substr(0, found));
					found = ImageNameR.find('.');
					string rightname = (ImageNameR.substr(0, found));

					string lname = foldername + leftname + "_" + rightname + "_first.pgm";
					string rname = foldername + leftname + "_" + rightname + "_second.pgm";
					char * lnamecstr = new char[lname.length() + 1];
					strcpy(lnamecstr, lname.c_str());
					char * rnamecstr = new char[rname.length() + 1];
					strcpy(rnamecstr, rname.c_str());

					Elas_process(DLeft, lnamecstr, rnamecstr, Mask); //for Mozhdeh's method
					//PatchMatch_process(DLeft, lname, rname); //for PatchMatch method (this one needs to be changed as it is super slow now even for low-res)

					if (is_coregistration == false)
					{
						PointCloud(ImageFolder, imleft, ileft, iright, this_camera,
							DLeft, SaveName,
							delta_trans, R_common, K_common1, K_common2, Mask, maxdist_im2cloud,
							PCVECT, pc_counter, imformat);
					}

					if (is_coregistration == true) {

						MatrixXd K3, K2, K1, R3, R2, R1, T3, T2, T1, A3, A2, A1, B3, B2, B1, C3, C2, C1;
						Matrix3b4 P3, P2, P1;
						double Rn;

						get_K_R_T_P(this_camera, ileft, false, 1, Rn, K1, R1, T1, P1);
						A1 = (P1.row(0).transpose()); B1 = (P1.row(1).transpose()); C1 = (P1.row(2).transpose());
						get_K_R_T_P(this_camera, iright, false, 1, Rn, K2, R2, T2, P2);
						A2 = (P2.row(0).transpose()); B2 = (P2.row(1).transpose()); C2 = (P2.row(2).transpose());

						MatrixXf Overlay_1x(DLeft.rows(), DLeft.cols());
						MatrixXf Overlay_1y(DLeft.rows(), DLeft.cols());
						MatrixXf Overlay_2x(DLeft.rows(), DLeft.cols());
						MatrixXf Overlay_2y(DLeft.rows(), DLeft.cols());

						for (int i = 0; i < DLeft.rows(); i++) {
							for (int j = 0; j < DLeft.cols(); j++) {
								double x1 = j;
								double y1 = i;
								double disp = DLeft(i, j);
								if (disp != -10 && Mask(i, j) == 255) {
									double x2 = x1 - disp;
									double y2 = y1;

									x1 += H1i_transform(0, 3);
									y1 += H1i_transform(1, 3);
									MatrixXd tmp(3, 1); tmp << x1, y1, 1;
									tmp = H1i_transform.block(0, 0, 3, 3) * tmp;
									x1 = (tmp(0, 0) / tmp(2, 0));
									y1 = (tmp(1, 0) / tmp(2, 0));

									x2 += H2i_transform(0, 3);
									y2 += H2i_transform(1, 3);
									tmp << x2, y2, 1;
									tmp = H2i_transform.block(0, 0, 3, 3) * tmp;
									x2 = (tmp(0, 0) / tmp(2, 0));
									y2 = (tmp(1, 0) / tmp(2, 0));

									Overlay_1x(i, j) = x1;
									Overlay_2x(i, j) = x2;
									Overlay_1y(i, j) = y1;
									Overlay_2y(i, j) = y2;
								}
							}
						}
						vector<int> imwrite_jpg_params;
						imwrite_jpg_params.push_back(IMWRITE_JPEG_QUALITY);
						imwrite_jpg_params.push_back(100);

						string imfullname = foldername + "Undist_" + leftname + ".tiff";
						string newname = foldername + "Lay_" + leftname + ".tiff";
						Mat Image1 = imread(imfullname, IMREAD_COLOR | IMREAD_ANYDEPTH | IMREAD_IGNORE_ORIENTATION);
						Mat Imcorrect1(Image1.size(), Image1.type());
						Mat X1(Image1.size(), CV_32FC1);
						Mat Y1(Image1.size(), CV_32FC1);
						eigen2cv(Overlay_1x, X1);
						eigen2cv(Overlay_1y, Y1);
						remap(Image1, Imcorrect1, X1, Y1, INTER_CUBIC, 0, Scalar(255, 255, 255));
						X1.release();
						Y1.release();
						Image1.release();
						imwrite(newname, Imcorrect1);//, imwrite_jpg_params
						Imcorrect1.release();

						imfullname = foldername + "Undist_" + rightname + ".tiff";
						newname = foldername + "Lay_" + rightname + ".tiff";
						Mat Image2 = imread(imfullname, IMREAD_COLOR | IMREAD_ANYDEPTH | IMREAD_IGNORE_ORIENTATION);
						Mat Imcorrect2(Image2.size(), Image2.type());
						Mat X2(Image2.size(), CV_32FC1);
						Mat Y2(Image2.size(), CV_32FC1);
						eigen2cv(Overlay_2x, X2);
						eigen2cv(Overlay_2y, Y2);
						remap(Image2, Imcorrect2, X2, Y2, INTER_CUBIC, 0, Scalar(255, 255, 255));
						X2.release();
						Y2.release();
						Image2.release();
						imwrite(newname, Imcorrect2);//, imwrite_jpg_params
						Imcorrect2.release();

						size_t pos = imleft.find("REG.TIF");
						string postfix, prefix;
						
						postfix= "GRE";
						prefix = imleft.substr(0, pos);
						TriFocal_Overlay(prefix, postfix, imformat, foldername, ImageList, this_camera,
							K1, K2, T1, T2, R1, R2,
							P1, P2,
							A1, A2, B1, B2, C1, C2,
							DLeft, Mask,
							Overlay_1x, Overlay_1y, Overlay_2x, Overlay_2y);

						postfix = "RED";
						prefix = imleft.substr(0, pos);
						TriFocal_Overlay(prefix, postfix, imformat, foldername, ImageList, this_camera,
							K1, K2, T1, T2, R1, R2,
							P1, P2,
							A1, A2, B1, B2, C1, C2,
							DLeft, Mask,
							Overlay_1x, Overlay_1y, Overlay_2x, Overlay_2y);

					}
				}
			}
		}
	}

	if (!PCVECT.empty()) {
		string foldername2 = string(ImageFolder);
		string fname = foldername2 + SaveName + to_string(pc_counter) + ".csv";
		WriteCSV(fname, PCVECT);
		pc_counter++;
		PCVECT.clear();
	}
	


	if (!Keep_undistored) {
		cout << "Removing the undistorted images ..." << endl;
		for (int i = 0; i < ImageList.size(); i++) {
			ImageGeometry_Remove(ImageFolder, ImageList.at(i), imformat);
		}
	}

	system("pause");
	return 1;
}
