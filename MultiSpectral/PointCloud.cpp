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

////////////////////////////////////
size_t WriteByBuffer(vector<char>& buffer, size_t offset, double number, char separator)
{
	auto end_pt = to_chars(buffer.data() + offset, buffer.data() + buffer.size(), number);

	if (end_pt.ec == errc::value_too_large)
	{
		throw runtime_error("failed writing");
	}

	if (end_pt.ptr >= buffer.data() + buffer.size())
	{
		throw runtime_error("failed writing");
	}

	*end_pt.ptr = separator;
	return static_cast<size_t>((end_pt.ptr + 1) - buffer.data());
};

void WriteCSV(const string& csv_filepath, const vector<PCpoint>& points) {
	ofstream ofs;
	ofs.exceptions(ifstream::badbit);
	ofs.open(csv_filepath, std::ios::out | std::ios::trunc); //

	constexpr char separ = ',';
	constexpr char separ_end = '\n';
	ofs << "X" << separ;
	ofs << "Y" << separ;
	ofs << "Z" << separ;
	ofs << "RED" << separ;
	ofs << "GREEN" << separ;
	ofs << "BLUE" << separ;
	ofs << "IM1" << separ;
	ofs << "IM2" << separ_end;


	constexpr size_t BUFFER_SIZE = 1 << 11; // 1.1kB
	vector<char> buffer(BUFFER_SIZE, '\0');

	for (const PCpoint& point : points)
	{
		
			size_t offset = 0;
			offset = WriteByBuffer(buffer, offset, point.x, separ);
			offset = WriteByBuffer(buffer, offset, point.y, separ);
			offset = WriteByBuffer(buffer, offset, point.z, separ);
			offset = WriteByBuffer(buffer, offset, point.r, separ);
			offset = WriteByBuffer(buffer, offset, point.g, separ);
			offset = WriteByBuffer(buffer, offset, point.b, separ);
			offset = WriteByBuffer(buffer, offset, point.im1, separ);
			offset = WriteByBuffer(buffer, offset, point.im2, separ_end);


			ofs << string_view(buffer.data(), offset);
		
	}

	ofs.close();
	return;
};


void PointCloud(char* FeatSavePP, string ImageNameL, int im1, int im2, CameraParam& im_cam,
	MatrixXd& DLeft, string SaveName,
	MatrixXd& delta_trans, Matrix3b3& R_common, MatrixXd& K_common1, MatrixXd& K_common2,
	Matrix<int, Dynamic, Dynamic>& Mask, double maxdist_im2cloud,
	vector<PCpoint>& pcvec,int& pc_counter,
	string imformat) {

	Matrix3b3 r;
	Matrix3b4 RT;
	Matrix3b4 Pmat1;
	MatrixXd Te2i_1;
	MatrixXd Te2i_2;
	Matrix3b3	Rxe2i_1;

	int im = im1;
	int sensor_id = (int)im_cam.sensor_ID(im, 0);
	int model_id = (int)im_cam.model_ID(im, 0);
	double PS, Cx, Cy, f, PSx, PSy, Cn1, Rn1;
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
	MatrixXd A(3, 3);
	A << PSx, 0.0, Cx,
		0.0, -PSy, Cy,
		0.0, 0.0, -f;
	MatrixXd Kmat1 = (A.inverse());
	Kmat1 = Kmat1 / abs(Kmat1(2, 2));
	r = im_cam.Bundestim_ori.block(im * 5 + 0, 0, 3, 3);
	Rxe2i_1 = r;
	Te2i_1 = ((im_cam.Bundestim_ori.row(im * 5 + 3)-delta_trans).transpose());
	RT.block(0, 0, 3, 3) = r;
	RT.block(0, 3, 3, 1) = -r * Te2i_1;
	Pmat1 = Kmat1 * RT;
	Te2i_2 = ((im_cam.Bundestim_ori.row(im2 * 5 + 3) - delta_trans).transpose());

	string foldername = string(FeatSavePP);
	int found = ImageNameL.find(imformat);
	string imfullname = foldername + "Undist_" + (ImageNameL.substr(0, found)); 
	imfullname += ".tiff";
	Mat Image = imread(imfullname, IMREAD_COLOR | IMREAD_ANYDEPTH);
	if (Image.empty())
	{
		cout << "Cannot read image: " << imfullname << std::endl;
		return;
	}
	Mat bgr[3]; 
	split(Image, bgr);
	Matrix<int,Dynamic,Dynamic> Blm, Grm, Rem;
	cv2eigen(bgr[0], Blm);
	cv2eigen(bgr[1], Grm);
	cv2eigen(bgr[2], Rem);
	

	RT.block(0, 0, 3, 3) = R_common;
	RT.block(0, 3, 3, 1) = -R_common * ( Te2i_1);
	MatrixXd Pmat1_common= K_common1* RT;
	RT.block(0, 0, 3, 3) = R_common;
	RT.block(0, 3, 3, 1) = -R_common * ( Te2i_2);
	MatrixXd Pmat2_common = K_common2 * RT;

	MatrixXd A_mat, L_mat;
	MatrixXd tmp_xy;
	MatrixXd tmp_X(4, 1);
	MatrixXd TMP31;
	MatrixXd some;
	MatrixXd disttemp;

	for (double y1 = 0; y1 < DLeft.rows(); y1++) {
		for (double x1 = 0; x1 < DLeft.cols(); x1++) {
			double disp = DLeft((int)y1, (int)x1);
			if (disp != -10 && Mask((int)y1, (int)x1)==255) {
				double x2 = x1 - disp;
				double y2 = y1;

				A_mat.setZero(4, 3);
				L_mat.setZero(4, 1);

				A_mat.row(0) << Pmat1_common(0, 0) - x1 * Pmat1_common(2, 0), Pmat1_common(0, 1) - x1 * Pmat1_common(2, 1), Pmat1_common(0, 2) - x1 * Pmat1_common(2, 2);
				A_mat.row(1) << Pmat1_common(1, 0) - y1 * Pmat1_common(2, 0), Pmat1_common(1, 1) - y1 * Pmat1_common(2, 1), Pmat1_common(1, 2) - y1 * Pmat1_common(2, 2);

				L_mat.row(0) << -Pmat1_common(0, 3) + Pmat1_common(2, 3)*x1;
				L_mat.row(1) << -Pmat1_common(1, 3) + Pmat1_common(2, 3)*y1;

				A_mat.row(2) << Pmat2_common(0, 0) - x2 * Pmat2_common(2, 0), Pmat2_common(0, 1) - x2 * Pmat2_common(2, 1), Pmat2_common(0, 2) - x2 * Pmat2_common(2, 2);
				A_mat.row(3) << Pmat2_common(1, 0) - y2 * Pmat2_common(2, 0), Pmat2_common(1, 1) - y2 * Pmat2_common(2, 1), Pmat2_common(1, 2) - y2 * Pmat2_common(2, 2);

				L_mat.row(2) << -Pmat2_common(0, 3) + Pmat2_common(2, 3)*x2;
				L_mat.row(3) << -Pmat2_common(1, 3) + Pmat2_common(2, 3)*y2;

				some = (A_mat.transpose());
				TMP31 = (some*A_mat).lu().solve(some*L_mat);
				disttemp = TMP31 - Te2i_1;
				if ( (disttemp.norm()) < maxdist_im2cloud) {
					tmp_X << TMP31(0, 0), TMP31(1, 0), TMP31(2, 0), 1.0;
					tmp_xy = Pmat1 * tmp_X;
					tmp_xy /= tmp_xy(2, 0);
					int x = round(tmp_xy(0, 0));
					int y = round(tmp_xy(1, 0));

					if (x > 0 && y > 0 && x < Image.cols && y < Image.rows) {

						PCpoint pcpt;
						pcpt.x = TMP31(0, 0) + delta_trans(0, 0);
						pcpt.y = TMP31(1, 0) + delta_trans(0, 1);
						pcpt.z = TMP31(2, 0) + delta_trans(0, 2);
						pcpt.r = Rem(y, x);
						pcpt.g = Grm(y, x);
						pcpt.b = Blm(y, x);
						pcpt.im1 = im1;
						pcpt.im2 = im2;
						pcvec.push_back(pcpt);
						
					}
				}
			}
		}
	}
	
	if (pcvec.size() > 1000000) {
		string foldername2 = string(FeatSavePP);
		string fname = foldername2 + SaveName + to_string(pc_counter)+ ".csv";
		WriteCSV(fname, pcvec);
		pc_counter++;
		pcvec.clear();
	}
};