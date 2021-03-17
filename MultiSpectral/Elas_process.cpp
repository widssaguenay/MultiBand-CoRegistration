
#include "Coregistration.h"

void Elas_process(MatrixXd& DLeft, const char* file_1, const char* file_2, Matrix<int, Dynamic, Dynamic>& Mask) {

	cout << "Processing: " << file_1 << ", " << file_2 << endl;


	image<uchar> *I1, *I2;
	I1 = loadPGM(file_1);
	I2 = loadPGM(file_2);

	int32_t width = I1->width();
	int32_t height = I1->height();

	const int32_t dims[3] = { width,height,width }; 
	float* D1_data = (float*)malloc(width*height * sizeof(float));
	float* D2_data = (float*)malloc(width*height * sizeof(float));

	Elas::parameters param(Elas::setting::ROBOTICS); //
	param.postprocess_only_left = false;
	param.disp_min = 0;
	param.disp_max = width;
	param.lr_threshold = 0.5;

	Elas elas(param);
	elas.process(I1->data, I2->data, D1_data, D2_data, dims);
	int b = remove(file_1); b = remove(file_2);


	DLeft.setZero((int)height, (int)width);
	int32_t counter = -1;
	for (int32_t i = 0; i<height; i++) {
		for (int32_t j = 0; j < width; j++) {
			counter++;
			DLeft(i, j)=(double)D1_data[counter];
			
		}
	}


	delete I1;
	delete I2;
	
	free(D1_data);
	free(D2_data);

	

	Mat X(height, width, CV_32FC1);
	MatrixXf MapXf = (DLeft.cast<float>());
	eigen2cv(MapXf, X);
	

	Mat Msk = Mat::zeros(height, width, CV_8UC1);
	for (int y1 = 0; y1 < DLeft.rows(); y1++) {
		for (int x1 = 0; x1 < DLeft.cols(); x1++) {
			double disp = DLeft(y1, x1);
			if (disp == -10) {
				Msk.at<uint8_t>(y1, x1) = (uint8_t)255;
			}
		}
	}


	/*vector<int> imwrite_jpg_params;
	imwrite_jpg_params.push_back(IMWRITE_JPEG_QUALITY);
	imwrite_jpg_params.push_back(100);
	imwrite("Mask.jpg", Msk, imwrite_jpg_params);*/


	Mat dst;
	inpaint(X, Msk, dst, 3, INPAINT_NS);

	cv2eigen(dst, MapXf);
	DLeft= (MapXf.cast<double>());
	MapXf.resize(1, 1);
	X.release();
	Msk.release();
	dst.release();

	//temporarily for test
	//char name[256] = ".\\Disparity.txt";
	//Write_Mat(name, DLeft, 1);
};
