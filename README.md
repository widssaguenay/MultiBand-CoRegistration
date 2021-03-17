# MultiSpectral-CoRegistration
Overlays images captured by cameras fixed in a rig, e.g. Parrot Sequoia and RedEdge MicaSense
This is an implementation of the following paper:

Shahbazi, M. and Cortes, C., 2019. Seamless Co-Registration of Images from Multi-Sensor Multispectral Cameras. The International Archives of Photogrammetry, Remote Sensing and Spatial Information Sciences, 42, pp.315-322.

## Requirements
The code is written in C++17 and tested only on Windows Visual Studio 2019
Dependencies include Eigen, LAPACK (prebuilt libraries are included in this repo), and OpenCV4.2

## Test Data
A folder including test data from a Sequoia multi-spectral camera is included.
The sample command line arguments for this dataset is:

MultiSpectral.exe --folder .\\Test\\ --model0 5 --model1 0 --model2 0 --model3 0 --model4 0 --model6 0 --imformat TIF --name PointCloud --size 2000 –keep 1 --maxdist 300 --coreg 1

