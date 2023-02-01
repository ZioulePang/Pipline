To compile on University Linux:

module load legacy-eng
module add gcc
module add qt/5.13.0
ln -s /usr/lib64/libGL.so.1 lib/libGL.so
qmake
make
--------------------------------------

To compile on OSX:
Use Homebrew to install qt

qmake -project QT+=opengl
qmake
make
--------------------------------------

To compile on Windows:
Install QT 5.13.0 
It can be found here: 
	https://download.qt.io/archive/qt/5.13/5.13.0/qt-opensource-windows-x86-5.13.0.exe.mirrorlist

Update graphics card driver
Double click LeedsGLRenderWindow.pro to open in QTCreator

Select the platform to compile to (32 or 64 bits)
Click details to select the build folder
Click Configure Project
--------------------------------------

To run on feng-linux / feng-gps:
./LeedsGLRenderWindow ../path_to/model.obj ../path_to/texture.ppm

---------------------------------------------

To run on OSX:
./LeedsGLRenderWindow.app/Contents/MacOS/RaytraceRenderWindowRelease  ../path_to/model.obj ../path_to/texture.ppm

---------------------------------------------
To run on Windows
./LeedsGLRenderWindow.exe ../path_to/model.obj ../path_to/texture.ppm
or
Click projects
Select "Run" on the left side menu under the active configuration
Add "../path_to/model.obj ../path_to/texture.ppm" to command line arguments
Click Run/Debug on the bottom left side.