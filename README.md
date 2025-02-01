# Learning Metric Fields for Fast Low-Distortion Mesh Parameterizations
--------------------------------------------------------------------------------------------------------------------------------------------------------
This code includes the implementation of the Eurographics 2025 paper "Learning Metric Fields for Fast Low-Distortion Mesh Parameterizations" authored by Guy Fargion and Ofir Weber. 

The use of this application is limited to academic use only!

The code is provided as-is and without any guarantees.

This GIF submodule is partially based on https://github.com/GuyFa/GIF

The CM submodule is partially based on https://github.com/Roipo/CompMajor

The LMF submodule is partially based on https://github.com/ThibaultGROUEIX/NeuralJacobianFields

For questions or comments about the code please contact:
Guy Fargion (guy.fargion@gmail.com)

----------------------------------------------------------------------------
The code should be platform independent though we never tested it on other than Windows OS.

## Installlation

### 1) GIF
A Visual Studio 2019 project is provided for easy compilation on Windows machines.

The following prerequisites are necessary for building and running the code:

1) Matlab R2022b

2) Boost 1.83.0 - We downloaded boost from here https://www.boost.org/. Our code only requires the headers of the Boost libraries. Hence, there is no need to build boost.

3) CGAL 5.6 - We installed CGAL 5.6 with the gmp and mpfr auxiliary libraries. There is no need to build CGAL as well.

4) GMM C++ template library version 4.2 (http://getfem.org/download.html).

5) PARDISO 8.2

6) Eigen 3.4.0


Other versions of the above listed tools might be compatible but weren't tested!

Now:

1) Install and build the above mentioned prerequisites.

2) Add the following environment variables to your system (see some possible paths):

GMM_INCLUDE_DIR		  (%your GMM folder path%)\gmm-4.2\include

MATLAB_64_DIR		    C:\Program Files\MATLAB\R2022b

CGAL_64_DIR		      C:\Program Files\CGAL-5.6

BOOST_64_DIR		    C:\Program Files\boost\boost_1_83_0

EIGEN_DIR  		      (%your Eigen folder path%)

PARDISO_BIN         (%path to the folder where PARDISO dll and lib are%)

PARDISO_LIC_PATH    (%path to the folder with PARDISO licence%)

OMP_NUM_THREADS  	  number of cores in your CPU (for PARDISO)

3) Add the the folder "MatlabScripts" to your Matlab path.

4) Make sure all the required dlls can be loacted by including the relevant paths into the system PATH variable.
For example:

PATH=
%MATLAB_64_DIR%\bin\win64;
%MATLAB_64_DIR%\extern\include;
%MATLAB_64_DIR%\extern\lib\win64\microsoft;
%BOOST_64_DIR%\libs;
%CGAL_64_DIR%\lib;
%CGAL_64_DIR%\auxiliary\gmp\lib;
%CGAL_64_DIR%\include\CGAL;
%PARDISO_BIN%;
%MATLAB_64_DIR%

### 2) CM

After installing GIF, CM should run as well with no further installation, as the environment variables needed for PARDISO already exist. Make sure you clone all the recursive submodules as well.

### 3) LMF
We used python 3.10.

All requirements appear in ```LMF/requirements.txt```.

The torch installation command:
```shell
pip3 install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu121
```

## Training
```LMF/script_train_default.py``` contains a script for both training and inference.

Each training/validation database folder is assumed to contain pairs of 3D disk-like meshes, and their GT result:

```
database
|model_0_source.obj
|model_0_target.obj
|- ...
|model_999_source.obj
|model_999_target.obj
```

and an additional file named ```data.file``` with the list of pairs:

```
{
"pairs": [
        [
            "model_0_source.obj",
            "model_0_target.obj"
        ],
        ...
        [
            "model_999_source.obj",
            "model_999_target.obj"
        ],
    ]
}
```

Now just set the appropriate paths in the arguments ```--root-dir-train``` and ```--root-dir-validation```, and run the script in ```train``` mode.

## Inference
Create the test database, in the same format as in the training/validation databases, set the argument ```--root-dir-test``` with the appropriate path, and run the script in ```test``` mode.


The full list of arguments is available in ```LMF/source_lmf/args_from_cli.py```.

The pre-trained network weights is available at https://livebiuac-my.sharepoint.com/:f:/g/personal/guy_fargion_live_biu_ac_il/Er6V-EH2Ni9Ksuh-wvDQduIBV8TFpDeCJyrzAlXi1IiRhA?e=98c3mo. Download the content of this shared folder into the folder ```LMF/network weights```.

## Object files
The training, validation and test sets along with their GTs are available at https://livebiuac-my.sharepoint.com/:f:/g/personal/guy_fargion_live_biu_ac_il/EvTjyBuS1S9Lk8BqmX1KSpYBYxU4MjVdga17W2H6wg8L9Q?e=MTq5cS

The models from the figures in the paper are available at https://livebiuac-my.sharepoint.com/:f:/g/personal/guy_fargion_live_biu_ac_il/EmQudZIFn6dMu-VPLGnwjfABY-A_5LX4OMvyaUWkfyyY5w?e=bAygeP
