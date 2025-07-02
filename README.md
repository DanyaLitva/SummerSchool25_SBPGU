## How to build this project:

**On windows, with Visual Studio 17 2022 and CMake**
1. Install and setup [MS-MPI](https://docs.microsoft.com/en-us/message-passing-interface/microsoft-mpi-release-notes) for Visual Studio
2. Clone this repository with `git clone https://github.com/DanyaLitva/SummerSchool25_SBPGU`
3. Enter the directory with this project and run commands below:
    ```
    mkdir build
    cd build
    cmake ../
    cmake --build . --config Release
    ```
To run tests, type `mpiexec -n 4 .\bin\test_project.exe`  
To run sample, type `mpiexec -n 4 .\bin\main.exe`


**On linux, with CMake**
1. Install mpi with `sudo apt install mpich`
2. Clone this repository with `git clone https://github.com/DanyaLitva/SummerSchool25_SBPGU`
3. Enter the directory with this project and run commands below:
    ```
    mkdir build
    cd build
    cmake ../
    cmake --build . --config Release
    ```
To run tests, type `mpirun -n 4 ./bin/test_project`  
To run sample, type `mpirun -n 4 ./bin/main`
