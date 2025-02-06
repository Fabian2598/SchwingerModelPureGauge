Simulation of the Schwinger Model without fermions. The model is simulated according to the action 

$$S_G\left(U\right) = \beta \sum_{\textbf{n}\in \Lambda} \Re \left(1 - U_{01}(\textbf{n})\right),$$
where $$U_{01} =U_{0}(\textbf{n}) \cdot U_{1}(\textbf{n}+\hat{0}) \cdot U_{0}^\dagger(\textbf{n}+\hat{1}) \cdot U_{1}^\dagger(\textbf{n}) $$ 

is the plaquette.

To compile create a new folder 

```
mkdir build
```

## Linux
If you are working on Linux you have to delete lines 4 and 5 in the CMakeLists.txt file. 
These set the compiler direction which is necessary for Windows. Then, in the `build` folder, run the following commands:

```
cmake ../
cmake --build .
```

This will create an executable for you to run. The lattice dimensions are fixed in the CMakeLists.txt, on lines 19-20. 
You can change the dimensions there as well as the executable name.

A running example with HMC is shown below

```
./SM_NsxNt.exe
-----------------------
|Pure Gauge U(1) theory|
-----------------------
Ns 16 Nt 16
Algorithm (0 Metropolis, 1 HMC): 1
----HMC----
Molecular dynamics steps: 10
Trajectory length: 1
beta min: 0.2
beta max: 2
Number of betas: 10
Thermalization: 1000
Measurements: 1000
Step (sweeps between measurements): 10
```

For `Algorithm: 0` the trajectory is ignored.

## Windows

The instructions are essentially the same. The CMakeLists.txt only needs the address of your C++ and C compiler on lines 4 and 5. 
Then, in the `build` folder, run the following commands:
```
cmake -G "MinGW Makefiles" -DCMAKE_CXX_COMPILER=C:\msys64\ucrt64\bin\g++ -DCMAKE_C_COMPILER=C:\msys64\ucrt64\bin\gcc ../
```

 This command depends on the compiler you are using. In this case, we are using MinGW. If you are using another compiler 
 you have to change the `-G` flag. The `-DCMAKE_CXX_COMPILER` and `-DCMAKE_C_COMPILER` flags are the address of the compiler.

Then you can run the executable

```
SM_NsxNt.exe
```

**Only the average plaquette value is measured, one can implement other observables.**