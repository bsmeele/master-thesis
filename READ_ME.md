Simulator for a memristor crossbar.

This simulator expects the Eigen library (https://eigen.tuxfamily.org/index.php?title=Main_Page) to be cloned in the `crossbar_model` subfolder

The memristor model is developed by https://ieeexplore-ieee-org.tudelft.idm.oclc.org/document/9181475
Large parts are a recreation of the Verilog-A code: https://www.emrl.de/JART.html#Artikel_3

The linear crossbar solver is developed by: https://ieeexplore-ieee-org.tudelft.idm.oclc.org/document/6473873

The intended interface for this simulator is the `CrossbarSimulator` class in `crossbar_simulator.h`
The `CrossbarSimulator` class has the following fields:
- `M`: the height of the crossbar
- `N`: the width of the crossbar
- `RRAM`: a M by N matrix containing instances of memristor classes. Initialized with default parameters (as defined in `JART_VCM_v1b_var.h`)
- `access_transistors`: a M by N matrix containing access transistors for each memristor. When set to false, a memristor is disconnected and acts as a resistor with infinite resistance for the purposes of simulation. Initialized as all true
- `Rswl1`: the access resistance for the left side of the crossbar. Initialized as 1 Ohm
- `Rswl2`: the access resistance for the right side of the crossbar. Initialized as 1 Ohm
- `Rsbl1`: the access resistance for the top side of the crossbar. Initialized as 1 Ohm
- `Rsbl2`: the access resistance for the bottom side of the crossbar. Initialized as 1 Ohm
- `Rwl`: the wordline resistance. Initialized as 1 Ohm
- `Rbl`: the bitline resistance. . Initialized as 1 Ohm
- `partial_G_ABCD`: used for crossbar calculations
- `linear_solver`: the solver used in crossbar calculations
The `CrossbarSimulator` class has the following methods:
- `SetRRAM()`: sets the memristors to their maximum or minimum internal state (thus setting to HRS or LRS) based on the provided weights
- `NonlinearSolve()`: calculates the nodal voltages for the provided input voltages
- `ApplyVoltage()`: calculates the output current for each individual memristor and updates the memristors for the provided input voltages and time step
- `CalculateIout()`: calculates the output current of the crossbar for the provided nodal voltages

The call order of `ApplyVoltage()` is roughly as follows:
- `NonlinearSolve()`
    - `FixedpointSolve()` (or other based on the method parameter, it is highly recommended to use the fixed point solver)
        - `memristor.GetResistance()`
        - `SolveCAM()`
- `memristor.ApplyVoltage()` (this function returns already returns the current for the individual memristor)
- Return the output current

The simulator can roughly be devided into three subcomponents:
- The memristor model in `memristor_model/JART_VCM_v1b_var.cpp`
- The linear crossbar solver in `crossbar_model/linear_crossbar_solver.cpp`
- The non-linear crossbar solver in `nonlinear_crossbar_solver.cpp`

For the `CrossbarSimulator` class as well as each subcomponent, an example use case in provided:
- `RRAM_validation.cpp` simulates applying voltage pulses to a crossbar and calculates the averaged individual memristor current as well as the averaged output current
- `main.cpp` simulates multiple non-linear crossbars with randomized weights and averages their execution time
- `memristor_model/main.cpp` simulates a single memristor to which a triangle wave is applied. `memristor_model/plot.py` can be used to plot the simulation results
- `crossbar_model/main.cpp` simulates multiple linear crossbars with randomized weights and averages their execution time

Some additional notes on `RRAN_validation.cpp`:
- This function expects a directory as a command line argument. Example: ./RRAM_validation.exe row_448-479_col_448-479_pos/
- This folder is expected to include:
    - input.bin, which includes a 10x32 matrix of input voltages
    - weight.bin, which includes a 32x32 matrix of weights
- After execution the following files will be written to the same folder:
    - outupt.bin, a 10x32x32 matrix containing the individual memristive currents for each input voltage
    - MAC.bin, a 10x32 matrix containing the crossbar output currents for each input voltage
- This function works as follows:
    - Read in the input data from the folder
    - Initialize a new crossbar
    - Set the weights in the crossbar based on the input weights
    - Define a voltage pulse based on the settings in `simulation_settings.h` (by default: 50 us pulse with 5 us rise and fall time, with 0.1V height and a time step of 1 us)
    - For each (of the 10) input voltages:
        - For each time step:
            - Use the `crossbar.ApplyVoltage()` function
            - Collect the output currents for each individual memristor
    - Collect the average individual memristor currents
    - Calculate the averaged output crossbar currents
    - Write the individual memristor curernts and crossbar currents to their appropriate files
- `RRAM_validation.bat` is an example script to apply this function to multiple folders

Various simulation related settings can be altered in `simulation_settings.h`
