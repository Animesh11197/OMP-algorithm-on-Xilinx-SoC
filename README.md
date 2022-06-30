# OMP-algorithm
This folder contains all the scipts and codes to implement OMP algorithm, both on MATLAB as well as on FPGA using High Level Synthesis.

CONTENTS:
1. OMP_MATLAB.m - contains the MATLAB script to implement OMP algorithm on MATLAB.
2. MAT_to_DAT.m - contains the MATLAB scipt for converting .mat datasets of different SNRs to create dataset files for our C code to read from.
3. OMP_source_codes - contains the source code for implemtation of OMP on FPGA.


STEPS for running OMP code on FPGA would be :
1. load the .MAT datasets into MATLAB workspace and then run the MAT_to_DAT.m script, this will create the .DAT data files in the same working directory of MATLAB.
2. You will need to move these .DAT files along to the FPGA.
3. Carefully set the target SoC and now build the source code files in Xilinx SDSoC and move the generated .elf file, BOOT.BIN file and the image.ub files into the FPGA.
4. After this we can reboot and run the specific executable file on Board.
