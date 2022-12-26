# MDCCRN
## WHAT IS MDCCRN
mdccrn is short for monoticity dependencies calculation in chemical reaction networks (CRNs).
It is a program that aims to improve accuracy for e.g. the metering of medication by finding transformations
of the classical coordinates of a crn that show monotonic behaviour.


## INSTALLATION INSTRUCTIONS
- Installation has only been tested on Ubuntu 16.04 and Slackware 14.2. It should however run on all major Linux distributions.
A thorough instruction to building the gui version as well as Windows support may be added in the future.
- Download the mdccrn.tar.gz or the mdccrn.zip file
- Extract the downloaded file (e.g. with `tar -xvzf mdccrn.tar.gz`)
- Go into the extracted dir 'mdccrn' (`cd mdccrn`)
- To test the program, go through the usage instructions below and then run it once with one of the provided example input files.


## (OPTIONAL) RECOMPILE FROM SOURCE
- The gui version can not yet be automatically compiled from source. The following therefore only creates the command line version of the program.
- Download and configure the GiNaC-library and the CLN-library (which GiNaC depends on) [https://www.ginac.de/Download.html].
- After downloading and extracting the mdccrn.tar.gz file (as described in the installation instructions above),
go into the extracted dir and use the make command (`make`).


## HOW TO USE MDCCRN
- Create an input file containing all necessary information of the chemical reaction network.
It must have the same structure as the provided example input files.
- Choose an output directory. All output files will be written to this directory. As these can be a lot of files, it is recommended to
create a dedicated directory for every chemical reaction network. This is however not necessary as long as an unique model name is chosen.
- Run the program; the executables are located in the bin folder. You may use the gui version or execute the program from within the terminal.
- The mode specifies the level of "thoroughness" (and performance-intensitivity) that the monotonicity is checked with.
Currently implemented are mode 0, which checks for monoticity among the species and reaction rates, and mode 1,
which ignores the reaction rates for performance gain. The default mode is 0. \
The syntax for the latter is as follows (from within the bin/linux/ folder): \
`./mdccrn_command_line <input file> <output directory> <model name> <mode (optional)>` \
Usage Example: \
`./mdccrn_command_line ../../input/input_example_1.txt ../../output model1`
- Note that the computational cost for calculating every graph with the required monoticity characteristics rises exponentially
when increasing the amount of species or reactions within the CRN. See below for details on performance.


## PERFORMANCE
The following data will give you an impression of how long the calculations may take.
It springs from testing various input examples on my laptop with 4GB RAM and an AMD A8-6410 APU.
The UNIX time command was used to determine execution time.

| name   | number reaction rates | number species | number reactions | number transformations checked | execution time | number consistent graphs |
| :----- | :-------------------: | -------------: | ---------------: | -----------------------------: | -------------: | -----------------------: |
| input1 | 2                     | 3              | 1                | 34                             | 0.02s          | 7                        |
| input2 | 4                     | 4              | 2                | 1322                           | 0.51s          | 45                       |
| input3 | 8                     | 5              | 5                | 166100                         | 1min and 27s   | 46                       |


## HOW TO INTERPRET THE OUTPUT
The output comprises of a .dot graph file and a .txt data file for each transformation fulfilling the system monoticity conditions. \
//TODO: detailled description


## IMPLEMENTATION OVERVIEW
The implementation provides:
- a crn class reading all necessary information from a given input file and providing it in form of processible data structures
- a class that takes a crn object and does all necessary calculation steps as well as a mechanism to write results to output files
For details on implementation, view the doxygen documentation and the source code.
