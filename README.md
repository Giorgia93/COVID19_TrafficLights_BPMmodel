# COVID19_TrafficLights_BPMmodel

Branching process reinfection model to simulate Delta COVID19 outbreak in New Zealand. Coded in MATLAB2019b.

To run:

1. Dowload all files and folders
2. (optional) Open mainSimLeakyParallel.m and change parameter values as needed
3. Run mainSimLeakyParallel.m

Main files:

- getParParallel.m Most model parameters are defined here.
- mainSimLeakyParallel.m Main run file. Parameters subject to change for sensitivity analyses are defined here
- runSimLeaky_LessVerbose.m Main loop of the branching process model
The main files call on spreadsheets in the "data" folder and on Matlab dependencies in the "dependencies folder". Running mainSimLeakyParallel.m will produce timeseries in the "results/timeseries" folder, a summary spreadsheet in the "results/summary" folder.
