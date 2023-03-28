%Highway Traffic Dynamics: Data-Driven Analysis and Forecast 
%Allan M. Avila 2019
%University of California Santa Barbara

%% Instructions for Koopman Modes.
% To generate Koopman modes call on the function:
% Modes=GenerateKoopmanModes(Data,Mode1,Mode2,Save)

%Inputs Required:
% Data is a string containing the name of the data set.

% Mode1 and Mode2 are integers indicating which modes to produce.
% Ordered by their period of oscilaliton from slowest to fastest.
% Mode1 can be < or = Mode2.

% Save is a logical (0 or 1) indicating the modes to be saved as jpeg's. 
% Except for the multliLane modes which are videos.

% Examples:2
clc; clear variables; close all;
Poi_Daily=GenerateKoopmanModes('Poi_Day_mean',1,3,0);







