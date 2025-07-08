function [] = setup(regionName)
%% Set initial parameters

%The name of the region
params.region = regionName;
%params.region = 'Luxembourg';


%Path of the data file
dataFile = ['./data/' params.region '.xlsx'];

%Population
params.N = 240000;

%The average ratio of total and detected cases. If there are considerable
%jumps expected, different values can be used for different time periods.
%In that case, the darkNumber should be a M x 2 matrix, where the first
%column contains the different ratios, and second column the start day of
%the corresponding value (given as the number of the day counted from the
%beginning of the data).
params.darkNumber = [1.8 1];

%Indices of special holidays with reduced testing resulting in lower than
%expected case numbers
specialHolidays = [];

%Start date of the data (used for plots), format 'DD/MM/YYYY'
if params.region == "Luxembourg"
    startDate = '25/02/2020';
else
    startDate = '04/10/2021';
end


%% Import data, calibrate model, and save parameters and data together

addpath('./SEIRWWfiles/')

%Default data file name is myRegion.
TT = readtable(dataFile);
YC = TT.cases';
YW = TT.ww';

%Determine c_t and plot label dates
[C, labs, firsts, longDates] = SEIRWWinit(YC,startDate,specialHolidays,params.darkNumber);

% Calibrate the model and save parameters
params = SEIRWWcalibrate(YC,YW,C,params);
save(['./parameters/params_' params.region '.mat'],'params','specialHolidays','startDate','labs','firsts')







