clear
close all
clc

basePath = getenv('USERPROFILE');  
addpath(genpath(fullfile(basePath, 'MATLAB', 'matlab2tikz')));

% Choose if you want do the sensitivity analysis
sensitivity_analysis = false;

% Choose which analysis you want to do
sens_analysis_R_matrix = true;
sens_analysis_dark_number = false;

% Choose which region you want study
% regionName = 'Luxembourg';
regionName = 'wwtp1';
% regionName = 'wwtp2';
% regionName = 'wwtp3';
% regionName = 'wwtp4';

paramFile = fullfile('parameters', ['params_' regionName '.mat']);

if sensitivity_analysis
    sens_analysis(regionName,sens_analysis_R_matrix,sens_analysis_dark_number)
end

% Choose the dark number 
if regionName == 'wwtp1'
    new_dark_number = 1.5;
elseif regionName == 'wwtp2'
    new_dark_number = 1.3;
elseif regionName == 'wwtp3'
    new_dark_number = 1.3;
elseif regionName == 'wwtp4'
    new_dark_number = 1.4;
else
    new_dark_number = 1.8;
end



if ~exist(paramFile, 'file') 
    latest_dark_number = new_dark_number;
else
    load(['./parameters/params_' regionName '.mat'])
    latest_dark_number = params.darkNumber(1);
end

%% Run the setup 
if ~exist(paramFile, 'file') || latest_dark_number ~= new_dark_number
    fprintf('File %s non trovato. Lancio setup...\n', paramFile);
    setup(regionName,new_dark_number);
else
    fprintf('File %s già presente. Setup non necessario.\n', paramFile);
end

%% Run the SEIR-WW-EKF. Calibration is not required every time new data is available.

%Path of the data file
dataFile = ['./data/' regionName '.xlsx'];

load(['./parameters/params_' regionName '.mat'])
addpath('./SEIRWWfiles/')
TT = readtable(dataFile);
date = datetime(TT.date);
YC = TT.cases';
YW = TT.ww';
YWip = WWinterpol(YW);

%Determine c_t and plot label dates. Remember to update the specialHolidays
%indices if needed.

[C, labs, firsts, longDates] = SEIRWWinit(YC, startDate, specialHolidays, params.darkNumber);
params.RW = params.RW0/10;
%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% 
%
%   SECTIONS BELOW CAN BE RUN ONE-BY-ONE
%
%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% 

RC_scale = 1e0;
RW_scale = 1e0;

% %% Estimate the wastewater data using only the case data
[Yest_case, ~, ~, ReffCase, ~] = SEIR_WW(params,YC,YW,C,[true,true],1000,firsts,labs,true,regionName,RW_scale,RC_scale);
% 
% %% Estimate the case numbers using only the wastewater data
[Yest_ww, ~, ~, ReffWW, ~] = SEIR_WW(params,YC,YW,C,[false,true],1000,firsts,labs,true,regionName,RW_scale,RC_scale);

%% Estimate the case numbers using both data
[Yest_both, ~, ~, Reff_both, ~] = SEIR_WW(params,YC,YW,C,[true,true],1000,firsts,labs,true,regionName,RW_scale,RC_scale);

%%
% Creazione figura
figure('Position', [100, 200, 1200, 450]);
hold on; grid on;
plot(ReffCase, ':k', 'LineWidth', 2);        % Reff da casi clinici (linea tratteggiata nera)
plot(ReffWW, 'r', 'LineWidth', 2);           % Reff da acque reflue (linea continua rossa)
plot([1, length(ReffCase)+10], [1, 1], 'k', 'LineWidth', 1);
if any(strcmp(params.region, {'Luxembourg'}))
    xlabel('Dates 2020–2021', 'FontSize', 16);
else
    xlabel('Dates 2020–2023', 'FontSize', 16);
end
ylabel('$R_{\mathrm{eff}}$', 'Interpreter','latex', 'FontSize', 16);
set(gca, 'FontSize', 14);
xticks(firsts);
xticklabels(labs);
legend({'Case data', 'Wastewater data'}, 'Location', 'northeast', 'FontSize', 14);

if ~exist('img', 'dir')
    mkdir('img');
end
outputFolder = fullfile('img', regionName);
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end
fileName = fullfile(outputFolder, ['grafico_Reff_' regionName '.tex']);
matlab2tikz(fileName, ...
    'showInfo', false);

correlation_Reff = corr(ReffCase', ReffWW')

%% Project forward in time

%At what date is the forward prediction done
predDay = length(YC);
numEstDays = 300;

%Use case/WW data or both
dataToUse = [true, false];

params.RW = params.RW0/10;
[~, XendC, PC] = SEIR_WW(params,YC,YW,C,dataToUse,predDay,firsts,labs,false,regionName);
[Y0, err] = SEIR_WW_FWD(XendC,C,PC,predDay+1,params,numEstDays);
Yc0 = cumsum([YC(1:predDay) Y0(1,:)]);

%% Creazione della figure
figure('Position',[100, 200, 1200, 450]); 
hold on; grid on;
plot(movmean([YC(1:predDay), Y0(1,:)], [6, 0]), 'r-', 'LineWidth', 2);
plot(movmean(YC, [6, 0]), 'k-', 'LineWidth', 2);
if any(strcmp(params.region, {'Luxembourg'}))
    xlabel('Dates 2020–2021', 'FontSize', 16);
else
    xlabel('Dates 2020–2023', 'FontSize', 16);
end
ylabel('Daily cases', 'FontSize', 16);
set(gca, 'FontSize', 14);
legend({'Projection', 'Data'}, 'Location', 'northeast', 'FontSize', 14);
box on;
if ~exist('img', 'dir')
    mkdir('img');
end
outputFolder = fullfile('img', regionName);
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end
fileName = fullfile(outputFolder, ['grafico_casi_giornalieri_' regionName '.tex']);
matlab2tikz(fileName, ...
    'showInfo', false);


figure('Position', [100, 200, 1200, 450]); 
hold on; grid on;
plot(Yc0, 'r-', 'LineWidth', 2);
plot(cumsum(YC), 'k-', 'LineWidth', 2);
if any(strcmp(params.region, {'Luxembourg'}))
    xlabel('Days starting from 25-Feb-2020', 'FontSize', 16);
else
    xlabel('Days starting from 04-Oct-2021', 'FontSize', 16);
end
ylabel('Cumulative cases', 'FontSize', 16);
set(gca, 'FontSize', 14);
set(gca, 'Layer', 'top');
legend({'Projection', 'Data', '±2σ Range', '±1σ Range'}, 'Location', 'northwest', 'FontSize', 14);
if ~exist('img', 'dir')
    mkdir('img');
end
outputFolder = fullfile('img', regionName);
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end
fileName = fullfile(outputFolder, ['grafico_casi_cumulativi_' regionName '.tex']);
matlab2tikz(fileName, ...
    'showInfo', false);


%% Correlation between cases and ww data
% Normalizza entrambi i vettori [0, 1]
cases_norm = (YC - min(YC)) / (max(YC) - min(YC));
ww_norm = (YW - min(YW)) / (max(YW) - min(YW));

% Calcola la correlazione di Pearson
R = corr(cases_norm, ww_norm);
% Plot
figure('Position', [100, 200, 1200, 450]); 
hold on; grid on;
plot(date, cases_norm, 'r-', 'LineWidth', 2);
plot(date, ww_norm, 'k-', 'LineWidth', 2);
xlabel('Date');
ylabel('Normalized value', 'FontSize', 16);
set(gca, 'FontSize', 14, 'Layer', 'top');
legend({'Normalized Cases', 'Normalized Wastewater'}, 'Location', 'best','FontSize', 14);
if ~exist('img', 'dir')
    mkdir('img');
end
outputFolder = fullfile('img', regionName);
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end
fileName = fullfile(outputFolder, ['img_data_' regionName '.tex']);
matlab2tikz(fileName, ...
    'showInfo', false);
