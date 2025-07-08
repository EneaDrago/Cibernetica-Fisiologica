clear
close all
clc

%Load the data and parameters
% regionName = 'Luxembourg';
regionName = 'wwtp1';
% regionName = 'wwtp2';
% regionName = 'wwtp3';
% regionName = 'wwtp4';
% regionName = 'wwtpTOT';

paramFile = fullfile('parameters', ['params_' regionName '.mat']);

if ~exist(paramFile, 'file')
    fprintf('File %s non trovato. Lancio setup...\n', paramFile);
    setup(regionName);
else
    fprintf('File %s già presente. Setup non necessario.\n', paramFile);
end


%% Run the SEIR-WW-EKF. Calibration is not required every time new data is available.
addpath(genpath('C:\Users\OEM\MATLAB\matlab2tikz'))

%Path of the data file
dataFile = ['./data/' regionName '.xlsx'];

load(['./parameters/params_' regionName '.mat'])
addpath('./SEIRWWfiles/')
TT = readtable(dataFile);
YC = TT.cases';
YW = TT.ww';

%Determine c_t and plot label dates. Remember to update the specialHolidays
%indices if needed.
[C, labs, firsts, longDates] = SEIRWWinit(YC, startDate, specialHolidays, params.darkNumber);

%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% 
%
%   SECTIONS BELOW CAN BE RUN ONE-BY-ONE
%
%%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% 
%% Estimate the wastewater data using only the case data
params.RW = params.RW0/10;
[Yest_case, ~, ~, ReffCase, ~] = SEIR_WW(params,YC,YW,C,[true,false],1000,firsts,labs,true,regionName);

%% Estimate the case numbers using only the wastewater data
[Yest_ww, ~, ~, ReffWW, ~] = SEIR_WW(params,YC,YW,C,[false,true],1000,firsts,labs,true,regionName);

%% Estimate the case numbers using only the interpolated wastewater data
% YWip = WWinterpol(YW);
% [Yest_wwip, ~, ~, ~, ~] = SEIR_WW(params,YC,YWip,C,[false,true],1000,firsts,labs,true,regionName);


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
numEstDays = 330;

%Use case/WW data or both
dataToUse = [true, false];

params.RW = params.RW0/10;
[~, XendC, PC] = SEIR_WW(params,YC,YW,C,dataToUse,predDay,firsts,labs,false,regionName);
[Y0, err] = SEIR_WW_FWD(XendC,C,PC,predDay+1,params,numEstDays);
Yc0 = cumsum([YC(1:predDay) Y0(1,:)]);
Xaux = XendC;
Xaux(7) = XendC(7) + PC(7,7).^.5;
Y1 = SEIR_WW_FWD(Xaux,C,PC,predDay+1,params,330);
Yc1 = cumsum([YC(1:predDay) Y1(1,:)]);
Xaux(7) = XendC(7) + 2*PC(7,7).^.5;
Y2 = SEIR_WW_FWD(Xaux,C,PC,predDay+1,params,330);
Yc2 = cumsum([YC(1:predDay) Y2(1,:)]);
Xaux(7) = XendC(7) - PC(7,7).^.5;
Ym1 = SEIR_WW_FWD(Xaux,C,PC,predDay+1,params,330);
Ycm1 = cumsum([YC(1:predDay) Ym1(1,:)]);
Xaux(7) = XendC(7) - 2*PC(7,7).^.5;
Ym2 = SEIR_WW_FWD(Xaux,C,PC,predDay+1,params,330);
Ycm2 = cumsum([YC(1:predDay) Ym2(1,:)]);

%% Creazione della figure
figure('Position',[400, 200, 560, 380]); 
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
xticks(firsts+firsts);
xticklabels(labs);
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

figure('Position', [400, 200, 560, 380]); 
hold on; grid on;
plot(Yc0, 'r-', 'LineWidth', 2);
plot(cumsum(YC), 'k-', 'LineWidth', 2);
h1 = fill([min(predDay,length(YC)):length(Yc0), fliplr(predDay:length(Yc0))], ...
          [Yc0(predDay), Ycm2(predDay+1:end), fliplr(Yc2(predDay+1:end)), Yc0(predDay)], ...
          [1, 0.77, 0.77], 'EdgeColor', 'none', 'FaceAlpha', 0.3);

h2 = fill([min(predDay,length(YC)):length(Yc0), fliplr(predDay:length(Yc0))], ...
          [Yc0(predDay), Ycm1(predDay+1:end), fliplr(Yc1(predDay+1:end)), Yc0(predDay)], ...
          [1, 0.67, 0.67], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
if any(strcmp(params.region, {'Luxembourg'}))
    xlabel('Dates 2020–2021', 'FontSize', 16);
else
    xlabel('Dates 2020–2023', 'FontSize', 16);
end
ylabel('Cumulative cases', 'FontSize', 16);
set(gca, 'FontSize', 14);
set(gca, 'Layer', 'top');
xticks(firsts+firsts);
xticklabels(labs);
legend({'Projection', 'Data', '±2σ Range', '±1σ Range'}, 'Location', 'northwest', 'FontSize', 14);
if ~exist('immagini', 'dir')
    mkdir('immagini');
end

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


%% Project forward in time predizione prima

%At what date is the forward prediction done
numEstDays = 330;
predDay = length(YC)-numEstDays+1;
numEstDays = 330;

%Use case/WW data or both
dataToUse = [true, false];

params.RW = params.RW0/10;
[~, XendC, PC] = SEIR_WW(params,YC,YW,C,dataToUse,predDay,firsts,labs,false,regionName);
[Y0, err] = SEIR_WW_FWD(XendC,C,PC,predDay+1,params,numEstDays);
Yc0 = cumsum([YC(1:predDay) Y0(1,:)]);
Xaux = XendC;
Xaux(7) = XendC(7) + PC(7,7).^.5;
Y1 = SEIR_WW_FWD(Xaux,C,PC,predDay+1,params,330);
Yc1 = cumsum([YC(1:predDay) Y1(1,:)]);
Xaux(7) = XendC(7) + 2*PC(7,7).^.5;
Y2 = SEIR_WW_FWD(Xaux,C,PC,predDay+1,params,330);
Yc2 = cumsum([YC(1:predDay) Y2(1,:)]);
Xaux(7) = XendC(7) - PC(7,7).^.5;
Ym1 = SEIR_WW_FWD(Xaux,C,PC,predDay+1,params,330);
Ycm1 = cumsum([YC(1:predDay) Ym1(1,:)]);
Xaux(7) = XendC(7) - 2*PC(7,7).^.5;
Ym2 = SEIR_WW_FWD(Xaux,C,PC,predDay+1,params,330);
Ycm2 = cumsum([YC(1:predDay) Ym2(1,:)]);

%% Creazione della figure 
figure('Position',[400, 200, 560, 380]); 
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
xticks(firsts+firsts);
xticklabels(labs);
legend({'Projection', 'Data'}, 'Location', 'northeast', 'FontSize', 14);
box on;

if ~exist('img', 'dir')
    mkdir('img');
end
outputFolder = fullfile('img', regionName);
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end
fileName = fullfile(outputFolder, ['projection_validation' regionName '.tex']);
matlab2tikz(fileName, ...
    'showInfo', false);

figure('Position', [400, 200, 560, 380]); 
hold on; grid on;
plot(Yc0, 'r-', 'LineWidth', 2);
plot(cumsum(YC), 'k-', 'LineWidth', 2);
h1 = fill([min(predDay,length(YC)):length(Yc0), fliplr(predDay:length(Yc0))], ...
          [Yc0(predDay), Ycm2(predDay+1:end), fliplr(Yc2(predDay+1:end)), Yc0(predDay)], ...
          [1, 0.77, 0.77], 'EdgeColor', 'none', 'FaceAlpha', 0.3);

h2 = fill([min(predDay,length(YC)):length(Yc0), fliplr(predDay:length(Yc0))], ...
          [Yc0(predDay), Ycm1(predDay+1:end), fliplr(Yc1(predDay+1:end)), Yc0(predDay)], ...
          [1, 0.67, 0.67], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
if any(strcmp(params.region, {'Luxembourg'}))
    xlabel('Dates 2020–2021', 'FontSize', 16);
else
    xlabel('Dates 2020–2023', 'FontSize', 16);
end
ylabel('Cumulative cases', 'FontSize', 16);
set(gca, 'FontSize', 14);
set(gca, 'Layer', 'top');
xticks(firsts+firsts);
xticklabels(labs);
legend({'Projection', 'Data', '±2σ Range', '±1σ Range'}, 'Location', 'northwest', 'FontSize', 14);
if ~exist('immagini', 'dir')
    mkdir('immagini');
end

if ~exist('img', 'dir')
    mkdir('img');
end
outputFolder = fullfile('img', regionName);
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end
fileName = fullfile(outputFolder, ['projection_validation_cumsum' regionName '.tex']);
matlab2tikz(fileName, ...
    'showInfo', false);