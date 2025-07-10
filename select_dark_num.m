function [] = select_dark_num(regionName)


%The name of the region
params.region = regionName;
%params.region = 'Luxembourg';


%Path of the data file
dataFile = ['./data/' params.region '.xlsx'];

%Population
if regionName == "wwtp1"
    params.N = 42931;
elseif regionName == "wwtp2"
    params.N = 68070;
elseif regionName == "wwtp3"
    params.N = 110871;
elseif regionName == "wwtp4"
    params.N = 60262;
else
    params.N = 240000;
end

%The average ratio of total and detected cases. If there are considerable
%jumps expected, different values can be used for different time periods.
%In that case, the darkNumber should be a M x 2 matrix, where the first
%column contains the different ratios, and second column the start day of
%the corresponding value (given as the number of the day counted from the
%beginning of the data).

J_min = inf;
dn_min = 0;
all_J = [];
all_nu = [];
all_gamma = [];
all_eps = [];
count = 1;

n_start = 1.1;
n_inc = 0.1;
n_end = 4;
for dn = n_start:n_inc:n_end
    params.darkNumber = [dn 1];

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
    datesRaw = TT{:,1};  % First column: dates
    YC = TT.cases';
    YW = TT.ww';
    YWip = WWinterpol(YW);

    % Convert to datetime if needed
    if ~isa(datesRaw, 'datetime')
        dates = datetime(datesRaw);
    else
        dates = datesRaw;
    end

    % Find the index of 01-Mar-2022
    targetDate = datetime('02-Aug-2021', 'InputFormat', 'dd-MMM-yyyy');
    index_calib = find(dates == targetDate);

    if isempty(index_calib)
        error('La data 01-Mar-2022 non è presente nel file.');
    end

    % Calibrate the model and save parameters
    YC_calibrate = YC(1:index_calib);
    YW_calibrate = YWip(1:index_calib);
    % C_calibrate = C(:,1:index_calib);

    %Determine c_t and plot label dates
    [C, labs, firsts, longDates] = SEIRWWinit(YC_calibrate,startDate,specialHolidays,params.darkNumber);


    [params,J] = SEIRWWcalibrate(YC_calibrate,YW_calibrate,C,params);
    save(['./parameters/params_' params.region '.mat'],'params','specialHolidays','startDate','labs','firsts')

    all_J = [all_J; J];
    if J < J_min
        dn_min = dn;
        J_min = J;
    end
    all_nu = [all_nu; params.nu];
    all_eps = [all_eps; params.WWexp];
    all_gamma = [all_gamma; params.gamma];

    fprintf('iterazione numero: %d\n', count);
    count = count +1;
end

disp("Il migliore dark number è: ")
disp(dn_min)

%% GRAFICI

%% Grafico J
figure('Position', [100, 200, 1200, 450]);
hold on; grid on;
plot(all_J, 'r-', 'LineWidth', 2);
set(gca, 'FontSize', 14, 'Layer', 'top');
ylabel('Cost function', 'FontSize', 16);
xlabel('Dark number', 'FontSize', 16)
x_vals = n_start:n_inc:n_end;
xticks(1:length(x_vals));
xticklabels(arrayfun(@(x) sprintf('%.1f', x), x_vals, 'UniformOutput', false));
if ~exist('img', 'dir')
    mkdir('img');
end
outputFolder = fullfile('img', regionName);
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end
fileName = fullfile(outputFolder, ['img_J_' regionName '.tex']);
matlab2tikz(fileName, ...
    'showInfo', false);

%% Grafico nu
figure('Position', [100, 200, 1200, 450]);
hold on; grid on;
plot(all_nu, 'r-', 'LineWidth', 2);
set(gca, 'FontSize', 14, 'Layer', 'top');
ylabel('\nu', 'FontSize', 16);
xlabel('Dark number', 'FontSize', 16)
x_vals = n_start:n_inc:n_end;
xticks(1:length(x_vals));
xticklabels(arrayfun(@(x) sprintf('%.1f', x), x_vals, 'UniformOutput', false));
if ~exist('img', 'dir')
    mkdir('img');
end
outputFolder = fullfile('img', regionName);
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end
fileName = fullfile(outputFolder, ['img_nu_' regionName '.tex']);
matlab2tikz(fileName, ...
    'showInfo', false);

%% Grafico epsilon
figure('Position', [100, 200, 1200, 450]);
hold on; grid on;
plot(all_eps, 'r-', 'LineWidth', 2);
set(gca, 'FontSize', 14, 'Layer', 'top');
ylabel('\epsilon', 'FontSize', 16);
xlabel('Dark number', 'FontSize', 16)
x_vals = n_start:n_inc:n_end;
xticks(1:length(x_vals));
xticklabels(arrayfun(@(x) sprintf('%.1f', x), x_vals, 'UniformOutput', false));
if ~exist('img', 'dir')
    mkdir('img');
end
outputFolder = fullfile('img', regionName);
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end
fileName = fullfile(outputFolder, ['img_epsilon_' regionName '.tex']);
matlab2tikz(fileName, ...
    'showInfo', false);

%% Grafico gamma
figure('Position', [100, 200, 1200, 450]);
hold on; grid on;
plot(all_gamma, 'r-', 'LineWidth', 2);
set(gca, 'FontSize', 14, 'Layer', 'top');
ylabel('\gamma', 'FontSize', 16);
xlabel('Dark number', 'FontSize', 16)
x_vals = n_start:n_inc:n_end;
xticks(1:length(x_vals));
xticklabels(arrayfun(@(x) sprintf('%.1f', x), x_vals, 'UniformOutput', false));
if ~exist('img', 'dir')
    mkdir('img');
end
outputFolder = fullfile('img', regionName);
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end
fileName = fullfile(outputFolder, ['img_gamma_' regionName '.tex']);
matlab2tikz(fileName, ...
    'showInfo', false);


end