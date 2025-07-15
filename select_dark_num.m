function [] = select_dark_num(regionName,sens_analysis_R_matrix,sens_analysis_dark_number)

%The name of the region
params.region = regionName;

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

if sens_analysis_dark_number

    figHandles = gobjects(1,6);
    for i = 1:6
        figHandles(i) = figure(i);
        clf;
        hold on;
        set(figHandles(i), 'Position', [100, 200, 1200, 450]);  % Imposta dimensione e posizione
    end

    n_start = 1.1; % First dark number
    n_inc = 0.1;
    n_end = 4;     % Last dark number
    
    for dn = n_start:n_inc:n_end
        params.darkNumber = [dn 1];
    
        fprintf('iterazione numero: %d\n', count);
        count = count +1;
    
        %Indices of special holidays with reduced testing resulting in lower than
        %expected case numbers
        specialHolidays = [];
    
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
    
        % Start date of the data (used for plots), format 'DD/MM/YYYY'
        if params.region == "Luxembourg"
            startDate = '25/02/2020';
            targetDate = datetime('02-Aug-2021', 'InputFormat', 'dd-MMM-yyyy');
            index_calib = find(dates == targetDate);
            if isempty(index_calib)
                error('La data 02-Aug-2021 non è presente nel file.');
            end
        else
            startDate = '04/10/2021';
            targetDate = datetime('01-Mar-2022', 'InputFormat', 'dd-MMM-yyyy');
            index_calib = find(dates == targetDate);
            if isempty(index_calib)
                error('La data 01-Mar-2022 non è presente nel file.');
            end
        end
    
        % Calibrate the model and save parameters
        YC_calibrate = YC(1:index_calib);
        YW_calibrate = YWip(1:index_calib);
    
        %Determine c_t and plot label dates
        [C, labs, firsts, longDates] = SEIRWWinit(YC_calibrate,startDate,specialHolidays,params.darkNumber);
    
    
        [params,J] = SEIRWWcalibrate(YC_calibrate,YW_calibrate,C,params);
    
        all_J = [all_J; J];
        if J < J_min
            dn_min = dn;
            J_min = J;
        end
        all_nu = [all_nu; params.nu];
        all_eps = [all_eps; params.WWexp];
        all_gamma = [all_gamma; params.gamma];
    
        if abs(dn - 1.8) < 1e-6
            params.RW = params.RW0/10;
            plot_color = 'g';
            [~, ~, ~, ~, ~] = SEIR_WW_sens(params,YC,YW,C,[true,false],1000,firsts,labs,true,regionName, figHandles,plot_color);
            [~, ~, ~, ~, ~] = SEIR_WW_sens(params,YC,YW,C,[false,true],1000,firsts,labs,true,regionName, figHandles,plot_color);
            disp(1.8)
        end
    
        if dn==3
            params.RW = params.RW0/10;
            plot_color = 'r';
            [~, ~, ~, ~, ~] = SEIR_WW_sens(params,YC,YW,C,[true,false],1000,firsts,labs,true,regionName, figHandles,plot_color);
            [~, ~, ~, ~, ~] = SEIR_WW_sens(params,YC,YW,C,[false,true],1000,firsts,labs,true,regionName, figHandles,plot_color);
            disp(3)
        end
    
    end
    
    params.darkNumber = [dn_min 1];
    [C, labs, firsts, ~] = SEIRWWinit(YC_calibrate,startDate,specialHolidays,params.darkNumber);
    [params,~] = SEIRWWcalibrate(YC_calibrate,YW_calibrate,C,params);
    plot_color = 'b';
    [~, ~, ~, ~, ~] = SEIR_WW_sens(params,YC,YW,C,[true,false],1000,firsts,labs,true,regionName, figHandles,plot_color);
    [~, ~, ~, ~, ~] = SEIR_WW_sens(params,YC,YW,C,[false,true],1000,firsts,labs,true,regionName, figHandles,plot_color);
    
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

if sens_analysis_R_matrix

    figHandles = gobjects(1,6);
    for i = 1:6
        figHandles(i) = figure(i);
        clf;
        hold on;
        set(figHandles(i), 'Position', [100, 200, 1200, 450]);  % Imposta dimensione e posizione
    end

    allHandles1 = gobjects(0); allLabels1 = {};
    allHandles2 = gobjects(0); allLabels2 = {};
    allHandles3 = gobjects(0); allLabels3 = {};

    RW_start = -3;
    RW_inc = 2;
    RW_end = 3;    

    RC_start = -3;
    RC_inc = 2;
    RC_end = 3;     
    
    for RW_scale = RW_start:RW_inc:RW_end
        for RC_scale = RC_start:RC_inc:RC_end
            params.darkNumber = [1.5 1];
        
            fprintf('iterazione numero: %d\n', count);
            count = count +1;
        
            %Indices of special holidays with reduced testing resulting in lower than
            %expected case numbers
            specialHolidays = [];
        
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
        
            % Start date of the data (used for plots), format 'DD/MM/YYYY'
            if params.region == "Luxembourg"
                startDate = '25/02/2020';
                targetDate = datetime('02-Aug-2021', 'InputFormat', 'dd-MMM-yyyy');
                index_calib = find(dates == targetDate);
                if isempty(index_calib)
                    error('La data 02-Aug-2021 non è presente nel file.');
                end
            else
                startDate = '04/10/2021';
                targetDate = datetime('01-Mar-2022', 'InputFormat', 'dd-MMM-yyyy');
                index_calib = find(dates == targetDate);
                if isempty(index_calib)
                    error('La data 01-Mar-2022 non è presente nel file.');
                end
            end
        
            % Calibrate the model and save parameters
            YC_calibrate = YC(1:index_calib);
            YW_calibrate = YWip(1:index_calib);
        
            %Determine c_t and plot label dates
            [C, labs, firsts, longDates] = SEIRWWinit(YC_calibrate,startDate,specialHolidays,params.darkNumber);
        
            [params,J] = SEIRWWcalibrate(YC_calibrate,YW_calibrate,C,params);
       
            params.RW = params.RW0/10;
            [~, ~, ~, ~, ~, ~, h1, l1, h2, l2, h3, l3] = SEIR_WW_sens_R(params,YC,YW,C,[true,true],1000,firsts,labs,true,regionName, figHandles,10^RW_scale,10^RC_scale);

            allHandles1(end+1) = h1; allLabels1{end+1} = l1;
            allHandles2(end+1) = h2; allLabels2{end+1} = l2;
            allHandles3(end+1) = h3; allLabels3{end+1} = l3;

        end
    
    end
    figure(figHandles(1));
    legend(allHandles1, allLabels1, 'Location','northeast', 'FontSize', 14);
    
    figure(figHandles(2));
    legend(allHandles2, allLabels2, 'Location','northeast', 'FontSize', 14);
    
    figure(figHandles(3));
    legend(allHandles3, allLabels3, 'Location','northwest', 'FontSize', 14);

end



end