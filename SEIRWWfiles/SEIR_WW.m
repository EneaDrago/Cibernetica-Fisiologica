function [Yest, Xend, P, Reff, errReff, Ysd] = SEIR_WW(params,YC,YW,C,useData,maxind,firsts,labs,pl,regionName)

%params: the coefficients in the SEIR ODE
%useData is 2D logical: 
%   useData(1): Use the daily case data?
%   useData(2): Use the wastewater data?


%%
if ~exist('pl')
    pl = false;
end

if length(maxind) == 1 && maxind > length(YC)
    maxind = length(YC);
end

%Indicator for plotting big figures
plotSpec = false;
if ischar(pl)
    plotSpec = true;
    if strcmp(pl,'big')
        plotloc = [100,200,1200,450];
    elseif strcmp(pl,'small')
        plotloc = [100,200,560,350];
    end
    pl = true;
end

%Process the WW data
WWinds = find(YW>-.5);
YW(WWinds) = 1e-5*YW(WWinds).^params.WWexp; 
minYW = min(YW(WWinds));
YCaux = sort(YC,'ascend');
excl = sum(YCaux < 0);
YWaux = sort(YW(WWinds),'ascend');
ccc = mean(YCaux(excl+1:excl+floor(length(YC)/10)))/mean(YCaux(excl+floor(length(YC)/10)+1:end));
aaa = mean(YWaux(1:floor(length(YWaux)/10)));
bbb = mean(YWaux(floor(length(YWaux)/10)+1:end));
if ccc*bbb < aaa
    YW = YW - min((aaa-ccc*bbb)/(1-ccc),minYW);
end


% Set parameters
alpha = params.alpha;
beta = params.beta;
tau = params.tau;
gamma = params.gamma;
nu = params.nu;   
eta = 1; % W compartment omitted
CC = params.modelErrorC; 
N = params.N;
params.sigma = 1; %Doesn't matter, W compartment omitted

%For sensitivity analysis
if isfield(params,'ctFactor')
    C = C*params.ctFactor;
    nu = nu*params.ctFactor;
end

%Outlier limit: If the discrepancy of WW data and modeled output is higher
%than OL_limit standard deviations, the Kalman correction is plateaued on
%the limit.
if isfield(params,'outlierLimit')
    OL_limit = params.outlierLimit;
else
    OL_limit = 4;
end


% Reaction stoichiometry (w.r.t. SEIRreaction-function)
% guardando colonna per colonna, mostra come cambiano le variabili di stato
% (è la matrice che, in ISI, chiamavamo F. E' la matrice A in x'=Ax)
AR = [-1  0  0  0  0  0;   %S
       1 -1  0  0  0  0;   %E
       0  1 -1  0  0  0;   %I
       1  0  0 -1  0  0;   %A
       0  1  0  0  0  0;   %N
       0  0  0  0  1 -1;   %W
       0  0  0  0  0  0];  %beta

minWW = 0;
iaux = find(YW(WWinds) < minWW);
YW(WWinds(iaux)) = minWW;

% Time step = 1/N_step (days)
N_step = 10;

% Initial error variance of beta
S_beta = params.S_beta;

% Variance of daily change of beta (initially)
Q_beta = params.Q_beta0;

%Initial state
% X(1): S(t)
% X(2): E(t)
% X(3): I(t)
% X(4): A(t)
% X(5): N(t)
% X(6): W(t)
% X(7): beta(t)
X = zeros(7,length(YC)+1);
X(:,1) = [N-params.E_init-params.I_init; params.E_init; params.I_init; params.E_init; params.E_init; 0; beta];


% Initial state error covariance
P = [params.varE_init + params.varI_init, -params.varE_init, -params.varI_init, -params.varE_init, -params.varE_init, 0, 0;
    -params.varE_init, params.varE_init, 0, params.varE_init, params.varE_init, 0, 0;
    -params.varI_init, 0, params.varI_init, 0, 0, 0, 0;
    -params.varE_init, params.varE_init, 0, params.varE_init, params.varE_init, 0, 0;
    -params.varE_init, params.varE_init, 0, params.varE_init, params.varE_init, 0, 0;
     0, 0, 0, 0, 0, 0, 0;
     0, 0, 0, 0, 0, 0, S_beta];


% Measurement error variance (for cases), assuming a Binomial distribution for the number of cases
RC = movmean(YC,[6,0]).*C(1:length(YC))./movmean(C(1:length(YC)),[6,0]).*(1-C(1:length(YC))) + 1;
if isfield(params,'Radditional')
    RC = RC + params.Radditional;
end


% Measurement error variance for WW data
RW = params.RW;

% Number of detected cases today depends linearly on the true number of new
% cases today
Ccase = [0, 0, 0,  0, 1, 0, 0];
Cww   = [0, 0, 0, nu, 0, 0, 0];

Ypred   = zeros(2,size(YC,2));
Yest    = zeros(2,size(YC,2));
errReff = zeros(1,length(YC));


%% ITERAZIONE KALMAN 
jaux = 1;
for jday = 1:max(maxind)
    
    % Reduce Q_beta after the first month. Higher Q_beta accounts for
    % errors in the initial estimate.
    if jday > 30.5
        Q_beta = params.Q_beta1; 
    end
    
    % Initialise prediction variables
    Xhat = X(:,jday);
    Phat = P;
    
    % Reset the "cases today" counter and corresponding covariance
    Xhat(5) = 0;
    Phat(5,:) = 0;
    Phat(:,5) = 0;
    
    % Time loop for one day (remind: each day is devided into N_step
    % substeps)
    for jt = 1:N_step
        [RR, Jf] = SEIRreaction(Xhat,N,alpha,tau,gamma,nu,eta,1/N_step);
        Xhat = Xhat + AR*RR;
        Q = CC*AR*diag(RR)*AR';
        Q(7,7) = Q_beta/N_step;
        Phat = (eye(7) + AR*Jf)*Phat*(eye(7) + Jf'*AR') + Q;
    end
    
    Call = [];
    R = [];
    outInds = [];
    Yday = [];
    
    % Predicted number of daily new cases and wastewater measurement
    Ypred(:,jday) = [C(jday)*Ccase; Cww]*Xhat + [0; minWW];
    
    %% Considero l'osservazione all'istante di tempo corrente
    % in base a quali misurazioni sono attive (e quali dati sono
    % disponibili), creo una Call

    % Check if there's case data for today
    % useData(1) = True => use case data
    % (if the data is missing, YC(jday) = -1)
    if  useData(1) && YC(jday) > -.5
        Call = [Call; C(jday)*Ccase];
        R = [R; RC(jday)];
        Yday = [Yday; YC(jday)];
        outInds = [outInds; 1];
    end
    
    % Check if there's wastewater data for today
    % useData(2) = True => use ww data
    WWii = 0;
    if useData(2) && YW(jday) > -.00005
        Call = [Call; Cww];
        R = [R; RW];
        Yday = [Yday; YW(jday)];
        outInds = [outInds; 2];
        WWii = length(Yday); 
    end
    
    
    % Check if there was new data on this time step
    if size(Call,1) > 0
    
        % Measurement covariance
        S = Call*Phat*Call' + diag(R);

        % Outlier detection and plateauing
        if WWii > 0
             discrepancy = (Yday(WWii) - Ypred(2,jday))/(S(WWii,WWii) + params.RW0 - params.RW)^.5;
              if abs(discrepancy) > OL_limit
                  Yday(WWii) = Ypred(2,jday) + OL_limit * sign(Yday(WWii)-Ypred(2,jday)) * (S(WWii,WWii) + params.RW0 - params.RW)^.5;
              end
         end
        
        % State update based on true and predicted number
        X(:,jday+1) = Xhat + Phat*Call'*S^-1*(Yday-Ypred(outInds,jday));

        % Covariance update
        P = Phat - Phat*Call'*S^-1*Call*Phat;
        
    else
        
        % In case of no new data, skip the update step
        P = Phat;
        X(:,jday+1) = Xhat;
    end

    % Ensure the states to be non-negative (typically not a problem)
    X(2,jday+1) = max(X(2,jday+1),0);
    X(3,jday+1) = max(X(3,jday+1),0);
    
    % Estimated number of daily new cases and wastewater measurement
    Yest(:,jday) = [C(jday)*Ccase; Cww]*X(:,jday+1) + [0; minWW];

    % Standard deviations for the outputs
    Ysd(1,jday) = C(jday)^2*Ccase*P*Ccase' + RC(jday);
    Ysd(2,jday) = Cww*P*Cww' + params.RW0;
    
    % Store the error variance of beta
    errReff(jday) = P(7,7)^.5*X(1,jday)/N/tau; 
    
    % Store the state estimate of requested times
    if jday == maxind(jaux)
        Xend(:,jaux) = X(:,jday+1);
        jaux = jaux + 1;
    end
end

% Calculate R_eff
Reff = X(7,2:end).*X(1,2:end)/N/tau;

% Plot if requested
if pl

    % se non sono state usate le misure dei CASE (i.e., solo WW)
    if ~useData(1)
       figure('Position', [100, 200, 1200, 450]);
       hold on; grid on;
       Yaux = movmean(Yest(1, :), [6, 0]);             % Stima filtrata dei casi
       Ysdaux = movmean(sqrt(Ysd(1, :)), [6, 0]);      % Deviazione standard smoothed
       fill([1:length(Yaux), fliplr(1:length(Yaux))], ...
           [max(Yaux - 2 * Ysdaux, 0), fliplr(Yaux + 2 * Ysdaux)], ...
           [1, 0.93, 0.93], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
       plot(Yaux, 'r-', 'LineWidth', 2);                          % Stima wastewater
       plot(movmean(YC, [6, 0]), ':k', 'LineWidth', 2);           % Dati reali (smoothed)
       set(gca, 'FontSize', 14, 'Layer', 'top');
       xticks(firsts);
       xticklabels(labs);
       ylabel('Daily new cases', 'FontSize', 16);
       if any(strcmp(params.region, {'Luxembourg'}))
           xlabel('Dates 2020–21', 'FontSize', 16);
       else
           xlabel('Dates 2021–23', 'FontSize', 16);
       end
       legend({'Estimated from wastewater data', 'Data', '2SD envelope'}, ...
           'Location', 'northeast', 'FontSize', 14);
       
       if ~exist('img', 'dir')
           mkdir('img');
       end
       outputFolder = fullfile('img', regionName);
       if ~exist(outputFolder, 'dir')
           mkdir(outputFolder);
       end
       fileName = fullfile(outputFolder, ['grafico_dailycases_ww' regionName '.tex']);
       matlab2tikz(fileName, ...
           'showInfo', false);
        
        % Correlation between case numbers and reconstructed case numbers
        CaseCorr = sum((YC-mean(YC)).*(Yest(1,:)-mean(Yest(1,:))))/norm(YC-mean(YC))/norm(Yest(1,:)-mean(Yest(1,:)))
        
        Ysm1 = movmean(YC,[6 0]);
        Ysm2 = movmean(Yest(1,:),[6 0]);
        
        CaseCorrSm = sum((Ysm1-mean(Ysm1)).*(Ysm2-mean(Ysm2)))/norm(Ysm1-mean(Ysm1))/norm(Ysm2-mean(Ysm2))
        

        % Creazione figura
        figure('Position', [100, 200, 1200, 450]);
        hold on; grid on;
        plot(cumsum(YC), 'k', 'LineWidth', 2);            % Dati reali
        plot(cumsum(Yest(1, :)), 'r', 'LineWidth', 2);     % Stima da acque reflue
        set(gca, 'FontSize', 14, 'Layer', 'top');
        xticks(firsts);
        xticklabels(labs);
        ylabel('Cumulative cases', 'FontSize', 16);
        if any(strcmp(params.region, {'Luxembourg'}))
            xlabel('Dates 2020–21', 'FontSize', 16);
        else
            xlabel('Dates 2021–23', 'FontSize', 16);
        end
        legend({'Data', 'Estimated'}, 'Location', 'northwest', 'FontSize', 14);

        if ~exist('img', 'dir')
            mkdir('img');
        end
        outputFolder = fullfile('img', regionName);
        if ~exist(outputFolder, 'dir')
            mkdir(outputFolder);
        end
        fileName = fullfile(outputFolder, ['grafico_casi_cumulativi_ww' regionName '.tex']);
        matlab2tikz(fileName, ...
           'showInfo', false);

    end

    % se non sono state usate le misure delle WW (i.e., solo CASE)
    % quindi, questo plot serve a valutare quanto bene i dati clinici (casi) riescano a spiegare la dinamica nelle acque reflue.
    if ~useData(2)

        figure('Position', [100, 200, 1200, 450]);
        hold on; grid on;
        Ymean = (1e5 * Yest(2, :) + minYW).^(1 / params.WWexp);
        Ydata = (1e5 * YW(WWinds) + minYW).^(1 / params.WWexp);
        Yhi = 1e5 * Yest(2, :) + minYW + 2 * 1e5 * sqrt(Ysd(2, :));
        Ylo = max(1e5 * Yest(2, :) + minYW - 2 * 1e5 * sqrt(Ysd(2, :)), 0);
        Yhi_transf = Yhi.^(1 / params.WWexp);
        Ylo_transf = Ylo.^(1 / params.WWexp);
        fill([1:length(Yest), fliplr(1:length(Yest))], ...
            [Ylo_transf, fliplr(Yhi_transf)], ...
            [1, 0.93, 0.93], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
        plot(Ymean, '-r', 'LineWidth', 2);                      % Stima
        plot(WWinds, Ydata, 'ok', 'MarkerFaceColor', 'k', 'MarkerSize', 5); % Dati osservati
        xlim([1, length(Yest)]);
        ylim([min(Ylo_transf(:)), max(Yhi_transf(:)) * 1.1]);
        set(gca, 'FontSize', 14, 'Layer', 'top');
        xticks(firsts);
        xticklabels(labs);
        ylabel('RNA in WW', 'FontSize', 16);
        if any(strcmp(params.region, {'Luxembourg'}))
            xlabel('Dates 2020–2021', 'FontSize', 16);
        else
            xlabel('Dates 2020–2023', 'FontSize', 16);
        end
        legend({'2SD envelope','Estimated from case data', 'Data'}, ...
            'Location', 'northeast', 'FontSize', 14);
        
        if ~exist('img', 'dir')
            mkdir('img');
        end
        outputFolder = fullfile('img', regionName);
        if ~exist(outputFolder, 'dir')
            mkdir(outputFolder);
        end
        fileName = fullfile(outputFolder, ['grafico_ww_rna_stimato' regionName '.tex']);
        matlab2tikz(fileName, ...
            'showInfo', false);

        Yaux = (1e5*Yest(2,:)+minYW).^(1/params.WWexp);
        WWaux = (1e5*YW(WWinds)+minYW).^(1/params.WWexp);
        WWcorr = sum((Yaux(WWinds)-mean(Yaux(WWinds))).*(WWaux-mean(WWaux)))/norm((Yaux(WWinds)-mean(Yaux(WWinds))))/norm((WWaux-mean(WWaux)))

    end
            
end

