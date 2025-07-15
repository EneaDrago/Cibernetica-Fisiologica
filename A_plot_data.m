% Carica i dati dal file Excel
filename = './data/wwtp1.xlsx';
T = readtable(filename);

% Estrai le colonne
date = datetime(T.date);
cases = T.cases;
ww = T.ww;

% Sostituisci -1 con NaN nei dati wastewater
% ww(ww == -1) = NaN;

% Normalizzazione su scala [0, 1] ignorando i NaN
cases_norm = (cases - min(cases)) / (max(cases) - min(cases));
ww_norm = (ww - min(ww,[],'omitnan')) / (max(ww,[],'omitnan') - min(ww,[],'omitnan'));

% Plot
figure; hold on; grid on;
plot(date, cases_norm, '-r', 'LineWidth', 2);
plot(date, ww_norm, '-b', 'LineWidth', 2);
legend({'Normalized Cases', 'Normalized Wastewater'}, 'Location', 'best');
xlabel('Date');
ylabel('Normalized value');
title('Normalized Comparison of Cases and Wastewater Signal');
