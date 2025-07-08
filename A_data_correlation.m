% Carica i dati dal file Excel
filename = './data/Toscana_wwtp1.xlsx';
T = readtable(filename);


% Estrai colonne
date = datetime(T.date);
cases = T.cases;
ww = T.ww;

% Considera solo i punti dove ww è valido (≠ -1)
valid_idx = (ww ~= -1) & ~isnan(ww);

% Estrai i dati validi
cases_valid = cases(valid_idx);
ww_valid = ww(valid_idx);
date_valid = date(valid_idx);

% Normalizza entrambi i vettori [0, 1]
cases_norm = (cases_valid - min(cases_valid)) / (max(cases_valid) - min(cases_valid));
ww_norm = (ww_valid - min(ww_valid)) / (max(ww_valid) - min(ww_valid));

% Calcola la correlazione di Pearson
R = corr(cases_norm, ww_norm);

% Visualizza il valore
fprintf('Correlazione (Pearson) usando solo i dati validi: %.4f\n', R);

% Plot
figure; hold on; grid on;
plot(date_valid, cases_norm, '-r', 'LineWidth', 2);
plot(date_valid, ww_norm, '-b', 'LineWidth', 2);
legend({'Normalized Cases', 'Normalized Wastewater'}, 'Location', 'best');
xlabel('Date');
ylabel('Normalized value');
title(sprintf('Correlation (valid data only) = %.4f', R));