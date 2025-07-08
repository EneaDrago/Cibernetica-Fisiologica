function [C, labs, firsts, longDates] = SEIRWWinit(YC,day1,specialHolidays,darkNumber)

%Create dates and labels for plots
monthLengths = [31 28 31 30 31 30 31 31 30 31 30 31];
longDates = {};
dates = {};
firsts = [1, 1 + cumsum(monthLengths), 366 + cumsum(monthLengths), 731 + cumsum(monthLengths), 1096 +  cumsum(monthLengths), 1461 + cumsum(monthLengths)];
firsts = firsts(1:end-1);
year = 2019;
for jm = 1:60
    for jd = 1:monthLengths(mod(jm-1,12)+1)
        dd = num2str(jd);
        if length(dd) == 1
            dd = ['0' dd];
        end
        mo = mod(jm-1,12)+1;
        mm = num2str(mo);
        if length(mm) == 1
            mm = ['0' mm];
        end
        longDates = {longDates{1:length(longDates)},[dd '/' mm '/' num2str(year)]};
        dates = {dates{1:length(dates)},[dd '/' mm]};
    end
    if mod(jm,12) == 0
        year = year + 1;
    end
end
longDates = {longDates{1:424}, '29/02/2020', longDates{425:end}};
dates = {dates{1:424}, '29/2', dates{425:end}};
firsts(15:end) = firsts(15:end) + 1;

ind0 = find(strcmp(longDates,day1));
if isempty(ind0)
    error('Start date not recognised! Should be given as DD/MM/YYYY from 01/01/2019 to 31/12/2023.')
end

excl = sum(firsts < ind0);
for jt = 1:length(firsts)  
    if mod(jt,2) == mod(excl+1,2)
        labs{jt} = dates{firsts(jt)};
    else
        if length(YC) < 250
            labs{jt} = dates{firsts(jt)};
        else   
            labs{jt} = ' ';
        end
    end
end

firsts = firsts(excl+1:end) - ind0 + 1;
labs = labs(excl+1:end);
longDates = longDates(ind0:end);


%Output coefficient (weekday-dependent)
C = ones(1,2000);

included = ones(size(C));
included(specialHolidays) = 0;
included(YC<0) = 0;

%Generate C for the first five weeks
for jd = 1:7
    C(jd:7:jd+28) = mean(YC(jd:7:jd+28))/mean(YC(1:35));
end

for jd = 36:length(YC)
    normC = (included(jd-7)+included(jd-14)+included(jd-21)+included(jd-28));
    C(jd) = (included(jd-7)*YC(jd-7) + included(jd-14)*YC(jd-14) + included(jd-21)*YC(jd-21) + included(jd-28)*YC(jd-28));
    
    if normC > 0 && C(jd) > 0
        C(jd) = C(jd)/normC;
        C(jd) = 28*C(jd)/sum(YC(jd-27:jd));
    else
        C(jd) = 1;
    end
    
end
for jd = length(YC)+1:length(C)
    C(jd) = (C(jd-7)+C(jd-14)+C(jd-21)+C(jd-28)+C(jd-35))/5;
end   
    
Cm = movmean(C,[6,0]);
C = C./Cm;   

for jd = 1:size(darkNumber,1)-1
    C(darkNumber(jd,2):darkNumber(jd+1,2)-1) = C(darkNumber(jd,2):darkNumber(jd+1,2)-1)/darkNumber(jd,1);
end
C(darkNumber(end,2):end) = C(darkNumber(end,2):end)/darkNumber(end,1);

C(specialHolidays) = .25*C(specialHolidays); 

   
       
