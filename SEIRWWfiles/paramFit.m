
function [Js, nu, RW0] = paramFit(params,YC,YW,C,gs,es,ns)

%Set parameters
params.gamma = gs;
params.WWexp = es;


%Estimate wastewater measurement variance
RW0 = RWest(YW,params.WWexp);
params.RW0 = RW0;
params.RW = params.RW0/10;

%If nu is positive use that
if ns > 0
    params.nu = ns;

%If nu is negative, it indicates that it is estimated by LS from the direct
%problem (WW reconstruction based on case data)
else
    
    %Find nu
    params.nu = 1;
    Yest = SEIR_WW(params,YC,YW,C,[true,false],1000);
    WWinds = find(YW > -.5);
    Y2 = 1e-5*YW(WWinds).^es;
    YCaux = sort(YC,'ascend');
    YWaux = sort(Y2,'ascend');
    ccc = mean(YCaux(1:floor(length(YC)/10)))/mean(YCaux(floor(length(YC)/10)+1:end));
    aaa = mean(YWaux(1:floor(length(YWaux)/10)));
    bbb = mean(YWaux(floor(length(YWaux)/10)+1:end));

    if ccc*bbb < aaa
        Y2 = Y2 - min((aaa-ccc*bbb)/(1-ccc),min(Y2));
    end

    XX = Yest(2,WWinds);
    nu = Y2*XX'*(XX*XX')^-1;
    params.nu = nu;

end

Yest = SEIR_WW(params,YC,YW,C,[false, true],1000);
Js = norm(Yest(1,:)-YC)^2; 

nu = params.nu;
