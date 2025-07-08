% Predictions using only case data
predsCase = zeros(length(WWinds),1);
params.RW = params.RW0/10;
[~, Xend, P] = SEIR_WW(params,YC,YW,C,[true, false],WWinds,[],[],false);
for jd = 1:length(WWinds)
    Yest = SEIR_WW_FWD(Xend(:,jd),C,P,WWinds(jd)+1,params,winLength);
    predsCase(jd) = sum(Yest(1,:));
end


% Predictions using only WW data
params.RW = params.RW0/10;
predsWW = zeros(length(WWinds),1);
[~, Xend, P] = SEIR_WW(params,YC,YW,C,[false, true],WWinds,[],[],false);
for jd = 1:length(WWinds)
    Yest = SEIR_WW_FWD(Xend(:,jd),C,P,WWinds(jd)+1,params,winLength);
    predsWW(jd) = sum(Yest(1,:));
end

% Predictions using only WW data interpolated
YWip = WWinterpol(YW);
params.RW = params.RW0/10; 
predsWWip = zeros(length(WWinds),1);
[~, Xend, P] = SEIR_WW(params,YC,YWip,C,[false, true],WWinds,[],[],false);
for jd = 1:length(WWinds)
    Yest = SEIR_WW_FWD(Xend(:,jd),C,P,WWinds(jd)+1,params,winLength);
    predsWWip(jd) = sum(Yest(1,:));
end

% Predictions using case & WW data
params.RW = params.RW0;
predsBoth = zeros(length(WWinds),1);
[~, Xend, P] = SEIR_WW(params,YC,YW,C,[true, true],WWinds,[],[],false);
for jd = 1:length(WWinds)
    Yest = SEIR_WW_FWD(Xend(:,jd),C,P,WWinds(jd)+1,params,winLength);
    predsBoth(jd) = sum(Yest(1,:));
end

% Calculate normalised errors
ersCase = (trueCases-predsCase)./trueCases.^.5;
ersWW = (trueCases-predsWW)./trueCases.^.5;
ersWWip = (trueCases-predsWWip)./trueCases.^.5;
ersBoth = (trueCases-predsBoth)./trueCases.^.5;



