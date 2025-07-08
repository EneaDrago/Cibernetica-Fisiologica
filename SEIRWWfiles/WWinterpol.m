function YW = WWinterpol(YW)

inds = find(YW>-.5);
    
for jt = 2:length(inds)
    YW(inds(jt-1)+1:inds(jt)-1) = YW(inds(jt-1)) +   (1:(inds(jt)-inds(jt-1)-1))/(inds(jt)-inds(jt-1))*(YW(inds(jt))-YW(inds(jt-1)));
end



