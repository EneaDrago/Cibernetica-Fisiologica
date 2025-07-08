function RW0 = RWest(YW,ccc)

inds = find(YW > -.5);
disc = zeros(length(inds)-4,1);

YW(inds) = YW(inds).^ccc*1e-5;

for jj = 1:length(disc)
    disc(jj) = abs(YW(inds(jj+2)) - sum(YW(inds(jj:jj+4)))/5);
end

RW0 = median(disc)^2;



