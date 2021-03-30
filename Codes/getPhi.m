function [PHI,min_ind] = getPhi(ERRT,lambda)
% lambda = 1*1e-3 works well
[r,c] = size(ERRT);
PHI = zeros(r,c);
for i=1:c
    for j=1:r
        PHI(j,i) = ERRT(j,i) + lambda*(j-1)*ERRT(end,i);
    end
end
min_ind = zeros(c,1);
for i=1:c
    [~,min_ind(i)] = min(PHI(:,i));
end
