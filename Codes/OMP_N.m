function [X,error] =OMP_N(A,Y,K)
% [x,error] = OMP_N(A,y,K)
% Finds columns in A that best describe data using OMP
% ****ASSUMES COLUMN NORMALIZED A******
% ---- Inputs
% A = dictionary
% y = observations
% K = sparsity
% ---- Outputs
% x = coefficient vector
% error = l2-squared error for reconstruction
%
% Mike Bianco 3/6/16
% Ruixian Liu 3/29/21


% initializing variables
[~,nA] = size(A); % number of columns in A
[~,nY] = size(Y); % number of data vectors

% Normalization
NORM = zeros(1,nA);
for i =1:nA
    NORM(1,i) = norm(A(:,i));
end
An  = A./repmat(NORM,size(A,1),1);

Xn = zeros(nA,nY);
if(K~=0)
for n = 1:nY
    y = Y(:,n);
    x = zeros(nA,1);
    ai = [];
    for m = 1:K
        r = y-An*x; % calculating residual vector
        aProj = abs(An'*r); % finding projections
        aProj(ai) = -1; %making sure no index is used twice (impose 'orthogonality')
        [~,I] = max(aProj); % max projection corresponds to optimal column
        ai = [ai,I]; % accumulating indices
        xp = pinv(An(:,ai))*y; % updating coefficients
        x(ai)=xp;
    end
    Xn(:,n)=x;
end
end
R = Y-An*Xn; % calculating residual matrix
X = Xn./NORM';
%error = mean(mean(abs(R)));
error = norm(R);