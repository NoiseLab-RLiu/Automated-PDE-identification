function [X,error] =Thresholding(A,y,K)

% initializing variables
[~,nA] = size(A); % number of columns in A
% Normalization
NORM = zeros(1,nA);
for i =1:nA
    NORM(1,i) = norm(A(:,i));
end
An  = A./repmat(NORM,size(A,1),1); % Normalizing
%r = y;
%R=[];
%k=0;
%while(k<=nA)
    %k=k+1;
aProj = pinv(An)*y;%abs(An'*y); % finding projections
[~,I] = maxk(abs(aProj),K);
xp = pinv(An(:,I))*y; % updating coefficients
x_th = zeros(nA,1);
x_th(I)=xp;
r = y-An*x_th; % calculating residual vector
    %R=[R,norm(r)];
error = norm(r);
%end
%     figure
%     plot(R)

% ERR = y-An*x_th; % calculating residual matrix
X = x_th./NORM';
% %error = mean(mean(abs(R)));
% error = norm(ERR);