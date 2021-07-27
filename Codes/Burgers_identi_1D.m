%%
load('Burgers_u1D.mat')
%dataset_index = 1;
Burgers_u_RK4 = Burgers_u1D{dataset_index}; % \nu=0.025,0.05,0.1
U=real(Burgers_u_RK4);
% dt=t(2)-t(1);
% dx=x(2)-x(1);
[N_x,M] = size(U);
%% Burgers equation extraction
num_used_bin = 1;
num_dict_cols = 7;

ERRT_weq_sp = zeros(num_dict_cols+1,num_used_bin); %Error(T), "+1" is for T=0 sparsity
PHI_weq = zeros(num_dict_cols+1,num_used_bin); %Error(T), "+1" is for T=0 sparsity
X_weq = zeros(num_dict_cols,num_used_bin);
lambda = 1e-0;


freq_num=1;

%%
LHS_res_freq = cell(freq_num,1);
LHS_err_freq = zeros(num_dict_cols+1,freq_num);
CorrectLHS_freq = zeros(freq_num,1);
method = 'FD'; % method used to calculate numeriacl derivative, 'FD' (finite diff) or 'SP' (spectral)

Uused = U;
Mused = M;
source_x=[];
index=1;

% drop indices
if(strcmp(method,'SP'))
    dropx=fix(N_x*0.2);
    dropt=fix(Mused*0.2);
else
    dropx=1;
    dropt=1;
end


% space derivatives
Ux=zeros(N_x,Mused);
Uxx=zeros(N_x,Mused);

for k = 1:Mused
    Utmp=Uused(:,k);
    Ux(:,k)=numder(Utmp, dx, 1,method);
    Uxx(:,k)=numder(Utmp, dx, 2,method);
end

% time derivatives
Ut=zeros(N_x,Mused);
Utt=zeros(N_x,Mused);
for i = 1:N_x
    Ut(i,:) = numder(Uused(i,:),dt,1,method);
    Utt(i,:) = numder(Uused(i,:),dt,2,method);
end

% vectorization
if(length(source_x)~=0)
    Ix = [dropx+1:source_x(1)-1,source_x(2)+1:N_x-dropx];%[21:49,52:80];%dropx+1:N_x-dropx;
else
    Ix = dropx+1:N_x-dropx;
end
It = dropt+1:Mused-dropt;
U1 = ones(length(Ix),length(It));

%veclen = (N_x-2*dropx)*(N_y-2*dropy)*(Mused-2*dropt);   
veclen = length(Ix)*length(It);
utv = reshape(Ut(Ix,It),veclen,1);
uttv = reshape(Utt(Ix,It),veclen,1);
uxv = reshape(Ux(Ix,It),veclen,1);
uxxv = reshape(Uxx(Ix,It),veclen,1);
uidv = reshape(Uused(Ix,It),veclen,1);
vec1 = reshape(U1,veclen,1);

Theta_e = [vec1,utv,uttv,uidv,uxv,uidv.*uxv,uxxv,uidv.*uxxv];%,srcv];
Theta_en = Theta_e;
for i=1:size(Theta_en,2)
    Theta_en(:,i) = Theta_en(:,i)/norm(Theta_en(:,i));
end

PHI_L = zeros(size(Theta_en,2),size(Theta_en,2));
LHS_res = zeros(size(Theta_en,2)-1,size(Theta_en,2));% OMP result when each atom in Theta_e treated as LHS
LHS_err = zeros(size(Theta_en,2),1); % error for OMP-residual of each atom in Theta_e treated as LHS
MININD = zeros(size(Theta_en,2),1);
for L=1:size(Theta_en,2)
    Theta = Theta_en;
    Theta(:,L) = []; % delete L-th atom
    lhs = Theta_en(:,L);
    CVERR = CrossValid5(Theta,lhs);
    CV_err = CVERR(:,end);
    [PHI_L(:,L),minind] = getPhi(CV_err,lambda); 

    [LHS_res(:,L),LHS_err(L)] = Thresholding(Theta,lhs,minind-1);
    LHS_err(L) = LHS_err(L)/norm([1;LHS_res(:,L)],2);

    MININD(L) = minind;
end
[a,b]=min(LHS_err);
Theta_correct = Theta_e;
Theta_correct(:,b) = [];
LHS_res_freq{index} = LHS_res;
LHS_err_freq(:,index) = LHS_err;
CorrectLHS_freq(index) = MININD(b)-1;
[X_weq(:,index),err] = Thresholding(Theta_correct,Theta_e(:,b),MININD(b)-1);
