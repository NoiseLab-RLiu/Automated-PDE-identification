%%
[N_x,N_y,M] = size(U);
Mused = M;
freq_intval = freq_src(2)-freq_src(1);
df = freq_intval;
%% Wave equation extraction
f_ind_min = freq_src(1)/df;
f_ind_max = freq_src(end)/df;
num_used_bin = length(f_ind_min:freq_intval/df:f_ind_max);
num_dict_cols = 11;

ERRT_heq = zeros(num_dict_cols+1,num_used_bin); %Error(T), "+1" is for T=0 sparsity
PHI_heq = zeros(num_dict_cols+1,num_used_bin); %Error(T), "+1" is for T=0 sparsity
X_heq = zeros(num_dict_cols,num_used_bin);
lambda = 1e-0;

%% center source within the region
c_Helmh = zeros((f_ind_max-f_ind_min)/(freq_intval/df)+1,1);

method='SP';
drop_param = 1;
if(strcmp(method,'SP')==1)
    drop_param = 0.4;
end


[dropx,dropy,~] = drop_entries(N_x,N_y,M,method,drop_param);
freq_num = length(freq_src);
LHS_res_freq = cell(freq_num,1);
LHS_err_freq = zeros(num_dict_cols+1,freq_num);
Sparsity_freq = zeros(freq_num,1);
for freq_ind=f_ind_min:freq_intval/df:f_ind_max
    index = (freq_ind - f_ind_min)/(freq_intval/df)+1; % index used in variables res and c
    f_ind = freq_ind+1; % index for positive freq
    Ufx = zeros(N_x,N_y);
    Ufy = zeros(N_x,N_y);
    Ufxy = zeros(N_x,N_y);
    Ufxx = zeros(N_x,N_y);
    Ufyy = zeros(N_x,N_y);
    Uft = zeros(N_x,N_y);
    Uftt = zeros(N_x,N_y);
    Ufid = zeros(N_x,N_y);

    for i=1:N_x %N_x = N_y
        Ufxx(i,:,1) = numder(Uf(i,:,f_ind),dx,2,method);%numder(Utmp(:,i), dx, 1,method);
        Ufx(i,:,1) = numder(Uf(i,:,f_ind),dx,1,method);
        Ufyy(:,i,1) = numder(Uf(:,i,f_ind),dy,2,method);
        Ufy(:,i,1) = numder(Uf(:,i,f_ind),dy,1,method);
    end 
    for i=1:N_y
        Ufxy(:,i,1) = numder(Ufx(:,i,1),dy,1,method);
    end

    
    Ufid(:,:,1) = Uf(:,:,f_ind);    
    U1 = ones(N_x,N_y);
    % define dictionary atoms
    if(length(source_x)~=0)
        Ix = [dropx+1:source_x(1)-1,source_x(2)+1:N_x-dropx];%[21:49,52:80];%dropx+1:N_x-dropx;
        Iy = [dropy+1:source_y(1)-1,source_y(2)+1:N_y-dropy];%[21:49,52:80];%dropy+1:N_y-dropy;
    else
        Ix = dropx+1:N_x-dropx;
        Iy = dropy+1:N_y-dropy;
    end
    veclen = length(Ix)*length(Iy);
    onev = reshape(U1(Ix,Iy,:),veclen,1);
    ufidv = reshape(Ufid(Ix,Iy,:),veclen,1);
    ufxv = reshape(Ufx(Ix,Iy,:),veclen,1);
    ufxxv = reshape(Ufxx(Ix,Iy,:),veclen,1); % The conj revises the error making by FiniteDiff function,
    ufyv = reshape(Ufy(Ix,Iy,:),veclen,1);    
    ufyyv = reshape(Ufyy(Ix,Iy,:),veclen,1);        % please see the Uf and Ufxx5, Ufyy5 in workspace to observe the errors. 
    ufxyv = reshape(Ufxy(Ix,Iy,:),veclen,1);
    uftv = reshape(Uft(Ix,Iy,:),veclen,1); % LHS for frequency f is: -k^2*Uf(f) 
    lhsv = reshape(Uftt(Ix,Iy,:),veclen,1); % LHS for frequency f is: -k^2*Uf(f) 
    % expanded theta
    Theta_e = [onev,ufidv,ufxv,ufyv,ufidv.*ufxv,ufidv.*ufyv,ufxxv,ufxyv,ufyyv,ufidv.*ufxxv,ufidv.*ufxyv,ufidv.*ufyyv];%uftv,
    
    ufxx_ind = 7;
    ufyy_ind = 9;
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
    Sparsity_freq(index) = MININD(b)-1;
   %[X_heq(:,index),err] = OMP_N(Theta_correct,Theta_e(:,b),MININD(b)-1);
    [X_heq(:,index),err] = Thresholding(Theta_correct,Theta_e(:,b),MININD(b)-1);
    if(b<ufxx_ind) % adjust indices for the change from \Phi to \Phi_d
        ufxx_ind=ufxx_ind-1;
        ufyy_ind=ufyy_ind-1;
    elseif(b<ufyy_ind)
        ufyy_ind=ufyy_ind-1;
    end

    omega = 2*pi*freq_ind*df;
    c_Helmh(index) = omega*(sqrt((abs(X_heq(ufxx_ind,index)) + abs(X_heq(ufyy_ind,index)))/2)); % c = (sqrt(c_x^2)+sqrt(c_y^2))/2

end

figure
stem(LHS_err)