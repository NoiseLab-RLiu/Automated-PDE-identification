%%
freq_intval = (freq_src(end)-freq_src(1))/(length(freq_src)-1);
if(isnan(freq_intval))
    freq_intval=eps;
end

%% Wave equation extraction
df = freq_intval;
f_ind_min = freq_src(1)/df;
f_ind_max = freq_src(end)/df;
num_used_bin = length(f_ind_min:freq_intval/df:f_ind_max);
num_dict_cols = 12;

ERRT_weq_sp = zeros(num_dict_cols+1,num_used_bin); %Error(T), "+1" is for T=0 sparsity
PHI_weq = zeros(num_dict_cols+1,num_used_bin); %Error(T), "+1" is for T=0 sparsity
X_weq = zeros(num_dict_cols,num_used_bin);
lambda = 1e-1;

[N_x,N_y,~] = size(U);
c = zeros((f_ind_max-f_ind_min)/(freq_intval/df)+1,1);

%%
method = 'FD'; % method used to calculate numeriacl derivative, 'FD' (finite diff) or 'SP' (spectral)
for freq_ind = f_ind_min:freq_intval/df:f_ind_max
    index = (freq_ind - f_ind_min)/(freq_intval/df)+1;
    Mused=300; % This is the video length used for ideal bandpassing
    Ufbp = zeros(N_x,N_y,Mused);  % band-passed Uf
    Ufbp(:,:,freq_ind+1) = Uf(:,:,freq_ind+1); 
    Ufbp(:,:,Mused-freq_ind+1) = Uf(:,:,Mused-freq_ind+1); 

    Urec = zeros(N_x,N_y,Mused); % recovered band-limited time-domain U
    for i=1:N_x
        for j=1:N_y
            Urec(i,j,:) = ifft(Ufbp(i,j,:));
        end
    end
    start_frame = 30; % Change the starting frame of selected period for PDE recovery
    Mused = 100; % This is the length used for PDE identification
    Uused = (Urec(:,:,start_frame+1:start_frame+Mused));

    [dropx,dropy,dropt] = drop_entries(N_x,N_y,Mused,method,1);%last parameter: if "FD", change to 1; otherwise 0.4
    
    % space derivatives
    Ux=zeros(N_x,N_y,Mused);
    Uxx=zeros(N_x,N_y,Mused);
    Uy = zeros(N_x,N_y,Mused);
    Uyy=zeros(N_x,N_y,Mused);
    Uxy=zeros(N_x,N_y,Mused);
    for k = 1:Mused
        Utmp=Uused(:,:,k);
        for i = 1:N_y
            Ux(:,i,k)=numder(Utmp(:,i), dx, 1,method);
            Uxx(:,i,k)=numder(Utmp(:,i), dx, 2,method);
        end
        for i = 1:N_x
            Uy(i,:,k)=numder(Utmp(i,:), dy, 1,method);
            Uyy(i,:,k)=numder(Utmp(i,:), dy, 2,method);
        end
        Utmp=Ux(:,:,k);
        for i = 1:size(Ux,1)
            uxy=numder(Utmp(i,:), dy, 1,method);
            Uxy(i,:,k)=uxy;
        end    
    end

    % time derivatives
    Ut=zeros(N_x,N_y,Mused);
    Utt=zeros(N_x,N_y,Mused);
    for i = 1:N_x
        for j=1:N_y
            Ut(i,j,:) = numder(Uused(i,j,:),dt,1,method);
            Utt(i,j,:) = numder(Uused(i,j,:),dt,2,method);
        end
    end
    
    % vectorization
    if(length(source_x)~=0)
        Ix = [dropx+1:source_x(1)-1,source_x(2)+1:N_x-dropx];%[21:49,52:80];%dropx+1:N_x-dropx;
        Iy = [dropy+1+source_y(1)-1,source_y(2)+1:N_y-dropy];%[21:49,52:80];%dropy+1:N_y-dropy;
    else
        Ix = dropx+1:N_x-dropx;
        Iy = dropy+1:N_y-dropy;
    end
    It = dropt+1:Mused-dropt;
    U1 = ones(length(Ix),length(Iy),length(It));
    
    %veclen = (N_x-2*dropx)*(N_y-2*dropy)*(Mused-2*dropt);   
    veclen = length(Ix)*length(Iy)*length(It);
    utv = reshape(Ut(Ix,Iy,It),veclen,1);
    uttv = reshape(Utt(Ix,Iy,It),veclen,1);
    uxv = reshape(Ux(Ix,Iy,It),veclen,1);
    uxxv = reshape(Uxx(Ix,Iy,It),veclen,1);
    uyv = reshape(Uy(Ix,Iy,It),veclen,1);
    uyyv = reshape(Uyy(Ix,Iy,It),veclen,1);
    uxyv = reshape(Uxy(Ix,Iy,It),veclen,1);
    uidv = reshape(Uused(Ix,Iy,It),veclen,1);
    vec1 = reshape(U1,veclen,1);
    
    Theta_e = [vec1,utv,uttv,uxv,uyv,uidv.*uxv,uidv.*uyv,uxxv,uxyv,uyyv,uidv.*uxxv,uidv.*uxyv,uidv.*uyyv];%,srcv];
    ufxx_ind = 8;
    ufyy_ind = 10;
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
        if(minind~=1)
            [LHS_res(:,L),LHS_err(L)] = OMP_N(Theta,lhs,minind-1);
        else
            LHS_res(:,L) = zeros(size(Theta_en,2)-1,1);
            LHS_err(L) = Inf;
        end
        MININD(L) = minind;
    end
    [a,b]=min(LHS_err);% select LHS
    Theta_correct = Theta_e;
    Theta_correct(:,b) = [];
    [X_weq(:,index),err] = OMP_N(Theta_correct,Theta_e(:,b),MININD(b)-1);
    if(b<ufxx_ind) % adjust indices for the change from \Phi to \Phi_d
        ufxx_ind=ufxx_ind-1;
        ufyy_ind=ufyy_ind-1;
    elseif(b<ufyy_ind)
        ufyy_ind=ufyy_ind-1;
    end
    c(index) = sqrt((abs(X_weq(ufxx_ind,index))+abs(X_weq(ufyy_ind,index)))/2);
end
