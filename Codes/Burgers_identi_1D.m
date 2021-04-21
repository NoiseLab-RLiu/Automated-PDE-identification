%%
%U=Burgers_u(:,:,2:2:150);
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
% maxval = max(max(max(Urec)));
% minval = min(min(min(Urec)));
% figure
% for t = 1:M
%     imagesc(squeeze(Urec(:,:,t)))
%     axis square
%     xlabel('x')
%     ylabel('y')
%     caxis([minval,maxval])
%     colorbar
%     pause(.1)
% end
% figure(1)
% plot(squeeze(Urec(50,50,:)))
% figure(3)
% plot(squeeze(Urec(50,:,10)))
% figure
% imagesc(Urec(:,:,10))
% axis square
%%
LHS_res_freq = cell(freq_num,1);
LHS_err_freq = zeros(num_dict_cols+1,freq_num);
CorrectLHS_freq = zeros(freq_num,1);
method = 'FD'; % method used to calculate numeriacl derivative, 'FD' (finite diff) or 'SP' (spectral)
tic

Uused = U;
Mused = M;
source_x=[];
index=1;
%Uused = Uused/10000;

% For PS 
dropx=fix(N_x*0.2);
dropt=fix(Mused*0.2);
% For FD
dropx=1;
dropt=1;


% space derivatives
Ux=zeros(N_x,Mused);
Uxx=zeros(N_x,Mused);

for k = 1:Mused
    Utmp=Uused(:,k);
    for i = 1:Mused
        Ux(:,k)=numder(Utmp, dx, 1,method);
        Uxx(:,k)=numder(Utmp, dx, 2,method);
    end
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

%Theta_e = [uftv,lhsv,onev,ufxv,ufyv,ufxxv,ufxyv,ufyyv,ufidv.*ufxv,ufidv.*ufyv,ufidv.*ufxxv,ufidv.*ufxyv,ufidv.*ufyyv];%uftv,
Theta_en = Theta_e;
for i=1:size(Theta_en,2)
    Theta_en(:,i) = Theta_en(:,i)/norm(Theta_en(:,i));
end
%     covmat = abs(Theta_en'*Theta_en)
%     figure 
%     imagesc(covmat)
%     axis square
%     colorbar
%     xticks([1:14])
%     yticks([1:14])
%     xticklabels({'1','u_t','u_{tt}','u_x','u_y','uu_x','uu_y','u_{xx}','u_{xy}','u_{yy}','uu_{xx}','uu_{xy}','uu_{yy}','u'})
%     yticklabels({'1','u_t','u_{tt}','u_x','u_y','uu_x','uu_y','u_{xx}','u_{xy}','u_{yy}','uu_{xx}','uu_{xy}','uu_{yy}','u'})
%     title('Correlation Matrix (absolute value)')
%     ax=gca;
%     ax.FontSize=15

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
    %if(minind~=1)
        [LHS_res(:,L),LHS_err(L)] = OMP_N(Theta,lhs,minind-1);
%         else
%             LHS_res(:,L) = zeros(size(Theta_en,2)-1,1);
%             LHS_err(L) = Inf;
%         end
    MININD(L) = minind;
end
[a,b]=min(LHS_err);
%b = 2;
%MININD(b)=3;
Theta_correct = Theta_e;
Theta_correct(:,b) = [];
LHS_res_freq{index} = LHS_res;
LHS_err_freq(:,index) = LHS_err;
CorrectLHS_freq(index) = MININD(b)-1;
[X_weq(:,index),err] = OMP_N(Theta_correct,Theta_e(:,b),MININD(b)-1);

%CVERR = CrossValid5(Theta,uttv);
% 
%ERRT_weq_sp(:,index) = CVERR(:,end);
%[PHI_weq(:,index),minind] = getPhi(ERRT_weq_sp(:,index),lambda); 
%     minind=4;
%     X_weq(:,index) = OMP_N(Theta,uttv,minind-1); %Note 1st entry of PHI_weq is for sparsity 0.
%     Uxx_comb = [Uxx_comb;uxxv];
%     Uyy_comb = [Uyy_comb;uyyv];
%     Utt_comb = [Utt_comb;uttv];
toc

figure
plot(log10(LHS_err))
figure
plot(log10(PHI_L(:,2)))

corrmat=Theta_en'*Theta_en;
figure
imagesc(abs(corrmat))
axis square
colorbar

% [omp_coef,err,R]=OMP_N(Theta_correct,Theta_e(:,b),3);
% 
% corrcoef(R,uidv.*uyyv)
% corrcoef(R,uxxv)
% corrr=[];
% Theta_c_en = Theta_en;
% Theta_c_en(:,b)=[];
% for i=1:12
%     corrmat=corrcoef(R,Theta_c_en(:,i));
%     corrr=[corrr,corrmat(2,1)];
% end
% figure
% stem(corrr)
% 
% pinv(Theta_c_en)*Theta_en(:,b)
%%
% figure
% plot((PHI_L(:,3)))
% figure
% subplot(2,2,1)
% imagesc(real(Ut(Ix,Iy,10)))
% axis square
% colorbar
% xlabel('y')
% ylabel('x')
% title('Ut')
% subplot(2,2,2)
% imagesc(real(Uxx(Ix,Iy,10)+Uyy(Ix,Iy,10)))%Uxx(Ix,Iy,1)+
% axis square
% colorbar
% xlabel('y')
% ylabel('x')
% title('Uxx+Uyy')
% subplot(2,2,3)
% imagesc(real(Utt(Ix,Iy,10)))
% axis square
% colorbar
% xlabel('y')
% ylabel('x')
% title('Utt')
% subplot(2,2,4)
% imagesc(real(Uy(Ix,Iy,10)))
% axis square
% colorbar
% xlabel('y')
% ylabel('x')
% title('Uy')
% 
% figure
% imagesc(real(0.25*Uxx(Ix,Iy,10)+0.25*Uyy(Ix,Iy,10)))
% axis square
% colorbar
% xlabel('y')
% ylabel('x')
% title('c^2(Uxx+Uyy)')
% % figure
% % imagesc(real(Utt(Ix,Iy,10)))
% % axis square
% % colorbar
% % xlabel('y')
% % ylabel('x')
% % title('c^2(Uxx+Uyy)')
% figure
% imagesc(real(U(Ix,Iy,10)))
% axis square
% colorbar
% xlabel('y')
% ylabel('x')
% title('U')
% 
% figure
% imagesc(real(Uxy(Ix,Iy,10)))
% axis square
% colorbar
% xlabel('y')
% ylabel('x')
% title('Uxy')
% % subplot(2,2,2)
% % imagesc(real(Ut(Ix,Iy,10)))
% % axis square
% % colorbar
% 
% figure
% imagesc(Utt(Ix,Iy,10)-0.25*(Uxx(Ix,Iy,10)+Uyy(Ix,Iy,10)))
% axis square
% colorbar
% xlabel('y')
% ylabel('x')
% title('Utt-c^2(Uxx+Uyy)')
%%
% % Uxx = zeros(100,100);
% % Ut = zeros(100,100);
% % Utt = zeros(100,100);
% % for i=1:100
% %     Uxx(:,i) = numder(Uused(:,i), dx, 2,'SP');
% %     Ut(i,:) = numder(Uused(i,:), dt, 1,'SP');
% %     Utt(i,:) = numder(Uused(i,:), dt, 2,'SP');
% % end
% % Ix = 20:80;%2:99;
% % It = 20:80;%2:99;
% % veclen = length(Ix)*length(It);
% % uxxv = reshape(Uxx(Ix,It),veclen,1);
% % uttv = reshape(Utt(Ix,It),veclen,1);
% % utv = reshape(Ut(Ix,It),veclen,1);
% % coef=OMP_N([utv,uxxv],uttv,2)
% figure
% imagesc(Urec(:,:,50))
% axis square
% 
% Utf = zeros(N_x,N_y,M);
% for i=1:N_x
%     for j =1:N_y
%         Utf(i,j,:)=fft(Ut(i,j,:));
%     end
% end
% 
% f_ind = 6
% figure
% subplot(2,2,1)
% imagesc(abs(Utf(:,:,f_ind)))
% axis square
% title('f(Ut)')
% colorbar
% subplot(2,2,2)
% imagesc(abs(1i*(2*pi*freq_ind/(Mused*dt))*Uf(:,:,f_ind)))
% axis square
% title('Ut in Uf')
% colorbar
% subplot(2,2,3)
% imagesc(abs(Uf(:,:,f_ind)))
% axis square
% title('Uf')
% colorbar
% subplot(2,2,4)
% imagesc(abs(-(2*pi*freq_ind/(Mused*dt))^2*Uf(:,:,f_ind)))
% axis square
% title('Utt in Uf')
% colorbar
% 
% if(dispersive==0) % only work for non-dispersive waves
%     coef_comb = OMP_N([Uxx_comb,Uyy_comb],Utt_comb,2);
%     c_comb_w = real(sqrt(coef_comb(1))+sqrt(coef_comb(2)))/2;
% end
% 
% 
% uxx_ind = 5;
% uyy_ind = 7;
% speed = zeros(num_used_bin,1);
% for i=1:num_used_bin
%     vel = (sqrt(real(X_weq(uxx_ind,i)))+sqrt(real(X_weq(uyy_ind,i))))/2;
%     speed(i) = vel;
% end
% 
% figure
% plot(freq_src(1:end),speed(1:end),'r*','MarkerSize',14)
% if(dispersive==0)
%     ylim([0.9*c(1),1.1*c(1)])
% end
% xlabel('Frequency (Hz)')
% ylabel('Phase speed (m/s)')
% title('c recovered by PDE recovery in time domain')
% ax = gca;
% ax.FontSize = 14; 
% 
% % figure
% % plot(Urec(50,:,50))
%% test feb 25
% unorm1=zeros(30,1);
% unorm2=zeros(30,1);
% for i=1:30
%     ures1= Urec(:,:,1:i);
%     ures2 = Urec(:,:,101:100+i);
%     unorm1(i) = norm(ures1(:));
%     unorm2(i) = norm(ures2(:));
% end
% figure
% plot(unorm1)
% hold on
% plot(unorm2)