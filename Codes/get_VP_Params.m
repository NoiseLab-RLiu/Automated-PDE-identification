N_x = 100;
N_y = 100;
M = 3001;
Mused = 300;
dx = .001;
dy = .001;
dt = 1/300000;
df = (1/dt)/Mused;

U = zeros(N_x,N_y,M);
for i=1:N_x
    for j=1:N_y
        U(i,j,:) = ifft(Xsus(:,i,j));
    end
end

% figure
% for i=1:300
%     imagesc(U(:,:,i))
%     axis square
%     pause(.1)
% end

maxval = max(max(max(U)));
minval = min(min(min(U)));
%% show frames
% figure
% for t=1:Mused
%     imagesc(U(:,:,t))
%     axis square
%     xlabel('x')
%     ylabel('y')
%     caxis([minval,maxval])
%     colorbar
%     title(strcat('frame ',num2str(t-1)))
%     pause(.1)
% end
%%
Uf = zeros(N_x,N_y,Mused);
for i=1:N_x
    for j=1:N_y
        Uf(i,j,:) = fft(U(i,j,1:Mused));
    end
end

freq_src = 0:df:Mused/2*df;