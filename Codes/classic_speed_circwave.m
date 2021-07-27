%% Parameters already defined in GenWave_analy.m
[~,~,Mused] = size(Uf);
Fmax = Mused-1;%Do not consider DC component
KM = 1000; %zero padding, make dk smaller
NN = Fmax;
kk_save = zeros(KM,KM,Fmax);
for nn = 1:Fmax
    fprintf('Computing Velocity for frequency index %06i / %06i \n', nn, Fmax)
    kk = (fft2(squeeze(Uf(:,:,nn+1)),KM,KM)); % without DC component, Uf is not real
    kk_save(:,:,nn) = kk; 
end



sf_vec = [];
for f=1:Fmax %1:Fmax+1
    s_vec = [];
    for s=0:KM/2-1
        summ = 0;
        for sx=0:s
            tmp = abs(kk_save(sx+1,round(sqrt(s^2-sx^2))+1,f));%+kk_save(sx+1,round(sqrt(s^2-sx^2))+1,Fmax+1-f));
            summ = summ+tmp;
        end
        s_vec = [s_vec,summ];
    end
    sf_vec=[sf_vec,s_vec'];
end

kk_save_sum = zeros(KM,KM);
for i=1:10
    kk_save_sum = kk_save_sum+fftshift(kk_save(:,:,i));
end

sf_vec = [];
for f=1:Fmax %1:Fmax+1
    s_vec = [];
    for s=0:KM-1
        summ = 0;
        for sx=0:s
            tmp = abs(fftshift(kk_save(sx+1,round(sqrt(s^2-sx^2))+1,f)));%+kk_save(sx+1,round(sqrt(s^2-sx^2))+1,Fmax+1-f));
            summ = summ+tmp;
        end
        s_vec = [s_vec,summ];
    end
    sf_vec=[sf_vec,s_vec'];
end

radius=[];
for f=1:length(freq_src)
    [~,kind] = max(sf_vec(1:KM/2,freq_src(f)/df));%sf_vec(1:500,f)
    radius = [radius,kind-1];
end

speeds = [];
for i=1:length(freq_src)
    speeds = [speeds,(freq_src(i))/radius(i)]; % whether df*(i-1) or df*i depends on the start freq used for kk_save
end


figure
plot(freq_src(1:end),speeds(1:end),'r*','MarkerSize',14)
%ylim([0.9*c,1.1*c])
xlabel('Frequency (Hz)')
ylabel('Phase speed (m/s)')
title('c recovered by wavenumber extraction method')
ax = gca;
ax.FontSize = 14; 