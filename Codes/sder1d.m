function [k,df]=sder1d(f,dx,beta)
% SDER1D(f,dx) spectral derivative of vector
nx=max(size(f));
f = reshape(f,nx,[]);
% initialize k
kmax=pi/dx;
dk = (2*pi/nx) * (1/dx);

if(mod(nx,2)==0)
    for i=1:fix(nx/2)
        k(i)=(i-1)*dk; 
        k(nx/2+i)=-kmax+(i-1)*dk; 
    end
else
    for i=1:floor(nx/2)
        k(i) = (i-1)*dk;
        k(nx-(i-1)) = -i*dk;
    end
    k(ceil(nx/2)) = (ceil(nx/2)-1)*dk;
end

k=sqrt(-1)*k;
% FFT and IFFT
ff=fft(f); 
ff=(k.^beta).*ff.'; 
df = ifft(ff);
%df=real(ifft(ff));