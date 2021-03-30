function deriv = numder(sig,step,order,option)
% function to generate numerical derivatives by finite difference(FD) or
% spectral method(SP). The spectral method requires a Tukey window for
% preprocessing.
%   
if(option=='FD')
    deriv = FiniteDiff(sig,step,order);
elseif(option=='SP')
    w = tukeywin(length(sig),0.4);
    [~,deriv] = sder1d(reshape(sig,[length(sig),1]).*w,step,order);
else
    sprintf('Error')
    deriv = NaN;
end
    
end

