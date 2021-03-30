function [dropx,dropy,dropt] = drop_entries(N_x,N_y,M,method,param)
% Calculate how many entries should be dropped at each end near the
% boundary. 'param' is difined directly by the number of dropped entries
% for FD and the cosine part portion within a Tukey window for spectral
% method
if(method=='FD')
    dropx = param;
    dropy = param;
    dropt = param;
elseif(method=='SP')
    dropx = ceil(param/2*N_x); 
    dropy = ceil(param/2*N_y);
    dropt = ceil(param/2*M);
else
    sprintf('error')
    [dropx,dropy,dropt] = deal(NaN,NaN,NaN);
end
end

