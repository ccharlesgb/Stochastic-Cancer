function [data] = anisosmooth(data, kappa, smoothing)
%ANISOSMOOTH Smooths data using anisotropic diffusion
%   kappa is the shape preserving sensitivity. SMALL kappa means HIGH
%   sensitivity to edges. (Preserves peak bandwidth)
%   smoothing is the amount of smoothing to apply. 50-100 provides a high
%   quality spectrum

    timestep = 0.5; %Highest timestep that offers stable solution
    frames = smoothing / timestep;
    size = length(data);
    %Default struct
    DefaultNode = struct();
    DefaultNode.I = 0;
    DefaultNode.I_dot = 0;
    DefaultNode.c = 0;
    DefaultNode.didx = 0;
    DefaultNode.dcdx = 0;
    DefaultNode.didx2 = 0; %Second order
    NA(1:size) = DefaultNode; %Populate array
    
    for m=1:size
        NA(m).I = data(m,2);
    end
    
    for f=1:frames
        for m=2:size - 1    
            NA(m).didx2 = (NA(m+1).I - 2.*NA(m).I + NA(m-1).I);
            NA(m).didx = (NA(m+1).I - NA(m).I); %Forward derivative
            NA(m).c = exp(-(abs(NA(m).didx)/kappa).^2);
            %disp(NA(m).c);
        end
        for m=2:size - 1
            NA(m).dcdx = (NA(m+1).c - NA(m).c); %Forward derivative
            NA(m).I_dot = NA(m).dcdx .* NA(m).didx + NA(m).c .* NA(m).didx2;
        end
        for m=1:size
           NA(m).I = NA(m).I + NA(m).I_dot * timestep;
        end
    end
    for m=1:size
        data(m,2) = NA(m).I;
    end
    
    maxval = max(data(:,2));
    if (maxval < 0.8)
       %fptrintf(1,'Warning smoothing has reduced peak intensities by %.1f%\n', maxval .* 100); 
    end
    %data = normalise(data); %Renormalise
end

% function plotI(NA)
%     for m=1:length(NA)
%         IS(m) = NA(m).I;
%     end
%     plot(IS);
% end

