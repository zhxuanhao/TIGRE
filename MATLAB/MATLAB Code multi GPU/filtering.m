function [ proj ] = filtering(proj,geo,angles,parker,weights)
%FILTERING Summary of this function goes here
%   Detailed explanation goes here
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% This file is part of the TIGRE Toolbox
% 
% Copyright (c) 2015, University of Bath and 
%                     CERN-European Organization for Nuclear Research
%                     All rights reserved.
%
% License:            Open Source under BSD. 
%                     See the full license at
%                     https://github.com/CERN/TIGRE/license.txt
%                     and
%                     https://www.mathworks.com/matlabcentral/fileexchange/view_license?file_info_id=35548
%
% Contact:            tigre.toolbox@gmail.com
% Codes:              https://github.com/CERN/TIGRE/
% Coded by:           Kyungsang Kim, modified by Ander Biguri 
%--------------------------------------------------------------------------


if parker
	% proj = permute(ParkerWeight(permute(proj,[2 1 3]),geo,angles,parker),[2 1 3]);
    % MODIFICATION, RB, 6/8/2017: We no longer permute before this
    % function, so do not permute before/after parker!
	proj = ParkerWeight(proj,geo,angles,parker);
end 

filt_len = max(64,2^nextpow2(2*geo.nDetector(1)));
[ramp_kernel] = ramp_flat(filt_len);

d = 1; % cut off (0~1)
[filt] = Filter(geo.filter, ramp_kernel, filt_len, d);
filt = repmat(filt',[1 geo.nDetector(2)]);

% MODIFICATION, RB: Let's use GPU to make filtering faster (it's about 4 x faster on
% Quadro P5000 as compared to CPU!)
filt=gpuArray(filt);
% MOD, RB, 6/8/2017: Also apply weights here on GPU; that gives us a
% significant speed-up (obviously the limiting assumption of constant
% detector offsets for all projecitons mnentioned in FDK.m and 
% in WeighAndFilterProjectionsForFDK.m is still valid). These
% simplificaitons give us speed, but reduce flexibility of TIGRE.
weights=gpuArray(weights);

% ADD, RB, 6/8/2017: multiplier to use in the loop below, to simplify the expression!
multiplier=1/2/geo.dDetector(1)*(2*pi/length(angles))/2*(geo.DSD/geo.DSO);

for ii=1:length(angles)
    
    % MODIFICATION, RB: Use gpuArray for speed
    %fproj = (zeros(filt_len,geo.nDetector(2),'single'));
   % fproj = zeros(filt_len,geo.nDetector(2),'single','gpuArray');
    
    % MODIFICATION, RB, 6/8/2017: Swapped dimensions, since we do not
    % permute proj before this funciton (to improve speed)
    fproj = zeros(geo.nDetector(2),filt_len,'single','gpuArray');
    
    % MOD, RB, 6/8/2017: Swapped indexes since we no longer permute proj
    % before this function
    % fproj(filt_len/2-geo.nDetector(1)/2+1:filt_len/2+geo.nDetector(1)/2,:) = proj(:,:,ii);
    % Let's also do WEIGHTING here to improve speed (done on GPU
    % presumably, and it DOES speed things up according to tests)
    fproj(:,filt_len/2-geo.nDetector(1)/2+1:filt_len/2+geo.nDetector(1)/2) = proj(:,:,ii).*weights;
   
    % MODIFICATION, RB, 6/8/2017: Added transpose (on GPU) before processing, since we
    % removed permute before and after filtering to improve speed.
    fproj = fft(fproj');   
    fproj = fproj.*filt; 
    fproj = real(ifft(fproj));
   
    %proj(:,:,ii) = fproj(end/2-geo.nDetector(1)/2+1:end/2+geo.nDetector(1)/2,:)/2/geo.dDetector(1)*(2*pi/  length(angles)   )/2*(geo.DSD/geo.DSO);
    % MODIFICATION, RB: Use gpuArray for speed    
    % MODIFICATION, RB, 6/8/2017: Added transpose (on GPU) just before gather, since we
    % removed permute before and after filtering to improve speed.
    % proj(:,:,ii) = gather((fproj(end/2-geo.nDetector(1)/2+1:end/2+geo.nDetector(1)/2,:)/2/geo.dDetector(1)*(2*pi/length(angles))/2*(geo.DSD/geo.DSO))');  
    % Also let's simplify the final expression by pre-computing the
    % multiplier before this loop:
    proj(:,:,ii) = gather((fproj(end/2-geo.nDetector(1)/2+1:end/2+geo.nDetector(1)/2,:)*multiplier)');  

end  % iterate through all projections

% MOD, RB, 6/8/2017: No longer needed, since we do NOT permute before this
% function call! Done to improve speed.
% proj=permute(proj,[2 1 3]);

end  % Function filtering

function [h, nn] = ramp_flat(n)
nn = [-(n/2):(n/2-1)]';
h = zeros(size(nn),'single');
h(n/2+1) = 1 / 4;
odd = mod(nn,2) == 1;
h(odd) = -1 ./ (pi * nn(odd)).^2;
end


function [filt] = Filter(filter, kernel, order, d)

f_kernel = abs(fft(kernel))*2;
filt = f_kernel(1:order/2+1)';
w = 2*pi*(0:size(filt,2)-1)/order;   % frequency axis up to Nyquist 

switch lower(filter)
    case 'ram-lak'
        % Do nothing
    case 'shepp-logan'
        % be careful not to divide by 0:
        filt(2:end) = filt(2:end) .* (sin(w(2:end)/(2*d))./(w(2:end)/(2*d)));
    case 'cosine'
        filt(2:end) = filt(2:end) .* cos(w(2:end)/(2*d));
    case 'hamming'  
        filt(2:end) = filt(2:end) .* (.54 + .46 * cos(w(2:end)/d));
    case 'hann'
        filt(2:end) = filt(2:end) .*(1+cos(w(2:end)./d)) / 2;
    otherwise
        disp(filter);
        error('Invalid filter selected.');
end

filt(w>pi*d) = 0;                      % Crop the frequency response
filt = [filt , filt(end-1:-1:2)];    % Symmetry of the filter
return

end