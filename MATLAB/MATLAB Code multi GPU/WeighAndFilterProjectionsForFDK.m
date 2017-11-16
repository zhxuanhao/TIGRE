function [proj_filt]=WeighAndFilterProjectionsForFDK(proj,geo,angles,varargin)

% RB, May 2017. Auxiliary function used by FDKLargeData(). It applies weights and ramp filters
% the projection stack for FDK. This is just code extracted from TIGER'S
% FDK() function.

%TODO docs FDK
% 
%
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
[filter,parker]=parse_inputs(proj,geo,angles,varargin);
geo.filter=filter;

%Input is data,geosize,angles



if size(geo.offDetector,2)==1
    offset=repmat(geo.offDetector,[1 length(angles)]);
else
    offset=geo.offDetector;
end

%% Weight
%proj=data
% MODIFICATION, RB, 01/13/2017: compute weights only once before the loop
% instead oinside the loop. This speeds weighting considerably (nearly 8 times)! The
% assumption is that detector offsets are CONSTANT for all projections,
% which may not be generally the case - SO THIS MODIFICATION CAN INTRODUCE
% PROBLEMS AND IS LIMITING. Should be OK for MRD's purposes though. This
% also affects us and vs creation (we assume constant for all
% projections).

us = ((-geo.nDetector(1)/2+0.5):1:(geo.nDetector(1)/2-0.5))*geo.dDetector(1) + offset(1,1);
vs = ((-geo.nDetector(2)/2+0.5):1:(geo.nDetector(2)/2-0.5))*geo.dDetector(2) + offset(2,1);
[uu,vv] = meshgrid(us,vs); %detector

%Create weight according to each detector element
% w = (geo.DSD)./sqrt((geo.DSD)^2+uu.^2 + vv.^2)';
% MODIFICATION, RB, 6/8/2017: We no longer permute below, so no transpose
% in weights!
w = (geo.DSD)./sqrt((geo.DSD)^2+uu.^2 + vv.^2);

% DO NOT PERMUTE, RB, 6/8/2017 - to improve speed.
% proj=permute(proj,[2 1 3]);

% MODIFICATION, 6/8/2017: We do weighting in the filtering() function on
% GPU, which is significantly faster! The limitation and assumption of
% offsets constancy for all projections remains true!

% for ii=1:length(angles)
% 
% % MODIFICATION, RB: Do NOT compute weights here: assume constancy
% % of detector offsens and compute them only once before this loop. This
% % also affects us and vs creation (we assume constant for all
% % projections).
%     
% %     us = ((-geo.nDetector(1)/2+0.5):1:(geo.nDetector(1)/2-0.5))*geo.dDetector(1) + offset(1,ii);
% %     vs = ((-geo.nDetector(2)/2+0.5):1:(geo.nDetector(2)/2-0.5))*geo.dDetector(2) + offset(2,ii);
% %     [uu,vv] = meshgrid(us,vs); %detector
% %     
% %     %Create weight according to each detector element
% %     w = (geo.DSD)./sqrt((geo.DSD)^2+uu.^2 + vv.^2);
%     
%     %Multiply the weights with projection data
%     % proj(:,:,ii) = proj(:,:,ii).*w';
%     proj(:,:,ii) = proj(:,:,ii).*w;  % RB: transpose is done before the loop!
% end

%% filter
% MOD, RB, 6/8/2017: Pass weights to the function to do weighting on GPU
% for speed!
proj_filt = filtering(proj,geo,angles,parker,w); % Not sure if offsets are good in here
%toc;
% proj_filt = permute(proj,[2 1 3]); % Not sure if offsets are good in here
%RMFIELD Remove fields from a structure array.
geo=rmfield(geo,'filter');


end

function [filter, parker]=parse_inputs(proj,geo,alpha,argin)
opts=     {'filter','parker'};
defaults=ones(length(opts),1);
% Check inputs
nVarargs = length(argin);
if mod(nVarargs,2)
    error('CBCT:SART:InvalidInput','Invalid number of inputs')
end

% check if option has been passed as input
for ii=1:2:nVarargs
    ind=find(ismember(opts,lower(argin{ii})));
    if ~isempty(ind)
        defaults(ind)=0;
    end
end

for ii=1:length(opts)
    opt=opts{ii};
    default=defaults(ii);
    % if one option isnot default, then extranc value from input
    if default==0
        ind=double.empty(0,1);jj=1;
        while isempty(ind)
            ind=find(isequal(opt,lower(argin{jj})));
            jj=jj+1;
        end
         if isempty(ind)
            error('CBCT:FDK:InvalidInput',['Optional parameter "' argin{jj} '" does not exist' ]); 
        end
        val=argin{jj};
    end
    
    switch opt
        % % % % % % % Verbose
        case 'parker'
            if default
                parker=0;
            else
                parker=val;
            end
           
        case 'filter'
            if default
                filter='ram-lak';
            else
                if  ~ischar( val)
                    error('CBCT:FDK:InvalidInput','Invalid filter')
                end
                filter=val;
            end
       
        otherwise
            error('CBCT:FDK:InvalidInput',['Invalid input name:', num2str(opt),'\n No such option in FAK()']);
    end
end

end