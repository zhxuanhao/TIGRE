function [noVolSeg, noProjSeg] = FindNoOfVolAndProjSegmentsForBackprojection(volDim, projDim, splitProjBetweenGPUs, noOfGPUs, gpuRAMinBytes)
% RB, May 2017.

% Function computes the number of necessary volume and projection set
% segments needed to backproject the data on the GPU (in other words:
% segments that fit into the GPU memory).

% volDim and projDim - dimensions of the full volume to be reconstructed and the
% full projection set in voxels/pixels. Both are vectors of three elements
% (3D volume and 3D projection stack, where the 3rd dimension is the number of projections)

% splitProjBetweenGPUs - if 1, we compute numbers of segments assuming that
% at the lowest (MEX) level we split projection data between GPUs and process in parallel. If 0, we
% assume that at the lowest level we split volume data between GPUs and
% process in parallel.

% The last two parameters are optional. If not specified, we assume we have
% Parallel Processing Toolbox and can get the number of GPUs and GPU RAM
% sizes using the toolbox functionality. Otherwise, the caller can specify
% both the number of GPUs and GPU RAM size.

% gpuRAMinBytes is assumed to be the NOMINAL (total) amount of GPU RAM known to the
% caller (e.g. 1.6106e+10 (1024^3*15) for 15 GB total GPU memory)

% The simplifying assumption made here is that IF MULTIPLE GPUs ARE PRESENT, 
% THEY HAVE THE SAME OR VERY SIMILAR AMOUNT OF RAM.

noVolSeg=1;
noProjSeg=1;

% The last two parameters are optional. If not specified, we assume we have
% Parallel Processing Toolbox

if nargin==4
    error('CBCT:FindNoOfVolAndProjSegmentsForBackprojection:InvalidInput: GPU RAM in bytes needs to be specified if the number of GPUs is given!');
end

if nargin < 4
    noOfGPUs=0;
    gpuRAMinBytes=0;
end

% Compute sizes of the full volume to be reconstructed and the
% full projection set IN BYTES. We compute them from the volume and
% projeciton dimensions (including number of projections) and knowing that
% we have 4 bytes per voxel/projection pixel (single).
volSize=volDim(1)*volDim(2)*volDim(3)*4;  % Singles are 4 bytes per element
projSize=projDim(1)*projDim(2)*projDim(3)*4;

% If noOfGPUs was specified as 0, we assume we have Parallel Processing Toolbox
% and want to figure out how many we have using the toolbox
% Otherwise we just use the passed nominal values and ASSUME SOME SAFETY
% MARGIN

if noOfGPUs>0  % Use user-specified number of GPUs and GPU RAM size
    N = noOfGPUs;
    disp([num2str(N) ' GPUs specified for use']);
    % We assume the user is passing nominal GPU RAM amnount. We know that
    % in reality less RAM is available, so assume some safety margin
    minAvailMem=0.75*gpuRAMinBytes;
    
else  % Figure things out using Parallel Processing toolbox
 
    % As of May 2017 (prototype v.1.0) we only look for Quadro P500 to
    % simplify things
    N = GetNumberOfSpecificGPUs('P5000');
    disp([num2str(N) ' GPU(s) available']);
    % Look through all GPUs and check the smallest amount of available memory
    % (to be conservative). I am assuming that all GPUs are the same (same
    % total memory), but may have slightly different amounts of available
    % memory
    % As of May 2017 (prototype v.1.0) we only look for memory on Quadro P500 to
    % simplify things
    minAvailMem=GetMinAvailGPUMemoryForSpecificGPUType('P5000');

    % Here also let's give ourselves some memory headroom (assume a little
    % bit less memory), although not as severe as in case of the nominal
    % RAM size passed by the caller
    minAvailMem=0.9*minAvailMem;

end  % If we have parallel processing toolbox

% We use an approach from pp. 44-45 of RB's Notebook #183: we compute the
% number of necessary segments that make the processed data elements smaller than the
% available GPU memory and that keeeps the ratio of volume and projection data in GPU memory
% close to 1 (balanced split between volume and projections). We also try to not split the smaller data element (volume or
% projections) if only possible (splitting only one element at his high level has less
% overhead).


% If splitProjBetweenGPUs==1, we use formula:
% volSize/noVolSeg + projSize/(N*noProjSeg) < minAvailMem (valid when we
% split PROJECTIONS between multiple GPUs)
% If splitProjBetweenGPUs==0, we use formula:
% volSize/(N*noVolSeg) + projSize/noProjSeg < minAvailMem 
% We assume that A = volSize/noVolSeg and B = projSize/(N*noProjSeg) 
% (N can be under volSize if splitProjBetweenGPUs==0) and
% that we try to keep the ratio of A/B close to 1, so A = B basically.
% Having these two equations, we can compute noVolSeg and noProjSeg (two
% unknowns) and round them properly to be on the safe side.
% Also, we first compute no of segments of the SMALLER data element, round
% it DOWN, then compute the number of segments of the larger data element
% using the rounded value of the first number of segments and round it UP to
% be on the safe side. This is done to try and NOT split the smaller data element
% into segments if only reasonable.

% We do some additional manipulation that produces integer numbers
% of segments and try to keep the smaller number of segments equal to 1
% if only possible. The code basically rounds DOWN the smaller number of
% segments (but not to less than 1!), then computes what the other number
% of segments has to be based on the available space.


if  splitProjBetweenGPUs==1 % We know we are splitting projections between GPUs at the lowest level
    noProjSeg=(2*projSize)/(N*minAvailMem);
    noVolSeg=(2*volSize)/minAvailMem;
    
    % Try rounding the smaller number of segments to 1 if reasonable, if not
    % round up, then compute the other number of segments based on how much
    % space is left.
    if noProjSeg<=noVolSeg
        noProjSeg=SmartRound(noProjSeg);  % This tries to keep number of segments at 1 if possible, otheriwse rounds up
        % We know how much space the projection data will take on GPU
        volGPU=minAvailMem-projSize/(N*noProjSeg);
        % volSize/noVolSeg = volGPU from our basic equation, so:
        noVolSeg=ceil(volSize/volGPU); % Just round up number of segments for larger element
    else
        noVolSeg=SmartRound(noVolSeg);
        % We know how much space the volume data will take on GPU, so
        projGPU=minAvailMem-volSize/noVolSeg;  % Space left for projection data
        % projSize/(N*noProjSeg)=projGPU from our basic equation, so:
        noProjSeg = ceil(projSize/(N*projGPU)); % For larger data element, we just round up number of segments
    end
    
else   % We are splitting VOLUME between the GPUs at the lowest level
    noProjSeg=(2*projSize)/minAvailMem;
    noVolSeg=(2*volSize)/(N*minAvailMem);
    
   % Try rounding the smaller number of segments to 1 if reasonable, if not
   % round up, then compute the other number of segments based on how much
   % space is left.
   if noProjSeg<=noVolSeg
        noProjSeg=SmartRound(noProjSeg);  % This tries to keep number of segments at 1 if possible, otheriwse rounds up
        % We know how much space the projection data will take on GPU
        volGPU=minAvailMem-projSize/noProjSeg;
        % volSize/(N*noVolSeg) = volGPU from our basic equation, so:
        noVolSeg=ceil(volSize/(N*volGPU)); % Just round up number of segments for larger element
   else
        noVolSeg=SmartRound(noVolSeg);
        % We know how much space the volume data will take on GPU, so
        projGPU=minAvailMem-volSize/(N*noVolSeg);  % Space left for projection data
        % projSize/noProjSeg=projGPU from our basic equation, so:
        noProjSeg = ceil(projSize/projGPU); % For larger data element, we just round up number of segments
   end
  
end  % We are splitting VOLUME between the GPUs at the lowest level

end  % FindNoOfVolAndProjSegmentsForBackprojection


%% Auxiliary Function

function roundedNoSeg=SmartRound(noOfSeg)

    if noOfSeg<=1 
        roundedNoSeg = 1;  % We can't have less than one segment
    else  
        roundedNoSeg=floor(noOfSeg);  % round down if it is larger
    end  

end  % function SmartRound


% function roundedNoSeg=SmartRound(noOfSeg)
% 
%     if noOfSeg<=1 || round(noOfSeg)==1
%         roundedNoSeg = 1;  % We can't have less than one segment, plus if it's close to 1, try to keep it at one
%     else  % noOfSeg is larger than 1.5
%         roundedNoSeg=ceil(noOfSeg);  % if it's larger than 1.5, round up
%     end  % if noOfSeg is larger than 1.5
% 
% end  % function SmartRound

