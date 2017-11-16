function [noVolSeg] = FindNoOfVolSegmentsForForwardProjection(volDim, projDim, noOfGPUs, gpuRAMinBytes)
% RB, May 2017.

% Function computes the number of necessary volume
% segments needed to forward project the data on the GPU (in other words:
% segments that fit into the GPU memory).

% volDim  - dimensions of the full volume to be reconstructed in voxels (vector of three elements)

% projDim - dimensions of ONE projection image (2D) in pixels (vector of TWO elements)
% No harm will be done if a 3D vector is passed (including number of
% projections), because this function uses only the first two elements of
% the vector anyway.

% The last two parameters are optional. If not specified, we assume we have
% Parallel Processing Toolbox and can get the number of GPUs and GPU RAM
% sizes using the toolbox functionality. Otherwise, the caller can specify
% both the number of GPUs and GPU RAM size.

% gpuRAMinBytes is assumed to be the NOMINAL (total) amount of GPU RAM known to the
% caller (e.g. 1.6106e+10 (1024^3*15) for 15 GB total GPU memory)

% NOTE: We do not really need noOfGPUs parameter in this function, but I
% leave it to keep the parameter list similar to
% FindNoOfVolAndProjSegmentsForBackprojection(). This parameter is also
% used as an indicator of presence or absence of Parallel Processing Toolbox.

% The simplifying assumption made here is that IF MULTIPLE GPUs ARE PRESENT, 
% THEY HAVE THE SAME OR VERY SIMILAR AMOUNT OF RAM.

noVolSeg=1;

% The last two parameters are optional. If not specified, we assume we have
% Parallel Processing Toolbox

if nargin==3
    error('CBCT:FindNoOfVolSegmentsForForwardProjection:InvalidInput: GPU RAM in bytes needs to be specified if the number of GPUs is given!');
end

if nargin < 3
    noOfGPUs=0;
    gpuRAMinBytes=0;
end

% Compute sizes of the full volume to be reconstructed and the
% SINGLE projection IN BYTES. We compute them from the volume and
% projeciton dimensions and knowing that
% we have 4 bytes per voxel/projection pixel (single).
volSize=volDim(1)*volDim(2)*volDim(3)*4;  % Singles are 4 bytes per element
projSize=projDim(1)*projDim(2)*4;

% If noOfGPUs was specified as 0, we assume we have Parallel Processing Toolbox
% and want to figure out how many we have using the toolbox
% Otherwise we just use the passed nominal values and ASSUME SOME SAFETY
% MARGIN

if noOfGPUs>0  % Use user-specified number of GPUs and GPU RAM size
    N = noOfGPUs; % N is not used here as of May 2017, but may be in the future
    disp([num2str(N) ' GPUs specified for use']);
    % We assume the user is passing the NOMINAL GPU RAM amnount. We know that
    % in reality less RAM is available, so assume some safety margin
    minAvailMem=0.75*gpuRAMinBytes;
    
else  % Figure things out using Parallel Processing toolbox

    % As of May 2017 (prototype v.1.0) we only look for Quadro P500 to
    % simplify things
    N = GetNumberOfSpecificGPUs('P5000');  % N is not used here as of May 2017, but may be in the future
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

% Compute into how many segments we need to split the data to fit it into
% the GPU memory. Unlike in
% FindNumberOfVolumeAndProjectionSegmentsForBackprojection(), we just use a
% simple formula, since we only have to take care of splitting the volume.

% Note that, unlike in
% FindNoOfVolAndProjSegmentsForBackprojection(), we do not
% consider number of GPUs N below. This is because during forward
% projection, we simply duplicate the volume (segment) between the multiple
% GPUs and each GPU processess a subset of projection images in parallel with other GPUs.
% However, at any given time, each GPU still holds the entire volume (segment) and only a
% single projeciton image when generating the projections.

noVolSeg=ceil((volSize+projSize)/minAvailMem);

end

