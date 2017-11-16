function [ projections ] = AxLargeData(vol, geo, angles, noOfGPUs, gpuRAMinBytes)

% RB, May 2017

% The last two parameters are optional. If not specified, we assume we have
% Parallel Processing Toolbox and can get the number of GPUs and GPU RAM
% sizes using the toolbox functionality. Otherwise, the caller can specify
% both the number of GPUs and GPU RAM size.

% The simplifying assumption made here is that IF MULTIPLE GPUs ARE PRESENT, 
% THEY HAVE THE SAME OR VERY SIMILAR AMOUNT OF RAM.

% Another important assumption is that WE HAVE ENOUGH RAM in the PC. I
% recommend 128 GB for reconstructing 2048^3 volumes with 512 - 1000
% 2048^2 projections. MATLAB creates some auxiliary intermediate copies of
% arrays, so a significant overhead above the 32 GB is needed to store 2048^3
% volume is reqiuired. In my test using MATLAB R2017a with 2048^3 volume and 800 x 2048^2
% projections, MATLAB memory usage reached 90 GB during some test runs
% (usually it was within 75 GB range).

% The last two parameters are optional. If not specified, we assume we have
% Parallel Processing Toolbox and can figure out these parameters using the
% toolbox functions

% gpuRAMinBytes is assumed to be the NOMINAL (total) amount of GPU RAM known to the
% caller (e.g. 1.6106e+10 (1024^3*15) for 15 GB total GPU memory)


if nargin==4
    error('CBCT:AxLargeData:InvalidInput: GPU RAM in bytes needs to be specified if the number of GPUs is given!');
end

if nargin < 4
    noOfGPUs=0;  % We'll use this to mean "Parallel Processing Toolbox is present"
    gpuRAMinBytes=0;
end

disp('Computing number of data segments for forward projection...');

% Reset all GPUs (only if Parallel Processing Toolbox is present; we know
% it is present if noOfGPUs was set to 0; see the condition above)
if noOfGPUs==0   % means we have Parallel Processing Toolbox
    ResetAllGPUs();
end % If Parallel Processing Toolbox is present

% Depending on the presence/absence of the Parallel Processing Toolbox,
% call the function below either without noOfGPUs and gpuRAMinBytes, or
% with these parameters (as specified by the caller)
if noOfGPUs==0  % means we DO have Parallel Processing Toolbox
    % We pass size of the volume (in voxels) and size of a SINGLE
    % projection image
    noVolSeg = FindNoOfVolSegmentsForForwardProjection(geo.nVoxel,geo.nDetector);
else  % user-specified number of GPUs and GPU RAM size
    noVolSeg = FindNoOfVolSegmentsForForwardProjection(geo.nVoxel,geo.nDetector, noOfGPUs, gpuRAMinBytes);
end  

disp([num2str(noVolSeg) ' volume segments are required']);

% Now we know how many volume segments we have to split our data in to
% create projections, so let's do it. 
% For each volume SEGMENT, we need to:
% 1. Project it onto all projections.
% 2. Add the projection results to the output projections (accumulate
%    projections from multiple volume segments in other words)
% This is whet is done in the code below for noVolSeg > 1
% tic;

disp('Processing data segments...');

% If we have only one volume segment, we can generate projections directly
% (just invoke the original TIGRE Ax)
if noVolSeg==1
    
    disp('Forward projecting full volume...');
    projections=Ax(vol,geo,angles,'interpolated');
    
else  % We need to split volume into segments!
    
    finalVolumeZSize=geo.nVoxel(3);  % need to remember this for the loop below
    % Compute sizes of our volume segments. Assume that we will 
    % split the volume HORIZONTALLY (along Z) only for simplicity
    % Round up below in case the volume size is not divisible by number of segments 
    % (we need integer segment size!). 
    volSegZSize=ceil(finalVolumeZSize/noVolSeg);
    origVolZOffset=geo.offOrigin;  % Need to remember it, since we'll modify it for each segment
    
    % Initialize our output projection stack
    % projections=zeros(geo.nDetector(1), geo.nDetector(2), length(angles), 'single');
    % BUG FIX, RB, 07/12/2017: The rows and columns were swapped above (bug not
    % noticed because projections were square in tests). The geo.nDetector vector keeps
    % COLUMNS (x size) first, and we need to specify ROWS first when allocating our projection stack:
    projections=zeros(geo.nDetector(2), geo.nDetector(1), length(angles), 'single');

    for volSeg=1:noVolSeg
        % Two indexes below define our Z slice range that will form our volume segment
        begZIndex=(volSeg-1)*volSegZSize+1;
        endZIndex=volSeg*volSegZSize;
        % Check bounds in case total Z number of slices is not divisible by
        % noVolSeg
        if endZIndex>finalVolumeZSize
            endZIndex=finalVolumeZSize;
        end
        sizeOfSegment=endZIndex-begZIndex+1;  % Last segment size may be different than others if total size is not divisible by noVolSeg
        % Now we have to update some geo parameters, since we are processing
        % just a SEGMENT of the original volume, not the entire thing
        geo.nVoxel(3)=sizeOfSegment;  % Z size of segment in VOXELS
        geo.sVoxel(3)=sizeOfSegment*geo.dVoxel(3); % Actual Z dimension of segment (in mm)
        % The Z origin shift below is critical: our segment is NOT centered at the
        % origin of the large volume, but offset vertically!
        volSegmentZOffset=((endZIndex+begZIndex-1)/2 - finalVolumeZSize/2)*geo.dVoxel(3);
        geo.offOrigin=origVolZOffset+[0;0;volSegmentZOffset];  % We only shift volume segment in Z

        % Here we are just accumulating projection images from different volume
        % segments (everything is linear, so this is OK).
        disp(['Forward projecting volume segment ' num2str(volSeg) '...']);
        projections=projections + Ax(vol(:,:,begZIndex:endZIndex),geo,angles,'interpolated');

    end  % iterating through VOLUME segments
    
end % If we needed to split the volume into segments

disp('AxLargeData() done!');

% toc;

end  % function

