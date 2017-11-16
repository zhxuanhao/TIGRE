function [ volume ] = FDKLargeData(projectionStack, geo, angles, noOfGPUs, gpuRAMinBytes)

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
    error('CBCT:FDKLargeData:InvalidInput: GPU RAM in bytes needs to be specified if the number of GPUs is given!');
end

if nargin < 4
    noOfGPUs=0;  % We'll use this to mean "Parallel Processing Toolbox is present"
    gpuRAMinBytes=0;
end

disp('Computing number of data segments for backprojection...');

% Reset all GPUs (only if Parallel Processing Toolbox is present; we know
% it is present if noOfGPUs was set to 0; see the condition above)
if noOfGPUs==0   % means we have Parallel Processing Toolbox
    ResetAllGPUs();
end % If Parallel Processing Toolbox is present

% Depending on the presence/absence of the Parallel Processing Toolbox,
% call the function below either without noOfGPUs and gpuRAMinBytes, or
% with these parameters (as specified by the caller)
if noOfGPUs==0  % means we DO have Parallel Processing Toolbox
    % We pass size of the volume (in voxels) and size of the FULL projeciton stack (in pixels)
    %[noVolSeg, noProjSeg] = FindNumberOfVolumeAndProjectionSegmentsForBackprojection(geo.nVoxel,[geo.nDetector; length(angles)])
    %[noVolSeg, noProjSeg] = FindNumberOfVolumeAndProjectionSegmentsForBackprojectionSV(geo.nVoxel,[geo.nDetector; length(angles)])
    % Last 1 below means that we are splitting PROJECTIONS between GPUs on the lowest level.
     [noVolSeg, noProjSeg] = FindNoOfVolAndProjSegmentsForBackprojection(geo.nVoxel,[geo.nDetector; length(angles)], 1);
 
else  % user-specified number of GPUs and GPU RAM size
    %[noVolSeg, noProjSeg] = FindNumberOfVolumeAndProjectionSegmentsForBackprojection(geo.nVoxel,[geo.nDetector; length(angles)], noOfGPUs, gpuRAMinBytes);
    %[noVolSeg, noProjSeg] = FindNumberOfVolumeAndProjectionSegmentsForBackprojectionSV(geo.nVoxel,[geo.nDetector; length(angles)], noOfGPUs, gpuRAMinBytes);
    [noVolSeg, noProjSeg] = FindNoOfVolAndProjSegmentsForBackprojection(geo.nVoxel,[geo.nDetector; length(angles)], 1, noOfGPUs, gpuRAMinBytes);
end  

disp([num2str(noVolSeg) ' volume segments and ' num2str(noProjSeg) ' projection segments are required']);

% Select P5000 GPU (important for fast projection filtering; otherwise we may default to some slower GPU)
SelectSpecificGPU('P5000');

%tic;
% Weigh and Pre-filter ALL projections before splitting anything
disp('Weighting and ramp filtering all projections...');
projFiltered=WeighAndFilterProjectionsForFDK(projectionStack,geo,angles);

% Now we know how many volume and projection segments we have to split our
% data into to reconstruct, so let's do it. 
% For each volume SEGMENT, we need to (conceptually, see the efficient implementation below):
% 1. Replicate it noProjSeg times.
% 2. Reconstruct each replicated volume segment copy with one of the noProjSeg segments of
%    projection set (need to apply vertical volume offset to account for
%    Z position of the volume segmemt!)
% 3. Add all replicated copies to obtain completely reconstructed volume
%    SEGMENT
% Then we combine (concatenate) all the segmemts into the final
% (large) reconstructed volume.

finalVolumeZSize=geo.nVoxel(3);

% Compute sizes of our volume and projeciton segments. Assume that we will 
% split the volume HORIZONTALLY (along Z) only for simplicity
% Round up below in case the sizes are not divisible by number of segments 
% (we need integer segment sizes!). 
volSegZSize=ceil(finalVolumeZSize/noVolSeg);
projSegSize=ceil(length(angles)/noProjSeg);
origVolZOffset=geo.offOrigin;  % Need to remember it, since we'll modify it for each segment

% Let's cover the case if we have 1 volume and 1 projection segment
% separately, to provide good throughput when we don't need to split
% anything

disp('Processing data segments...');

if noVolSeg==1 && noProjSeg==1  % We are NOT splitting anything

    disp('Backprojecting full volume and full projecion set...');
    volume = Atb(projFiltered, geo, angles);

else  % we are splitting either volume or projections, or both

    volume=zeros(geo.nVoxel(1), geo.nVoxel(2), geo.nVoxel(3), 'single');
    
    % Iterate through all volume segments
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

        % If we only have a single projeciton set segment, we just need to generate
        % our volume segment without any extra operations 
        % (e.g. no need to select subsets of projections and angles). Do it in a seperate
        % condition - it is faster than what is done for multiple projection
        % set segments, most likely due to extra copy of projections array being made
        % if the range of projection indexes is specified in call to Atb.

        if noProjSeg==1  % We aren't splitting the projection set
            disp(['Backprojecting volume segment ' num2str(volSeg) '...']);
            volume(:,:,begZIndex:endZIndex) = Atb(projFiltered, geo, angles);
        else  % More than one projection set segments need to be used

            % OK, now we need to "replicate" our volume segment so that we can partially reconstruct
            % it with each projection segment, then we will add the replicated volume segments
            % to obtain the fully reconstructed segment, and put in in the proper
            % location in the large volume. We can actually accomplish it by just
            % accumulating backprojeciton results in an auxiliary volume segment,
            % then putting it in the proper location in the large volume. Even
            % simpler (this is what I'll do below), we can accumulate backprojection results 
            % DIRECTLY in the large volume, without any need for auxiliary segment!

            for projSeg=1:noProjSeg
                % Two indexes below define our projection number range that will
                % form our projection set segment to be processed with the volume
                % segment 
                begProjIndex=(projSeg-1)*projSegSize+1;
                endProjIndex=projSeg*projSegSize;

                % Check bounds in case total number of projections is not divisible by
                % noProjSeg
                if endProjIndex>length(angles)
                    endProjIndex=length(angles);
                end

                % Accumulate our partial results DIRECTLY into the large volume
                % (everything is linear, so this is correct).    
                disp(['Backprojecting volume segment ' num2str(volSeg) ' and projection segment ' num2str(projSeg) '...']);
                volume(:,:,begZIndex:endZIndex) = volume(:,:,begZIndex:endZIndex) + Atb(projFiltered(:,:,begProjIndex:endProjIndex), geo, angles(begProjIndex:endProjIndex));

            end  % iterating through PROJECTION segments

        end  % If we have more than one projection segment

    end  % iterating through VOLUME segments
    
end % if we are splitting volume and/or projection data

disp('FDKLargeData() done!');

% toc;

end  % function

