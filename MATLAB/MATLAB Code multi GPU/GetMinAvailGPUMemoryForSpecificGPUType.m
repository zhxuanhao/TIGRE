function [ availMemory ] = GetMinAvailGPUMemoryForSpecificGPUType(gpuNameSegment)

% Function looks for specific GPUs (with names containing the string given
% in gpuNameSegment) on the system and returns the
% minimum available memory (if more than one GPU is found).
% This funciton assumes that Parallel Processing Toolbox is installed!

noOfGPUs=0;  % We will count GPU cards of specific type in this function

availMemory=1e+99;  % Set it to something very large

for i=1:gpuDeviceCount
    currDev=gpuDevice(i);
    currAvailMem=currDev.AvailableMemory;
    if ~isempty(strfind(currDev.Name, gpuNameSegment))  % We found GPU we wanted (name containing specific segment)
        noOfGPUs=noOfGPUs+1;
        if  currAvailMem<availMemory  % Only check for memory on GPUs we are interested in!
            availMemory=currAvailMem;
        end 
    end
end  % iterating through GPUs

if noOfGPUs==0  % We have not found any GPUs of type we wanted
    availMemory=0;
end

end % function

