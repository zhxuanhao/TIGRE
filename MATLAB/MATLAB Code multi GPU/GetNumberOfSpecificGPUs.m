function [ noOfGPUs ] = GetNumberOfSpecificGPUs(gpuNameSegment)

% This function assumes that Parallel Processing Toolbox is installed!
% Function returns number of GPUs of specific type (given by the string
% in gpuNameSegment) GPUs present in the system

noOfGPUs=0;  % We will count specific GPU cards in this function

for i=1:gpuDeviceCount
    currDev=gpuDevice(i);
    if ~isempty(strfind(currDev.Name, gpuNameSegment))  % We found GPU we wanted (name containing specific segment)
        noOfGPUs=noOfGPUs+1;
    end
end  % iterating through GPUs

end  % Function

