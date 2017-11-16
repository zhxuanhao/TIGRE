function [ selectedDev ] = SelectSpecificGPU(gpuNameSegment)

% This function assumes that Parallel Processing Toolbox is installed!
% Function below RESETS ALL GPUs, then selects the last (last ID)
% encountered GPU with name containing gpuNameSegment.
% Function selects (sets as current) the sepecific GPU (given by the string
% in gpuNameSegment) GPUs present in the system. If multiple GPUs of the
% same type are present, the LAST encountered will be returned.

GPUFound=0;

for i=1:gpuDeviceCount  % We iterate through ALL GPUs
    selectedDev=gpuDevice(i);    % This also RESETS the specifies GPU
    if ~isempty(strfind(selectedDev.Name, gpuNameSegment))  % We found GPU we wanted (name containing specific segment)
        GPUFound=i;  % Remember index of the device (we'll remember the last index if multiple devices of the same type are present)
    end
end  % iterating through GPUs

if GPUFound>0
    selectedDev=gpuDevice(GPUFound);
else
    selectedDev = []; % Return empty device if specific type not found
end

end  % Function

