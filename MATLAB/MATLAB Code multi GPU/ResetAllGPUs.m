function ResetAllGPUs()

% Function resets all GPUS present in the system.
% NOTE: Function assumes that the Parallel Processin Toolbox is installed!

for i=1:gpuDeviceCount
    reset(gpuDevice(i));
end  % for


end  % Function

