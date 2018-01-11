function [proj,angles]=loadNikkonData(fpath,angles,nsamples,whitelevel)

%% compute index of desired data
maxdata=length(angles);
steps=floor(maxdata/nsamples);
idx=1:steps:steps*nsamples;

%% get filename
% assuming TIF and 4 digits.

myfiles = dir(fullfile([fpath,'\', '*.tif']));
%% load images
proj=[];
for ii=1:nsamples
     proj(:,:,ii)=single(imread([fpath, '\' myfiles(idx(ii)).name]));
end
%% Beer lambert
proj=-log(proj/single(whitelevel));
angles=angles(idx);