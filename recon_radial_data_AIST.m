%addpath('fessler\irt')
%addpath('fessler\irt\systems')

cd 'fessler\irt';
run setup.m
cd ..\..

rawdatfile = 'E:\ResearchMRI\data\meas_MID06716_FID05964_SAG_MPnRAGE_RADIAL_GOLDEN.dat';

k1 = mapVBVD_v2(rawdatfile);
%k1 = mapVBVD(rawdatfile);

params.nCol = k1.image.NCol; % 384
params.nSpokes = k1.image.NLin; % 45120
params.nCha = k1.image.NCha; % 20
params.baseResolution = k1.hdr.MeasYaps.sKSpace.lBaseResolution; % 192
params.fovx_mm = k1.hdr.MeasYaps.sSliceArray.asSlice{1}.dReadoutFOV*2/2; % 192

N = [params.baseResolution, params.baseResolution, params.baseResolution];

% -------------------------------------------------------------------------
% Create trajectory
% -------------------------------------------------------------------------
[trajRAD_mm, trajRAD, polarAngle, azimuthalAngle] = calcRadialTrajGA3D(params);

k_rad = reshape(trajRAD, params.nCol, params.nSpokes, 3);

k_rad = permute(k_rad,[1 3 2]);
% -------------------------------------------------------------------------
% Create NUFFT object
% -------------------------------------------------------------------------
nufft_args = {N, [4 4 4], 2*N, N/2, 'table', 2^10, 'minmax:kb'}; 

mask = true(N);

f.basis = {'dirac'};

Gm_rad = Gmri(trajRAD_mm, mask, 'fov', params.fovx_mm, 'basis', f.basis, 'nufft', nufft_args);

% -------------------------------------------------------------------------
% Calculate DCF weights 
% -------------------------------------------------------------------------
[wi, actmaxerr, psf_final] = genNufftWeightsPipe(2*pi*trajRAD, N, 0.1); % trajRAD should be scaled from -pi to pi 

% -------------------------------------------------------------------------
% Read image data
% -------------------------------------------------------------------------
kspace = k1.image.unsorted();

kspace = permute(kspace, [1 3 2]); % reorder to [COL SPOKES CHA] 

% Do the reconstruction
im_rad = zeros([N, params.nCha]);

for iC = 1:params.nCha

   fprintf('Reconstructing channel %d of %d...\n', iC, params.nCha);
    
    ima_rad = Gm_rad'*(vec(kspace(:,:,iC)).*wi);
  
    im_rad(:,:,:,iC) = reshape(ima_rad, N);
    
end
    
% -------------------------------------------------------------------------
% Combine channels
% -------------------------------------------------------------------------
im_rad_sos = sos(im_rad);
