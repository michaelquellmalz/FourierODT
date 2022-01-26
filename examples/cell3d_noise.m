%% Phase Retrieval in Diffraction Tomography
%% Setup the Optical System
% Import the data.

m = 160;
freq = 1;
typ = 'cell3d'; % 'ball'

nk = 1.5*m;
mt = 1.5*m;
rM = 40;
os = OpticalSystem(m, mt, nk, rM, 3);
os.grid_factor = 1.5;
os.oversampling = 1;
%% 
% Create a test image.

f = os.createProbe(typ);
ff = f(:,:,os.nf/2+2);
imagesc(ff.')
os.z3(1,1,os.nf/2+2)
%%
limits = [-.05,1.05];
caxis(limits)
axis off;  axis tight
%print(typ,  '-dpng')
colorbar
title('Ground Truth f');
fprintf('xmin=%f, xmax=%f\n',min(os.z1(:)),max(os.z1(:)))
%% Setup the Diffraction Model
% Build the diffraction model corresponding to the defined system. Generating 
% the 3D data with dm.forward_convolution takes more than 1 day, so we save it 
% and load it. But using only dm.forward would be some inverse crime.

dm = DiffractionModel(os);
filename = ['data_sets/' typ '_utot_m' num2str(m) '_mt' num2str(mt) '_nk' num2str(nk) '.mat'];
if exist(filename,'file')
  load(filename)
else
tic
utot_true = dm.forward_convolution(f);
toc
save(filename,'utot_true','-v7.3')
end
%%
noiseLevel = 0.05;
rng(8)
name = typ;
utot = utot_true;
if (noiseLevel > 0)
  name = [typ '-noise' num2str(noiseLevel)];
  noiseUniform = randn(size(utot_true)) + 1i*randn(size(utot_true));
  noiseUniform = noiseUniform / norm(noiseUniform(:));
  utot = utot_true + noiseLevel * noiseUniform * norm(utot_true(:));
end

imagesc(abs(utot(:,:,os.mt/4)))
axis off
print([name '-sinogram.png'],  '-dpng')
colorbar
title('Sinogram from forward trafo with Born approximation')
fprintf('point meta min=%.4f, point meta max=%.4f\n',...
  min(abs(utot(:))), max(abs(utot(:))))
%% 
% Test the backward transform with Backpropagation.

fprintf('Backpropagation\n'); 
dm.backwardMethod = BackwardCGNE;
dm.backwardMethod.CGNEopt.maxIter = 0;
dm.backwardMethod.weightType = 'nonuniform';
tic; f1 = dm.backward(utot); toc
ff = f1(:,:,os.nf/2+2);
imagesc(ff.')
axis off; axis equal; axis tight
caxis(limits)
print([name '-rec-bp.png'],  '-dpng')
colorbar;
fprintf('mse: %e\npsnr: %f\nssim: %f\n', ...
    immse(f1,f), psnr(f1,f), ssim(f1,f));
%% 
% Smoothing the backpropagation with TV.

fprintf('Backpropagation & TV Denoising\n');
tvOpt = TVoptions;
tvOpt.lambda = .05;
tvOpt.FBPDopt.maxIter = 20;
tvOpt.FBPDopt.backtrack = BacktrackYH;
tic; f1tv = tvDenoising(f1, tvOpt); toc
ff = f1tv(:,:,os.nf/2+2);
imagesc(ff.')
axis off; axis equal; axis tight
caxis(limits)
print([name '-rec-bp-tv.png'],  '-dpng')
colorbar;
fprintf('mse: %e\npsnr: %f\nssim: %f\n', ...
    immse(f1tv,f), psnr(f1tv,f), ssim(f1tv,f));
%% 
% Test the backward transform with CG.

fprintf('Backward CG\n');
dm.backwardMethod = BackwardCGNE;
dm.backwardMethod.CGNEopt.maxIter = 10;
dm.backwardMethod.weightType = 'nonuniform';
tic; f2 = dm.backward(utot); toc
ff = f2(:,:,os.nf/2+2);
imagesc(ff.')
axis off; axis equal; axis tight
caxis(limits)
print([name '-rec-cg.png'],  '-dpng')
colorbar; 
fprintf('mse: %e\npsnr: %f\nssim: %f\n', ...
    immse(f2,f), psnr(f2,f), ssim(f2,f));
%% 
% Smoothing the CG reconstruction with TV.

fprintf('Backward CG & TV Denoising\n');
tvOpt = TVoptions;
tvOpt.lambda = .05;
tvOpt.FBPDopt.maxIter = 20;
tvOpt.FBPDopt.backtrack = BacktrackYH;
tic; f2tv = tvDenoising(f2, tvOpt); toc
ff = f2tv(:,:,os.nf/2+2);
imagesc(ff.')
axis off; axis equal; axis tight
caxis(limits)
print([name '-rec-cg-tv.png'],  '-dpng')
colorbar;
fprintf('mse: %e\npsnr: %f\nssim: %f\n', ...
    immse(f2tv,f), psnr(f2tv,f), ssim(f2tv,f));
%% 
% Test the backward transform with FBPD and TV regularization.

fprintf('Backward PD\n');
dm.backwardMethod = BackwardFBPD;
dm.backwardMethod.lambda = .05;
dm.backwardMethod.FBPDopt.sigma = 1;
dm.backwardMethod.FBPDopt.tau = 1;
dm.backwardMethod.weightType = 'nonuniform';
dm.backwardMethod.FBPDopt.maxIter = 100;
dm.backwardMethod.FBPDopt.backtrack = BacktrackYH;
tic; f3 = dm.backward(utot); toc
ff = f3(:,:,os.nf/2+2);
imagesc(ff.')
axis off; axis equal; axis tight
caxis(limits)
print([name '-rec-tv.png'],  '-dpng')
colorbar
fprintf('mse: %e\npsnr: %f\nssim: %f\n', ...
    immse(f3,f), psnr(f3,f), ssim(f3,f));
%% Phase Retrieval Methods
% Take the absulute values of the total fields.

uabs = abs(utot);
%% 
% Test error reduction with CG. 

fprintf('Error Reduction CG\n');
dm.backwardMethod = BackwardCGNE;
dm.backwardMethod.CGNEopt.maxIter = 5;
dm.backwardMethod.weightType = 'nonuniform';
PRopt = HIOoptions;
PRopt.supportRadius = 40;
PRopt.maxIter = 10;
tic; f4 = phaseRetrieval(dm, uabs, PRopt); toc
ff = f4(:,:,os.nf/2+2);
imagesc(ff.')
axis off; axis equal; axis tight
caxis(limits)
print([name '-pr-cg.png'],  '-dpng')
colorbar
title('Error Reduction CG');
fprintf('mse: %e\npsnr: %f\nssim: %f\n', ...
    immse(f4,f), psnr(f4,f), ssim(f4,f));
%% 
% Smoothing the CG retrieval with TV.

fprintf('Error Reduction & TV Denoising\n');
tvOpt = TVoptions;
tvOpt.lambda = 0.05;
tvOpt.FBPDopt.backtrack = BacktrackYH;
tvOpt.FBPDopt.maxIter = 20;
tic; f4tv = tvDenoising(f4, tvOpt); toc
ff = f4tv(:,:,os.nf/2+2);
imagesc(ff.')
axis off; axis equal; axis tight
caxis(limits)
print([name '-pr-cg-tv.png'],  '-dpng')
colorbar
title('Error Reduction & TV Denoising');
fprintf('mse: %e\npsnr: %f\nssim: %f\n', ...
    immse(f4tv,f), psnr(f4tv,f), ssim(f4tv,f));
%% 
% Test error reduction with FBPD and TV regularization. 100 iterations show 
% still some artifact in the center, which is much smaller than before, but 200 
% iterations work fine.

fprintf('Error Reduction PD\n');
dm.backwardMethod = BackwardFBPD;
dm.backwardMethod.lambda = .05;
dm.backwardMethod.FBPDopt.sigma = 1;
dm.backwardMethod.FBPDopt.tau = 1;
dm.backwardMethod.FBPDopt.maxIter = 5;
dm.backwardMethod.weightType = 'nonuniform';
dm.backwardMethod.FBPDopt.backtrack = BacktrackYH;
dm.backwardMethod.resume = true;
PRopt = HIOoptions;
PRopt.supportRadius = 40;
PRopt.maxIter = 200;
PRopt.f0 = f4;
tic; f5 = phaseRetrieval(dm, uabs, PRopt); toc
ff = f5(:,:,os.nf/2+2);
imagesc(ff.')
axis off; axis equal; axis tight
caxis(limits)
print([name '-pr-tv.png'],  '-dpng')
colorbar
title('Error Reduction PD');
fprintf('mse: %e\npsnr: %f\nssim: %f\n', ...
    immse(f5,f), psnr(f5,f), ssim(f5,f));
