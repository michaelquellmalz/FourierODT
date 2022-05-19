%% Phase Retrieval in Diffraction Tomography
%% Setup the Optical System

m = 160;
freq = 1;
typ = 'heart1'; % 'hearts' 'heart1'

nk = 1.5*m;
mt = 1.5*m;
rM = 40;
os = OpticalSystem(m, mt, nk, rM);
os.grid_factor = 1.5;
os.oversampling = 1;
%% 
% Create a test image.

f = os.createProbe(typ);
surf(os.z1,os.z3,f,'EdgeColor','none')
view(2)
limits = [-.05,.55];
caxis(limits)
axis off; axis equal; axis tight
print(typ,  '-dpng')
colorbar
title('Ground Truth f');
axis on
fprintf('xmin=%f, xmax=%f\n',min(os.z1(:)),max(os.z1(:)))
%% Setup the Diffraction Model
% Build the diffraction model corresponding to the defined system.

dm = DiffractionModel(os);
filename = ['data_sets/' typ '_utot_m' num2str(m) '_mt' num2str(mt) '_nk' num2str(nk) '_rM' num2str(rM) '.mat'];
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
if (noiseLevel > 0)
  name = [typ '-noise' num2str(noiseLevel)];
end
noiseUniform = randn(size(utot_true)) + 1i*randn(size(utot_true));
noiseUniform = noiseUniform / norm(noiseUniform(:));
utot = utot_true + noiseLevel * noiseUniform * norm(utot_true(:));
imagesc(abs(utot))
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
newplot; surface(os.z1,os.z3,f1,'EdgeColor','none')
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
tvOpt.lambda = .1;
tvOpt.FBPDopt.maxIter = 50;
tvOpt.FBPDopt.backtrack = BacktrackYH;
tic; f1tv = tvDenoising(f1, tvOpt); toc
newplot; surface(os.z1,os.z3,real(f1tv),'EdgeColor','none')
axis off; axis equal; axis tight
caxis(limits)
crop([name '-rec-bp-tv.png'],0,0)
colorbar;
fprintf('mse: %e\npsnr: %f\nssim: %f\n', ...
    immse(f1tv,f), psnr(f1tv,f), ssim(f1tv,f));
%% 
% Test the backward transform with CG.

fprintf('Backward CG\n');
dm.backwardMethod = BackwardCGNE;
dm.backwardMethod.CGNEopt.maxIter = 5;
dm.backwardMethod.weightType = 'nonuniform';
tic; f2 = dm.backward(utot); toc
newplot; surface(os.z1,os.z3,f2,'EdgeColor','none')
axis off; axis equal; axis tight
caxis(limits)
print([name '-rec-cg.png'],  '-dpng') %default -r144
colorbar; 
fprintf('mse: %e\npsnr: %f\nssim: %f\n', ...
    immse(f2,f), psnr(f2,f), ssim(f2,f));
%% 
% Smoothing the CG reconstruction with TV.

fprintf('Backward CG & TV Denoising\n');
tvOpt = TVoptions;
tvOpt.lambda = .1;
tvOpt.FBPDopt.maxIter = 50;
tvOpt.FBPDopt.backtrack = BacktrackYH;
tic; f2tv = tvDenoising(f2, tvOpt); toc
newplot; surface(os.z1,os.z3,real(f2tv),'EdgeColor','none')
axis off; axis equal; axis tight
caxis(limits)
print([name '-rec-cg-tv.png'],  '-dpng')
colorbar;
fprintf('mse: %e\npsnr: %f\nssim: %f\n', ...
    immse(f2tv,f), psnr(f2tv,f), ssim(f2tv,f));
%% 
% Test the backward transform with FBPD and TV regularization.

tic
lambda = logspace(-3,1,20);
resi = zeros(size(lambda));
tvnorm_f = zeros(size(lambda));
dm.backwardMethod = BackwardFBPD;
dm.backwardMethod.weightType = 'nonuniform';
dm.backwardMethod.FBPDopt.maxIter = 50;
dm.backwardMethod.FBPDopt.backtrack = BacktrackYH;
for j = 1:length(lambda)
dm.backwardMethod.lambda = lambda(j);
f3 = dm.backward(utot);
f3l = dm.objectToLayer(f3);
uborn = utot - exp(1i * os.k0 * os.rM);
Fu = dm.planeToLayer(uborn);
Fu = dm.LayerRemoveWeights(Fu);

resi(j) = sum(dm.extractWeights.*abs(f3l - Fu).^2,'all');
tvnorm_f(j) = (sum(sqrt(sum(gradient(f3).^2,4)),'all'));
end
toc

loglog(resi,tvnorm_f)
writematrix([lambda(:),resi(:),tvnorm_f(:)],[name '-rec-tv-lcurve.txt'],...
  'Delimiter',' ')

dm.backwardMethod = BackwardFBPD;
dm.backwardMethod.lambda = 0.1;
dm.backwardMethod.weightType = 'nonuniform';
dm.backwardMethod.FBPDopt.maxIter = 50;
dm.backwardMethod.FBPDopt.backtrack = BacktrackYH;
tic; f3 = dm.backward(utot); toc
newplot; surface(os.z1,os.z3,f3,'EdgeColor','none')
axis off; axis equal; axis tight
caxis(limits)
print([name '-rec-tv.png'],  '-dpng')
colorbar
title('Backward PD');
fprintf('mse: %e\npsnr: %f\nssim: %f\n', ...
    immse(f3,f), psnr(f3,f), ssim(f3,f));

dm.backwardMethod = BackwardFBPD;
dm.backwardMethod.lambda = 1e-3;
dm.backwardMethod.weightType = 'nonuniform';
dm.backwardMethod.FBPDopt.maxIter = 50;
dm.backwardMethod.FBPDopt.backtrack = BacktrackYH;
tic; f3 = dm.backward(utot); toc
newplot; surface(os.z1,os.z3,f3,'EdgeColor','none')
axis off; axis equal; axis tight
caxis(limits)
print([name '-rec-tv-lambda' num2str(dm.backwardMethod.lambda) '.png'],  '-dpng')
colorbar
title(['Backward PD, lambda = ' num2str(dm.backwardMethod.lambda)]);
fprintf('mse: %e\npsnr: %f\nssim: %f\n', ...
    immse(f3,f), psnr(f3,f), ssim(f3,f));

dm.backwardMethod = BackwardFBPD;
dm.backwardMethod.lambda = 1e-2;
dm.backwardMethod.weightType = 'nonuniform';
dm.backwardMethod.FBPDopt.maxIter = 50;
dm.backwardMethod.FBPDopt.backtrack = BacktrackYH;
tic; f3 = dm.backward(utot); toc
newplot; surface(os.z1,os.z3,f3,'EdgeColor','none')
axis off; axis equal; axis tight
caxis(limits)
print([name '-rec-tv-lambda' num2str(dm.backwardMethod.lambda) '.png'],  '-dpng')
colorbar
title(['Backward PD, lambda = ' num2str(dm.backwardMethod.lambda)]);
fprintf('mse: %e\npsnr: %f\nssim: %f\n', ...
    immse(f3,f), psnr(f3,f), ssim(f3,f));

dm.backwardMethod = BackwardFBPD;
dm.backwardMethod.lambda = 1e-1;
dm.backwardMethod.weightType = 'nonuniform';
dm.backwardMethod.FBPDopt.maxIter = 50;
dm.backwardMethod.FBPDopt.backtrack = BacktrackYH;
tic; f3 = dm.backward(utot); toc
newplot; surface(os.z1,os.z3,f3,'EdgeColor','none')
axis off; axis equal; axis tight
caxis(limits)
print([name '-rec-tv-lambda' num2str(dm.backwardMethod.lambda) '.png'],  '-dpng')
colorbar
title(['Backward PD, lambda = ' num2str(dm.backwardMethod.lambda)]);
fprintf('mse: %e\npsnr: %f\nssim: %f\n', ...
    immse(f3,f), psnr(f3,f), ssim(f3,f));

dm.backwardMethod = BackwardFBPD;
dm.backwardMethod.lambda = 1e0;
dm.backwardMethod.weightType = 'nonuniform';
dm.backwardMethod.FBPDopt.maxIter = 50;
dm.backwardMethod.FBPDopt.backtrack = BacktrackYH;
tic; f3 = dm.backward(utot); toc
newplot; surface(os.z1,os.z3,f3,'EdgeColor','none')
axis off; axis equal; axis tight
caxis(limits)
print([name '-rec-tv-lambda' num2str(dm.backwardMethod.lambda) '.png'],  '-dpng')
colorbar
title(['Backward PD, lambda = ' num2str(dm.backwardMethod.lambda)]);
fprintf('mse: %e\npsnr: %f\nssim: %f\n', ...
    immse(f3,f), psnr(f3,f), ssim(f3,f));

dm.backwardMethod = BackwardFBPD;
dm.backwardMethod.lambda = 1e1;
dm.backwardMethod.weightType = 'nonuniform';
dm.backwardMethod.FBPDopt.maxIter = 50;
dm.backwardMethod.FBPDopt.backtrack = BacktrackYH;
tic; f3 = dm.backward(utot); toc
newplot; surface(os.z1,os.z3,f3,'EdgeColor','none')
axis off; axis equal; axis tight
caxis(limits)
print([name '-rec-tv-lambda' num2str(dm.backwardMethod.lambda) '.png'],  '-dpng')
colorbar
title(['Backward PD, lambda = ' num2str(dm.backwardMethod.lambda)]);
fprintf('mse: %e\npsnr: %f\nssim: %f\n', ...
    immse(f3,f), psnr(f3,f), ssim(f3,f));
%% Phase Retrieval Methods
% Take the absulute values of the total fields.

uabs = abs(utot);
%% 
% Test error reduction with CG.

fprintf('Error Reduction CG\n')
tic
dm.backwardMethod = BackwardCGNE;
dm.backwardMethod.CGNEopt.maxIter = 5;
dm.backwardMethod.weightType = 'nonuniform';
PRopt = HIOoptions;
PRopt.supportRadius = 40;
PRopt.maxIter = 10;
%PRopt.f0 = f;
tic; f4 = phaseRetrieval(dm, uabs, PRopt); toc
newplot; surface(os.z1,os.z3,f4,'EdgeColor','none')
axis off; axis equal; axis tight
caxis(limits)
print([name '-pr-cg.png'],  '-dpng')
colorbar
title('Error Reduction CG');
toc
fprintf('mse: %e\npsnr: %f\nssim: %f\n', ...
    immse(f4,f), psnr(f4,f), ssim(f4,f));
%% 
% Smoothing the CG retrieval with TV.

fprintf('Error Reduction CG & TV Denoising\n');
tvOpt = TVoptions;
tvOpt.lambda = 0.1;
tvOpt.FBPDopt.backtrack = BacktrackYH;
tvOpt.FBPDopt.maxIter = 50;
tic; f4tv = tvDenoising(f4, tvOpt); toc
newplot; surface(os.z1,os.z3,f4tv,'EdgeColor','none')
axis off; axis equal; axis tight
caxis(limits)
print([name '-pr-cg-tv.png'],  '-dpng')
colorbar
title('Error Reduction CG & TV Denoising');
fprintf('mse: %e\npsnr: %f\nssim: %f\n', ...
    immse(f4tv,f), psnr(f4tv,f), ssim(f4tv,f));
%% 
% Test error reduction with FBPD and TV regularization.

tic
lambda = logspace(-5,1,20);
resi = zeros(size(lambda));
tvnorm_f = zeros(size(lambda));
for j = 1:length(lambda)
  dm.backwardMethod = BackwardFBPD;
  dm.backwardMethod.lambda = lambda(j);
  dm.backwardMethod.FBPDopt.maxIter = 10;
  dm.backwardMethod.weightType = 'nonuniform';
  dm.backwardMethod.FBPDopt.backtrack = BacktrackYH;
  dm.backwardMethod.resume = true;
  PRopt = HIOoptions;
  PRopt.supportRadius = 40;
  PRopt.maxIter = 50;
  PRopt.f0 = f4;
  f5 = phaseRetrieval(dm, uabs, PRopt);
  
  resi(j) = sum((abs(dm.forward(f5)) - uabs).^2,'all');
  tvnorm_f(j) = sum(sqrt(sum(gradient(f5).^2,4)),'all');
end
toc

loglog(resi,tvnorm_f)
writematrix([lambda(:),resi(:),tvnorm_f(:)],[name '-pr-tv-lcurve.txt'],...
  'Delimiter',' ')
%%
fprintf('Error Reduction PD\n');
tic
dm.backwardMethod = BackwardFBPD;
dm.backwardMethod.lambda = 0.05;
dm.backwardMethod.FBPDopt.maxIter = 10;
dm.backwardMethod.weightType = 'nonuniform';
dm.backwardMethod.FBPDopt.backtrack = BacktrackYH;
dm.backwardMethod.resume = true;
PRopt = HIOoptions;
PRopt.supportRadius = 40;
PRopt.maxIter = 50;
PRopt.f0 = f4;
tic; f5 = phaseRetrieval(dm, uabs, PRopt); toc
newplot; surface(os.z1,os.z3,f5,'EdgeColor','none')
axis off; axis equal; axis tight
caxis(limits)
print([name '-pr-tv.png'],  '-dpng')
colorbar
title('Error Reduction PD');
fprintf('mse: %e\npsnr: %f\nssim: %f\n', ...
    immse(f5,f), psnr(f5,f), ssim(f5,f));
toc

for lambda2 = 10.^(-5:1)
  dm.backwardMethod = BackwardFBPD;
  dm.backwardMethod.lambda = lambda2;
  dm.backwardMethod.FBPDopt.maxIter = 10;
  dm.backwardMethod.weightType = 'nonuniform';
  dm.backwardMethod.FBPDopt.backtrack = BacktrackYH;
  dm.backwardMethod.resume = true;
  PRopt = HIOoptions;
  PRopt.supportRadius = 40;
  PRopt.maxIter = 50;
  PRopt.f0 = f4;
  f5 = phaseRetrieval(dm, uabs, PRopt);
  newplot; surface(os.z1,os.z3,f5,'EdgeColor','none')
  axis off; axis equal; axis tight
  caxis(limits)
  print([name '-pr-tv-lambda' num2str(lambda2) '.png'],  '-dpng')
  colorbar
  title('Error Reduction PD');
  fprintf('lambda = %g,  mse: %e,  psnr: %f,  ssim: %f\n  ', ...
     lambda2, immse(f5,f), psnr(f5,f), ssim(f5,f));
end
