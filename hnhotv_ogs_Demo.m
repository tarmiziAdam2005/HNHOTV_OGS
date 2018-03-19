clc
clear all;
close all;

imageName = 'peppers.bmp';

Img = double(imread(imageName)); %Your Image goes here

N = numel(Img);

[row, col] = size(Img);

row = int2str(row);
col = int2str(col);

imageSize = [row 'x' col];

K     =   fspecial('average',1); % For denoising. (K is the identity operator)
f = imfilter(Img,K,'circular');


sigma = 20; % Noise level ( in paper, 10, 20 and 30 are tested)

fprintf('The noise std of the observed image: %g.\n', sigma); 


f = f +  sigma * randn(size(Img)); %Add a little noise

%**************Initialize parameters for denoising*****************

opts.lam       = 0.35; % Play around with this !! 
opts.omega     = 3.0;  % and also this !!

opts.grpSz         = 3; % OGS group size
opts.Nit           = 400;
opts.Nit_inner     = 5;
opts.tol           = 1.0e-5;
opts.rho           = 0.009;  %initial rho value. rho will be updated in each iteration of the algorithm.
opts.p             = 0.1;

regType = 'ani';

out = hnhotv_ogs(f, Img, K, opts);

figure;
imshow(out.sol,[]);
title(sprintf('HNHOTV-OGS Denoised (PSNR = %3.3f dB,SSIM = %3.3f, cputime %.3f s) ',...
                       psnr_fun(out.sol,Img),ssim_index(out.sol,Img), out.cpuTime));
              
figure;
imshow(uint8(f));
title(sprintf('Noisy (PSNR = %3.3f dB, SSIM = %3.3f)',psnr_fun(f,Img), ssim_index(f,Img)));

figure;
imshow(Img,[]);
title('Original');


