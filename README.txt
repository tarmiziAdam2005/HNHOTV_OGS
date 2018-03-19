
Ver 1.1a.



This is a short note on the usage of HNHOTV-OGS paper published in the following journal:

-> Multidimensional Systems and Signal Processing
-> https://link.springer.com/article/10.1007/s11045-018-0567-3




Created by Tarmizi Adam on 15 March 2018. Copyright 2018 by Tarmizi Adam. Permission is granted herein for use in academic purpose only.

***************************************************************


1. Introduction

-This is a simple document on how to run denoising(codes can also be used for deblurring) experiments on images corrupted with additive Gaussian noise.

2. Instructions

-Make sure all the related ".m" files are in one folder. The main solver for HNHOTV-OGS algoritm is "hnhotv_ogs.m". Inside the main solver, the overlapping group sparse (OGS) solver "gstvdm.m" is called.

To start a simple denoising demo, run the matlab script "hnhotv_ogs_Demo.m".

3. Some guides

- In the paper, we tested 8 images and three noise levels. For each image and its corresponding noise levels, we manually find the best regularization parameter \lambda and \omega to ensure the best restoration quality. Then each image is denoised ten times and the average is taken.

- To somewhat simplify the user, we include just several values of the regularization parameter \lambda and \omega here for ease. The user is still encouraged to manually fine tune the these parameters and play around with them.


***Noise Level 10***

butterfly2 (450 x 450):        \lambda = 0.85, \omega = 2.5.
Einstein   (256 x 256):        \lambda = 0.85, \omega = 1.8.

***Noise Level 20***

peppers1   (256 x 256):        \lambda = 0.35, \omega = 3.0.
Boats      (512 x 512):        \lambda = 0.35, \omega = 1.5.

***Noise Level 30***

peppers1   (256 x 256):           \lambda = 0.23, \omega = 3.8.
butterfly2 (450 x 450):           \lambda = 0.27, \omega = 6.0.


4. Citations

Article paper of this work is at: 

https://link.springer.com/article/10.1007/s11045-018-0567-3


If you happen to use our code, please cite the work accordingly Thank you.
     








