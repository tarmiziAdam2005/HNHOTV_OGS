function out = hnhotv_ogs(f,Img,K,opts)

%This program solves the proposed Hybrid Non-convex Higher order 
%Overlapping group sparse Total variation image restoration.

[row, col] = size(f);
u               = f;

grpSz           = opts.grpSz; %Group size
Nit             = opts.Nit;
Nit_inner       = opts.Nit_inner;
tol             = opts.tol;
lam             = opts.lam; % The regularization parameter.
rho             = opts.rho;
omega           = opts.omega;
alpha           = 0.07;
stepLength      = 1.3; %   (gamma)
p               = opts.p;
relError        = zeros(Nit,1); % Compute error relative to the previous iteration.
psnrGain        = relError;     % PSNR improvement every iteration
ssimGain        = relError;
rhoVal        = relError;
xx            = relError;
%**************Initialize Lagrange Multipliers***********************

% Two Lagrange multipliers for the v sub-problems (The OGS term)

mu1        = zeros(row,col); % Multiplier for Dux
mu2        = zeros(row,col); % Multiplier for Duy

% Four Lagrange multipliers for the w sub-problems (The 2nd order
% non-convex term)

mu3        = zeros(row,col); %Multiplier for Duxx
mu4        = mu3;            %Multiplier for Duxy
mu5        = mu3;            %Multiplier for Duyx
mu6        = mu3;            %%Multiplier for Duyy

% One Lagrange multiplier for the z sub-problem (simple projection, box
% constraint)
mu7        = zeros(row, col);

%*************** v sub-problem variable initialization ******************

v1         = zeros(row, col); % v1 solution of the v sub-problem for Dux
v2         = v1; % v2 solutiono fhte v sub-problem for Duy

%************** w sub-problem variable initialization *******************

w1         = zeros(row, col);
w2         = w1;
w3         = w1;
w4         = w1;

%************* z sub-problem variable initialization ********************
z          = zeros(row, col);


eigK            = psf2otf(K,[row col]); %In the fourier domain
eigKtK          = abs(eigK).^2;

eigDtD  = abs(fft2([1 -1], row, col)).^2 + abs(fft2([1 -1]', row, col)).^2;
eigDDtDD = abs(psf2otf([1 -2 1],[row col])).^2 + abs(psf2otf([1 -1;-1 1],[row col])).^2 ...
            + abs(psf2otf([1 -1;-1 1],[row col])).^2 + abs(psf2otf([1;-2;1],[row col])).^2;




[D,Dt]      = defDDt(); %Declare forward finite difference operators
[DD,DDt]    = defDDt2;  %Declare 2nd order forward finite difference operators


[Dux, Duy] = D(u);
[Duxx,Duxy,Duyx,Duyy] = DD(u);

curNorm = sqrt(norm(Duxx(:) - w1(:),'fro'))^2 + sqrt(norm(Duxy(:) - w2(:),'fro'))^2 + sqrt(norm(Duyx(:) - w3(:),'fro'))^2 + sqrt(norm(Duyy(:) - w4(:),'fro'))^2;


Ktf = imfilter(f,K,'circular');
tg = tic;
    for k = 1:Nit
        
        %***************** v sub-problem (Group sparse problem)***********
        v1 = gstvdm(Dux - mu1 , grpSz , 1/rho, Nit_inner);
        v2 = gstvdm(Duy - mu2 , grpSz , 1/rho, Nit_inner);
        
        
        %***** w sub-problem (2nd order non-convex problem, IRL1 algo)**
        
        x1      = Duxx + mu3/rho;
        x2      = Duxy + mu4/rho;
        x3      = Duyx + mu5/rho;
        x4      = Duyy + mu6/rho;
        
        wgt1 = omega*p./(abs(Duxx)+ 0.00001).^(1-(p)); %IRL1 Weight update
        wgt2 = omega*p./(abs(Duxy) + 0.00001).^(1-(p));% IRL1 Weight update
        wgt3 = omega*p./(abs(Duyx) + 0.00001).^(1-(p));% IRL1 Weight update
        wgt4 = omega*p./(abs(Duyy) + 0.00001).^(1-(p));% IRL1 Weight update
               
       w1      =shrink(x1,wgt1.*omega/rho);
       w2      =shrink(x2,wgt2.*omega/rho);
       w3      =shrink(x3,wgt3.*omega/rho);
       w4      =shrink(x4,wgt4.*omega/rho);
            
        %********* u sub-problem (Least squares problem, use FFT's)********
        u_old   = u;
        %rhs  = lam*Ktf + rho*Dt(v1 + mu1, v2 + mu2) + rho*DDt(w1 - mu3/rho, w2 - mu4/rho, w3 - mu5/rho, w4 - mu6/rho) + rho*(z-mu7/rho);%hybrid
        rhs  = lam*Ktf + rho*Dt(v1 + mu1, v2 + mu2) + rho*DDt(w1 - mu3/rho, w2 - mu4/rho, w3 - mu5/rho, w4 - mu6/rho) + rho*(z-mu7);% or this (22/11/2016)
        lhs  = lam*eigKtK + rho*eigDtD + rho*eigDDtDD + rho; %hybrid
        
        
        u       = fft2(rhs)./lhs;
        u       = real(ifft2(u));
        
        [Dux, Duy] = D(u);
        [Duxx,Duxy,Duyx,Duyy] = DD(u);
        
        % *****************z sub-problem (projection)********************
        
        z = min(255,max(u + mu7,0));
        
        
        %******* mu update (Lagrange multiplier update)******************
        
        mu1     = mu1 + (v1 - Dux);
        mu2     = mu2 + (v2 - Duy);
        
        mu3     = mu3 + rho*(Duxx - w1);
        mu4     = mu4 + rho*(Duxy - w2);
        mu5     = mu5 + rho*(Duyx - w3);
        mu6     = mu6 + rho*(Duyy - w4);
        
        mu7     = mu7 + (u - z);
        
        gs1 = sqrt(conv2(abs(Dux).^2, ones(grpSz),'same'));
        gs2= sqrt(conv2(abs(Duy).^2, ones(grpSz),'same'));
            
        %***** Some statistics ***
        relError(k)    = norm(u - u_old,'fro')/norm(u, 'fro');
        psnrGain(k)    = psnr_fun(u,Img);
        ssimGain(k)    = ssim_index(u,Img);
       
         if relError(k) < tol
            break;
         end
        
         
        normOld = curNorm;
        curNorm = sqrt(norm(Duxx(:) - w1(:),'fro'))^2 + sqrt(norm(Duxy(:) - w2(:),'fro'))^2 + sqrt(norm(Duyx(:) - w3(:),'fro'))^2 + sqrt(norm(Duyy(:) - w4(:),'fro'))^2;
        rhoVal(k) = rho;
        
        if curNorm > alpha*normOld
            rho = stepLength*rho;
        end 
        
    end
    
 tg = toc(tg);
    
 out.sol                 = u;
 out.relativeError       = relError(1:k);
 out.rhoValue            = rhoVal(1:k);
 out.cpuTime             = tg;
 out.psnrGain            = psnrGain(1:k);
 out.ssimGain            = ssimGain(1:k);
 out.psnrRes             = psnr_fun(u, Img);
 out.ssimRes             = ssim_index(u, Img);
 out.OverallItration     = size(out.relativeError,2); %No of itr to converge
 out.xxVal               = xx(1:k);
end


function [D,Dt] = defDDt()
D  = @(U) ForwardDiff(U);
Dt = @(X,Y) Dive(X,Y);
end

function [Dux,Duy] = ForwardDiff(U)
 Dux = [diff(U,1,2), U(:,1,:) - U(:,end,:)];
 Duy = [diff(U,1,1); U(1,:,:) - U(end,:,:)];
end

function DtXY = Dive(X,Y)
  % Transpose of the forward finite difference operator
  % is the divergence fo the forward finite difference operator
  DtXY = [X(:,end) - X(:, 1), -diff(X,1,2)];
  DtXY = DtXY + [Y(end,:) - Y(1, :); -diff(Y,1,1)];   
end

function [DD,DDt] = defDDt2
        % defines finite difference operator D
        % and its transpose operator
        DD  = @(U) ForwardD2(U);
        DDt = @(Duxx,Duxy,Duyx,Duyy) Dive2(Duxx,Duxy,Duyx,Duyy);
 end
    
  function [Duxx Duxy Duyx Duyy] = ForwardD2(U)
        %
        Duxx = [U(:,end) - 2*U(:,1) + U(:,2), diff(U,2,2), U(:,end-1) - 2*U(:,end) + U(:,1)];
        Duyy = [U(end,:) - 2*U(1,:) + U(2,:); diff(U,2,1); U(end-1,:) - 2*U(end,:) + U(1,:)];
        %
        Aforward = U(1:end-1, 1:end-1) - U(  2:end,1:end-1) - U(1:end-1,2:end) + U(2:end,2:end);
        Bforward = U(    end, 1:end-1) - U(      1,1:end-1) - U(    end,2:end) + U(    1,2:end);
        Cforward = U(1:end-1,     end) - U(1:end-1,      1) - U(  2:end,  end) + U(2:end,    1);
        Dforward = U(    end,     end) - U(      1,    end) - U(    end,    1) + U(    1,    1);
        % 
        Eforward = [Aforward ; Bforward]; Fforward = [Cforward ; Dforward];
        Duxy = [Eforward, Fforward]; Duyx = Duxy;
        %
  end
    
   function Dt2XY = Dive2(Duxx,Duxy,Duyx,Duyy)
        %
        Dt2XY =         [Duxx(:,end) - 2*Duxx(:,1) + Duxx(:,2), diff(Duxx,2,2), Duxx(:,end-1) - 2*Duxx(:,end) + Duxx(:,1)]; % xx
        Dt2XY = Dt2XY + [Duyy(end,:) - 2*Duyy(1,:) + Duyy(2,:); diff(Duyy,2,1); Duyy(end-1,:) - 2*Duyy(end,:) + Duyy(1,:)]; % yy
        %
        Axy = Duxy(1    ,    1) - Duxy(      1,    end) - Duxy(    end,    1) + Duxy(    end,    end);
        Bxy = Duxy(1    ,2:end) - Duxy(      1,1:end-1) - Duxy(    end,2:end) + Duxy(    end,1:end-1);
        Cxy = Duxy(2:end,    1) - Duxy(1:end-1,      1) - Duxy(  2:end,  end) + Duxy(1:end-1,    end);
        Dxy = Duxy(2:end,2:end) - Duxy(  2:end,1:end-1) - Duxy(1:end-1,2:end) + Duxy(1:end-1,1:end-1);
        Exy = [Axy, Bxy]; Fxy = [Cxy, Dxy];
        %
        Dt2XY = Dt2XY + [Exy; Fxy];
        Dt2XY = Dt2XY + [Exy; Fxy];
   end
   
function z = shrink(x,r)
z = sign(x).*max(abs(x)- r,0);
end

