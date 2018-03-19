function psnr = psnr_fun(x, y)

if nargin<3
    m1 = max( abs(x(:)) );
    m2 = max( abs(y(:)) );
    
    if min(m1,m2)>50 % revised by LJ
        vmax = 255;
    else
        vmax = 1;
    end
end

mse= mean(abs(x(:) - y(:)).^2);
psnr = 10 * log10 (vmax^2/mse );