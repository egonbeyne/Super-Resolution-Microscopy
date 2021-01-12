function psf_block = init_reconstruct(newpixelsize)

psf_pixels = 7;
xsize_psf = psf_pixels;
ysize_psf = psf_pixels;

% + 2*edge to prevent the psf box to exceed the boundary of the image
xsubsize = (xsize_psf);   
ysubsize = (ysize_psf);    
X_psf = linspace(1, xsubsize, xsubsize);
Y_psf = linspace(1, ysubsize, ysubsize);

[Xi_psf, Yi_psf] = meshgrid(X_psf, Y_psf);
c_psf = cat(2,Xi_psf',Yi_psf');
d_psf = reshape(c_psf,[],2);

xi_psf = d_psf(:, 1);
yi_psf = d_psf(:, 2);
%

psf = PSF(((psf_pixels+1)/2),((psf_pixels+1)/2),10/newpixelsize,1,0,xi_psf,yi_psf);
psf_block = reshape(psf,[psf_pixels,psf_pixels]);

function [out] = PSF(xs, ys, sg, int, b, x, y)

out = (((int/(2*pi*sg^2))*exp(-((x-xs).^2 + (y-ys).^2)/(2*sg^2)))+(b));
end

end