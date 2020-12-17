function [] = reconstruct(npixels,scaling,loc_data)

xsize_psf = npixels;
ysize_psf = npixels;

xsubsize = xsize_psf*scaling;
ysubsize = ysize_psf*scaling;

psf_im = [zeros(npixels*scaling, npixels*scaling)];

X_psf = linspace(1, xsubsize, xsubsize);
Y_psf = linspace(1, ysubsize, ysubsize);

% Grid coordinates for every point in the ROI
[Xi_psf, Yi_psf] = meshgrid(X_psf, Y_psf);
c_psf = cat(2,Xi_psf',Yi_psf');
d_psf = reshape(c_psf,[],2);

xi_psf = d_psf(:, 1);
yi_psf = d_psf(:, 2);

xs_psf = loc_data(:,1);
ys_psf = loc_data(:,2);
sg_psf = loc_data(:,3);
b_psf = loc_data(:,4);

%calculates the psf function for every localized molecule
for i = 1:length(loc_data)   
psf = PSF(xs_psf(i)*scaling,ys_psf(i)*scaling,sg_psf(i)*scaling/3,10,b_psf(i)*scaling,xi_psf,yi_psf);
psf_im = psf_im + reshape(psf,[xsize_psf*scaling,ysize_psf*scaling]);

end

imagesc(psf_im);
end

function [out] = PSF(xs, ys, sg, int, b, x, y)

out = (((int/(2*pi*sg^2))*exp(-((x-xs).^2 + (y-ys).^2)/(2*sg^2)))+(b));
end