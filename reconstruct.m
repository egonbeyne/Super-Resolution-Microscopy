%Reconstruct an (npixels x npixels) image of gaussian blobs on a subsize grid with pixelsize:
%(nanometer) based on (locations)

function [] = reconstruct(npixels,nanometer,locations)


xsize_psf = npixels;
ysize_psf = npixels;

a = 150;    %nanometer size of original pixels 
scaling = a/nanometer;  %scaling to achieve new nano_size subgrid
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

xs_psf = locations(:,1);
ys_psf = locations(:,2);
b_psf = locations(:,3);

%calculates the psf function for every localized molecule
for i = 1:length(locations)   
    
psf = PSF(xs_psf(i)*scaling,ys_psf(i)*scaling,1,1,b_psf(i),xi_psf,yi_psf);

%stack all psf's in one image
psf_im = psf_im + reshape(psf,[xsize_psf*scaling,ysize_psf*scaling]);

end

imagesc(psf_im);
end

function [out] = PSF(xs, ys, sg, int, b, x, y)

out = (((int/(2*pi*sg^2))*exp(-((x-xs).^2 + (y-ys).^2)/(2*sg^2)))+(b));
end