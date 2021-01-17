%Reconstruct an (npixels x npixels) image of gaussian blobs on a subsize grid with pixelsize:
%(nanometer) based on (locations)

function [total_image] = reconstruct(prev_image,locations,oldpixelsize,newpixelsize,npixels,dx)

scaling = oldpixelsize/newpixelsize;
psf_pixels = 7;
xsize_psf = psf_pixels;
ysize_psf = psf_pixels;

xsubsize = (xsize_psf);   
ysubsize = (ysize_psf);    
X_psf = linspace(1, xsubsize, xsubsize);
Y_psf = linspace(1, ysubsize, ysubsize);

[Xi_psf, Yi_psf] = meshgrid(X_psf, Y_psf);
c_psf = cat(2,Xi_psf',Yi_psf');
d_psf = reshape(c_psf,[],2);

xi_psf = d_psf(:, 1);
yi_psf = d_psf(:, 2);

dx = abs(dx);
LB = 0.3; %to set a minimum for dx
if (dx < LB)
    dx = LB;
end

psf = PSF(((psf_pixels+1)/2),((psf_pixels+1)/2),dx*10/newpixelsize,1,0,xi_psf,yi_psf);
psf_block = reshape(psf,[psf_pixels,psf_pixels]);

function [out] = PSF(xs, ys, sg, int, b, x, y)

out = (((int/(2*pi*sg^2))*exp(-((x-xs).^2 + (y-ys).^2)/(2*sg^2)))+(b));
end

%can be cleaned every molecule
psf_im = zeros((npixels+psf_pixels+2)*scaling, (npixels+psf_pixels+2)*scaling);

%set equal to image from previous reconstruction
total_image = prev_image;

%get locations
xs_psf = locations(:,2);
ys_psf = locations(:,1);
%b_psf = locations(:,3);

%add up all the molecules in one image
for i = 1:size(locations,1)   

%put the psf block at right place in image +psf_pixels to ensure the safety edge
psf_im((round((xs_psf(i)+1)*scaling)+1):(round((xs_psf(i)+1)*scaling)+(psf_pixels)),...
    (round((ys_psf(i)+1)*scaling)+1):(round((ys_psf(i)+1)*scaling)+(psf_pixels))) = psf_block;

%total picture is formed
total_image = total_image + psf_im;

%clean sheet for new molecule
psf_im = zeros((npixels+psf_pixels+2)*scaling, (npixels+psf_pixels+2)*scaling);

end


end
