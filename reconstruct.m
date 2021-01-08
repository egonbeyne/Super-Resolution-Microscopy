%Reconstruct an (npixels x npixels) image of gaussian blobs on a subsize grid with pixelsize:
%(nanometer) based on (locations)

function [total_image] = reconstruct(prev_image,locations,oldpixelsize,newpixelsize,npixels)

psf_pixels = 9;    %5 is the middle pixel for location of psf
xsize_psf = psf_pixels;
ysize_psf = psf_pixels;

scaling = round(oldpixelsize/newpixelsize);  %scaling to achieve new nanometer-size subpixels
% + 2*edge to prevent the psf box to exceed the boundary of the image
xsubsize = (xsize_psf);   
ysubsize = (ysize_psf);    

%can be cleaned every molecule
psf_im = zeros((npixels+psf_pixels)*scaling, (npixels+psf_pixels)*scaling);

%set equal to image from previous reconstruction
total_image = prev_image;

%make psf-block
X_psf = linspace(1, xsubsize, xsubsize);
Y_psf = linspace(1, ysubsize, ysubsize);

[Xi_psf, Yi_psf] = meshgrid(X_psf, Y_psf);
c_psf = cat(2,Xi_psf',Yi_psf');
d_psf = reshape(c_psf,[],2);

xi_psf = d_psf(:, 1);
yi_psf = d_psf(:, 2);
%

%get locations
xs_psf = locations(:,1);
ys_psf = locations(:,2);
%b_psf = locations(:,3);

%calculates the psf function for every localized molecule in locations
for i = 1:size(locations,1)   

psf = PSF(((psf_pixels+1)/2),((psf_pixels+1)/2),10/newpixelsize,1,0,xi_psf,yi_psf);

%reshape psf data into block of subpixel scale
psf_block = reshape(psf,[psf_pixels,psf_pixels]);

%put the psf block at right place in image +psf_pixels to ensure the safety edge
psf_im((round(xs_psf(i)*scaling))+1:(round(xs_psf(i)*scaling)+(psf_pixels)),...
    (round(ys_psf(i)*scaling))+1:(round(ys_psf(i)*scaling)+(psf_pixels))) = psf_block;

%total picture is formed
total_image = total_image + psf_im;

%clean sheet for new molecule
psf_im = zeros((npixels+psf_pixels)*scaling, (npixels+psf_pixels)*scaling);

end


end

function [out] = PSF(xs, ys, sg, int, b, x, y)

out = (((int/(2*pi*sg^2))*exp(-((x-xs).^2 + (y-ys).^2)/(2*sg^2)))+(b));
end