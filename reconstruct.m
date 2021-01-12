%Reconstruct an (npixels x npixels) image of gaussian blobs on a subsize grid with pixelsize:
%(nanometer) based on (locations)

function [total_image] = reconstruct(prev_image,locations,oldpixelsize,newpixelsize,npixels,psf_block)

scaling = oldpixelsize/newpixelsize;
psf_pixels = 7;

%can be cleaned every molecule
psf_im = zeros((npixels+psf_pixels)*scaling, (npixels+psf_pixels)*scaling);

%set equal to image from previous reconstruction
total_image = prev_image;

%get locations
xs_psf = locations(:,1);
ys_psf = locations(:,2);
%b_psf = locations(:,3);

%add up all the molecules in one image
for i = 1:size(locations,1)   

%put the psf block at right place in image +psf_pixels to ensure the safety edge
psf_im((round(xs_psf(i)*scaling))+1:(round(xs_psf(i)*scaling)+(psf_pixels)),...
    (round(ys_psf(i)*scaling))+1:(round(ys_psf(i)*scaling)+(psf_pixels))) = psf_block;

%total picture is formed
total_image = total_image + psf_im;

%clean sheet for new molecule
psf_im = zeros((npixels+psf_pixels)*scaling, (npixels+psf_pixels)*scaling);

end


end

