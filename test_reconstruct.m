oldpixelsize = 150;
newpixelsize = 10;
npixels = 256;
psf_pixels = 9;

%initializing empty picture

total_image = zeros((npixels+psf_pixels)*round(oldpixelsize/newpixelsize), (npixels+psf_pixels)*round(oldpixelsize/newpixelsize));

%painting new molecules when locations are given:

locations = [0,0,0.857442420686184,112.272440132808,1;  %new location measurements
            10,10,0.857442420686184,112.272440132808,1];
total_image = reconstruct(total_image,locations,oldpixelsize,newpixelsize,npixels); 

locations = [20,0,0.857442420686184,112.272440132808,1;  %new location measurements
            30,10,0.857442420686184,112.272440132808,1];
total_image = reconstruct(total_image,locations,oldpixelsize,newpixelsize,npixels);

locations = [40,0,0.857442420686184,112.272440132808,1;  %new location measurements
            256,256,0.857442420686184,112.272440132808,1];
total_image = reconstruct(total_image,locations,oldpixelsize,newpixelsize,npixels);

%show final result       
paint(total_image)

function [total_image] = reconstruct(prev_image,locations,oldpixelsize,newpixelsize,npixels)

psf_pixels = 9;    %5 is the middle pixel for location of psf
xsize_psf = psf_pixels;
ysize_psf = psf_pixels;

scaling = oldpixelsize/newpixelsize;  %scaling to achieve new nanometer-size subpixels
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

psf = PSF(((psf_pixels+1)/2),((psf_pixels+1)/2),1,1,0,xi_psf,yi_psf);

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

function [] = paint(image)
A = image;
lowestValue = min(A(A(:)>0));
highestValue = max(A(:));
imagesc(A);
cmap = jet(256);
colormap(cmap);
caxis(gca,[lowestValue-2/256, highestValue]);
% Make less than lowest value black:
cmap(1,:)=[0,0,0];
colormap(cmap)
colorbar

end
% function [total_image,psf_im,xi_psf,yi_psf,scaling] = init_image(prev_image,oldpixelsize,newpixelsize,npixels,start)
% psf_pixels = 9;    %5 is the middle pixel for location of psf
% xsize_psf = psf_pixels;
% ysize_psf = psf_pixels;
% 
% scaling = oldpixelsize/newpixelsize;  %scaling to achieve new nanometer-size subpixels
%  
% %can be cleaned every molecule
% psf_im = zeros((npixels+psf_pixels)*scaling, (npixels+psf_pixels)*scaling);
% 
% %clean sheet when starting 
% if start == 1
%   total_image = zeros((npixels+psf_pixels)*scaling, (npixels+psf_pixels)*scaling);
% else
%   total_image = prev_image;
% end
% 
% X_psf = linspace(1, xsize_psf, xsize_psf);
% Y_psf = linspace(1, ysize_psf, ysize_psf);
% 
% % Grid coordinates for every point in the ROI
% [Xi_psf, Yi_psf] = meshgrid(X_psf, Y_psf);
% c_psf = cat(2,Xi_psf',Yi_psf');
% d_psf = reshape(c_psf,[],2);
% 
% xi_psf = d_psf(:, 1);
% yi_psf = d_psf(:, 2);
% 
% end