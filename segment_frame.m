
function [centroids] = segment_frame(imageData,oldpixelsize)

b = 1200/oldpixelsize;  %block size in pixels for removing noise
gauss_sigma = 100/oldpixelsize;

%Removing the median value locally on a [b x b] grid for noise elimination
fun = @(block_struct) block_struct.data - median(block_struct.data(:));

data = blockproc(imageData,[b,b],fun);

%filtering the image with a gaussian with sigma should be 1.5~2
image_f = imgaussfilt(data,gauss_sigma);

%std of the data for thresholding, high std should have higher threshold
sigma = std(data(:));
thr = 3*sigma; 

% thresholding pixel values based on median
segmented_im = image_f.*(image_f>thr);

%finding regional maxima as locations of molecules
if (sum(any(segmented_im)) ~= 0)
  regional_max = imregionalmax(segmented_im);
  [y,x] = find(regional_max);
  centroids = [x,y];
else
   centroids = NaN(1,2);    %no molecules detected
end

end