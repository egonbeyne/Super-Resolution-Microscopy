
function [centroids] = segment_frame(imageData)
%Removing the median value for background variation
imageData = imageData - median(imageData(:));

%filtering the image with a gaussian with sigma should be 1.5~2
image_f = imgaussfilt(imageData,3/2);

%std of the data for thresholding
%sigma_f = std(image_f(:));
sigma = std(imageData(:));

%high std should have higher threshold
thr = 3*sigma;  %should be sigma_f?

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