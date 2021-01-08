% Opening the raw .tif image
dir = 'C:\Users\kaan_\EE\minor\Final project\matlab code\data\eye2.0\sequence\00004.tif';
t = Tiff(dir,'r');

% Reading the image
imageData = double(read(t));
close(t);

%Removing the median value for background variation %%this will be expanded
%into local medians for background variations
imageData = imageData - median(imageData(:));

%filtering the image with a gaussian with sigma should be 1.5~2
image_f = imgaussfilt(imageData,3/2);

%std of the data for thresholding
sigma_f = std(image_f(:));
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

mol = length(centroids(:,1)); %number of molecules
bsize = 10;                    %size of window around molecule
 
figure(1)
 %plotting segments
 imshow(uint16(imageData)*100, 'InitialMagnification', 'fit')
 
 if (isnan(centroids))
     %do nothing
 else
   hold on
   plot(centroids(:,1),centroids(:,2), 'b*')
   for i = 1:mol
       rectangle('Position', [centroids(i,1)-bsize/2, centroids(i,2)-bsize/2, bsize, bsize],'EdgeColor','r')
   end
   hold off
 end
 
figure(2)
 %plotting segments
 imshow(uint16(image_f)*100, 'InitialMagnification', 'fit')
 if (isnan(centroids))
     %do nothing
 else
   hold on
   plot(centroids(:,1),centroids(:,2), 'b*')
   for i = 1:mol
       rectangle('Position', [centroids(i,1)-bsize/2, centroids(i,2)-bsize/2, bsize, bsize],'EdgeColor','r')
   end
   hold off
 end
 

 
% medval = median(imageData(:));           %median value in data 
%imageData =  imageData - medval;        %substracting the median value

% s = regionprops(segmented_im);

%%distance of each point to every other point
%  distances = (pdist2([x,y],[x,y]) < bsize)- eye(length(x));