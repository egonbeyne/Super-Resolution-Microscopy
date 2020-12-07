

% Opening the raw .tif image
dir = 'C:\Users\Egon Beyne\Desktop\Super-Resolution Microscopy\Data\sequence\00002.tif';
t = Tiff(dir, 'r');

% Reading the image
imageData = double(read(t));
close(t);

%threshold value in between the average and maximum intensity (sens = [0,1])
sensitivity = 0.7; 
%call segmentation function to estimate molecule locations
centroids = segment_frame(imageData,sensitivity);

N = length(centroids(:,1)); %number of molecules

%plotting segments
imshow(uint16(imageData)*100, 'InitialMagnification', 'fit');
hold on
plot(centroids(:,1), centroids(:,2), 'b*')
bsize = 10;
for i = 1:N
    rectangle('Position', [centroids(i,1)-bsize/2, centroids(i,2)-bsize/2, bsize, bsize],'EdgeColor','r')
end
hold off

