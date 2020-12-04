
function [centroids] = segment_frame(imdata,sens)
R = size(imdata,1);
C = size(imdata,2);
segmented_im = zeros(R,C);

%determine threshold based on mean and max value
meanval = mean(imdata(:));                      %mean value in data
maxval = max(imdata(:));                        %max value in data
thr = meanval + (1-sens)*(maxval - meanval);    %threshold being somewhere between average and max intensity

%thresholding pixel values
for i = 1:R
    for j = 1:C
        if imdata(i,j) > thr
             segmented_im(i,j) = imdata(i,j);
        else
             segmented_im(i,j) = 0;
        end
    end
end

%finding centre of molecules (logical values)
s = regionprops(logical(segmented_im));         %estimation of the centre of a connected region
centroids = cat(1, s.Centroid);                 %array of molecule locations

end