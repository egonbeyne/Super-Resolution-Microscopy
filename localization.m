% Constants
clear all;
a           = 150e-9;   % [m] Pixel size
xsize       = 4;
ysize       = 4;
sensitivity = 0.7;  %threshold value in between the average and maximum intensity (sens = [0,1])
dataLoc     = 'C:\Users\kaan_\EE\minor\Final project\matlab code\data\eye2.0\sequence\';
%GtLoc       = 'C:\Users\Egon Beyne\Desktop\Super-Resolution Microscopy\Data\Ground-truth\';
Nfiles      = 41;   %Number of datafiles

% Empty array to store all localizations
tot_loc     = [];
tot_psf     = [];

for iter = 1:Nfiles

% Open the image
imageData = OpenIm(dataLoc, iter);

% Segmenting the data
centroids = segment_frame(imageData,sensitivity);

% Localization
localizations = Fit_Gaussian(centroids, imageData, xsize, ysize, iter);

xs = localizations(:, 1)*a;
ys = localizations(:, 2)*a;

% Store localization and frame number in one big array
tot_loc = [tot_loc; [localizations(:, 1), localizations(:, 2),localizations(:,5)]];
tot_psf = [tot_psf; [localizations(:,1),localizations(:,2),localizations(:,3), localizations(:,4)]];


% Cramer-Rao lower bound calculation
N     = sum(imageData, 'all');               % [-] Number of signal photons 
sigg  = mean(localizations(3))*a;                              % Converting std. of PSF to m from pixels
sige2 = (sigg^2) + ((a^2)/12);         
tau   = 2*pi*(sigg^2)*mean(localizations(4))/(N*(a^2));    % [-] Dimensionless background parameter

% Cramer rao lower bound
dx   = sqrt(sige2*(1 + (4*tau) + sqrt(2*tau/(1 + (4*tau))))/N);

% Mortensen lower bound
dx_ls = sqrt((sigg^2 + ((a^2)/12))*((16/9) + (4*tau))/N)/1e-9;

% Compare the results to the ground truth
%comp = norm(groundtruth(localizations, GtLoc, iter));

%{
%plotting segments
imshow(uint16(imageData)*100, 'InitialMagnification', 'fit');
hold on
axis on
plot(centroids(:,1), centroids(:,2), 'b*')
plot(xs/a + 0.5, ys/a + 0.5, 'r+', 'MarkerSize', 8, 'LineWidth', 1);
bsize = 10;
N = length(centroids(:,1));
for i = 1:N
    rectangle('Position', [centroids(i,1)-bsize/2, centroids(i,2)-bsize/2, bsize, bsize],'EdgeColor','r')
    
end
hold off
%}

end

%%%%%%%%
%plotting final result
reconstruct(256,20,tot_psf)

function [localizations] = Fit_Gaussian(centroids, imageData, xsize, ysize, iter)

% Empty array to append localizations to
localizations = [];

X = linspace(1, xsize, xsize);
Y = linspace(1, ysize, ysize);

% Grid coordinates for every point in the ROI
[Xi, Yi] = meshgrid(X, Y);

% Making a 2 dimensional array for all points in the grid
c = cat(2,Xi',Yi');
d = reshape(c,[],2);

xi = d(:, 1);
yi = d(:, 2);

% Coordinates for the origin of the local coordinate system in each region
% of interest
origin = round(centroids) - [xsize/2, ysize/2];

L = size(centroids);

localizations = zeros(L(1), 5);

for i = 1:L(1)
    
    % Extract coordinates of the origin of one region of interest
    roi = origin(i, :);
       
    % Select pixel data in the region of interest
    localData = imageData(roi(2):(roi(2)+xsize-1), roi(1):(roi(1)+ysize-1));

    % rearranging ROI
    Ii  = transpose(localData);
    I = Ii(:);

    % Calculating the weights (assuming poissoning noise, based on
    % (Jiaqing,2016))
    w = max(I.^-1, .005);

    % Define a point spread function
    PSF = @(xs, ys, sg, int, b, x, y)...
            (((int/(2*pi*sg^2))*exp(-((x-xs).^2 + (y-ys).^2)/(2*sg^2)))+(b));

    % Weighted least squares
    LS = fit([xi, yi], I, PSF, 'startpoint', [3,3,1,5000, 100],... 'lower', [0, 0, 0, 0, 0], ...
        'Robust', 'LAR',...
        'Algorithm', 'Trust-Region',...
        'weight', w,...
        'TolX', [10^-2],...
        'MaxIter', 5);

    % Fitting parameters [xs, ys, sigma, Intensity]
    param = coeffvalues(LS);
    
    xs  = param(1) + roi(1) - 1.5; % +0.5 is necessary because pixel count starts at one in the local frame
    ys  = param(2) + roi(2) - 1.5; % and because the center of a pixel is where the count starts
    sg  = param(3);
    b   = param(5);
    
    localizations(i, :) = [xs, ys, sg, b, iter];
    
    
end


end

function [comp] = groundtruth(localizations, groundtruth_loc, iter)

% Format the iteration number to match the name of the data
num_gt  = '00000';
num_gt(5-strlength(string(iter)) + 1:5) = string(iter);
dir_gt = append(groundtruth_loc, num_gt ,'.csv');

% Open the csv file with ground truth
GT = readmatrix(dir_gt);

% Select the position ground truth
pos = GT(:, 3:4);

% Sort the result, such that the same points are compared
sort = sortrows(pos, 1);

% Select the position only 
pos_loc = localizations(:, 1:2);

% Sort the localization and convert to nm
loc  = sortrows(pos_loc, 1)*150;

% Get the number of datapoints present and the number found
N_pr = size(sort);
N_fd = size(loc);

% Missed detection
if N_pr(1)>N_fd(1)
    % Iterate until all missed detections are found
    for i = 1:(N_pr-N_fd)
        % Find the difference in location between datapoints. The missed
        % detection is probably where they are closest
        diff = sort(2:end, :) - sort(1:(end-1), :);

        % Total difference
        tot = sqrt((diff(:, 1).^2) + (diff(:, 2).^2));

        % Get the index of the minimum 
        [~, I] = min(tot);

        % Remove the entry
        sort(I, :) = [];
    end
end

% Calculate the squared error
rsq  = (sort - loc).^2;

L_rsq = size(rsq);

% Compute the mean squared error
comp = sum(sqrt(rsq))/L_rsq(1);


end

function imageData = OpenIm(location, iter)

% Format the iteration number to match the name of the data
num  = '00000';
num(5-strlength(string(iter)) + 1:5) = string(iter);
dir = append(location, num ,'.tif');

%'C:\Users\Egon Beyne\Desktop\Super-Resolution Microscopy\Data\sequence\00002.tif';
% Opening the raw .tif file
t = Tiff(dir, 'r');

% Reading the image
imageData = double(read(t));


close(t);

end


