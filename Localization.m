% Opening the raw .tif image
dir = 'C:\Users\Egon Beyne\Desktop\Super-Resolution Microscopy\Data\sequence\00002.tif';
t = Tiff(dir, 'r');

% Reading the image
imageData = double(read(t));

% Performing a linear stretch
stretch = imageData;

% Manually selecting a region of interest
ROI  = stretch(136:145, 108:115);

% Localization
localizations = Fit_Gaussian(centroids, imageData);

a     = 150e-9;                                  % [m] Pixel size

xs = localizations(:, 1)*a;
ys = localizations(:, 2)*a;

% Cramer-Rao lower bound calculation

N     = sum(imageData, 'all');               % [-] Number of signal photons (here per molecule, not sure  if this is correct)
sigg  = mean(localizations(3))*a;                              % Converting std. of PSF to m from pixels
sige2 = (sigg^2) + (a^2/12);         
tau   = 2*pi*(sigg^2)*mean(localizations(4))/(N*(a^2));    % [-] Dimensionless background parameter

% Cramer rao lower bound
dx   = sqrt(sige2*(1 + (4*tau) + sqrt(2*tau/(1 + (4*tau))))/N);

% Compare the results to the ground truth
comp = groundtruth(localizations)

% Plot the ROI and estimated location
%imshow(uint16(ROI)*100, 'InitialMagnification', 'fit');
%axis on
%hold on
%plot(param(1) , param(2) , 'r+', 'MarkerSize', 8, 'LineWidth', 1);

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


close(t);


function [localizations] = Fit_Gaussian(centroids, imageData)

% Empty array to append localizations to
localizations = [];

% Grid with local coordinates inside the ROI
xlen = 4;
ylen = 4;

X = linspace(1, xlen, xlen);
Y = linspace(1, ylen, ylen);

% Grid coordinates for every point in the ROI
[Xi, Yi] = meshgrid(X, Y);

% Making a 2 dimensional array for all points in the grid
c = cat(2,Xi',Yi');
d = reshape(c,[],2);

xi = d(:, 1);
yi = d(:, 2);

% Coordinates for the origin of the local coordinate system in each region
% of interest
origin = round(centroids) - [xlen/2, ylen/2];

L = size(centroids);

localizations = zeros(L(1), 4);

for i = 1:L(1)
    
    % Extract coordinates of the origin of one region of interest
    roi = origin(i, :);
       
    % Select pixel data in the region of interest
    localData = imageData(roi(2):(roi(2)+xlen-1), roi(1):(roi(1)+ylen-1));

    % rearranging ROI
    Ii  = transpose(localData);
    I = Ii(:);

    % Calculating the weights (assuming poissoning noise, based on
    % (Jiaqing,2016))
    w = max(I.^-1, .005);

    % Define a point spread function
    PSF = @(xs, ys, sg, int, b, x, y)...
            (((int/(2*pi*sg^2))*exp(-((x-xs).^2 + (y-ys).^2)/(2*sg^2)))+b);

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
    
    localizations(i, :) = [xs, ys, sg, b];
    

end


end

function [comp] = groundtruth(localizations)
% Open the csv file with ground truth
GT = readmatrix('C:\Users\Egon Beyne\Desktop\Super-Resolution Microscopy\Data\Ground-truth\00002.csv');

% Select the position ground truth
pos = GT(:, 3:4);

% Sort the result, such that the same points are compared
sort = sortrows(pos, 1);

% Select the position only 
pos_loc = localizations(:, 1:2);

% Sort the localization and convert to nm
loc  = sortrows(pos_loc, 1)*150;

% Calculate the squared error
rsq  = (sort - loc).^2;

L_rsq = size(rsq)

% Compute the mean squared error
comp = sum(sqrt(rsq))/L_rsq(1);


end





