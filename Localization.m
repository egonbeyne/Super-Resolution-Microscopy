% Opening the raw .tif image
dir = 'C:\Users\Egon Beyne\Desktop\Super-Resolution Microscopy\Data\sequence\00001.tif';
t = Tiff(dir, 'r');

% Reading the image
imageData = double(read(t));

% Performing a linear stretch
stretch = imageData;

% Manually selecting a region of interest
ROI  = stretch(136:145, 108:115);

% Make a grid of the local coordinates in the region of interest
[ylen, xlen] = size(ROI);

X = linspace(1, xlen, xlen);
Y = linspace(1, ylen, ylen);

% Grid coordinates for every point in the ROI
[Xi, Yi] = meshgrid(X, Y);

% Making a 2 dimensional array for all points in the grid
c = cat(2,Xi',Yi');
d = reshape(c,[],2);

xi = d(:, 1);
yi = d(:, 2);

% rearranging ROI
Ii  = transpose(ROI);
I = Ii(:);

% Calculating the weights (assuming poissoning noise, based on (Jiaqing,2016)
w = max(I.^-1, .005);

% Define a point spread function
PSF = @(xs, ys, sg, int, b, x, y)...
        (((int/(2*pi*sg^2))*exp(-((x-xs).^2 + (y-ys).^2)/(2*sg^2)))+b);

% Weighted least squares
LS = fit([xi, yi], I, PSF, 'startpoint', [5,5,1,5000, 1], 'lower', [0, 0, 0, 0, 0], ...
    'Robust', 'LAR',...
    'Algorithm', 'Trust-Region',...
    'weight', w,...
    'TolX', 10^-5,...
    'MaxIter', 100);

% Fitting parameters [xs, ys, sigma, Intensity]
param = coeffvalues(LS);

% 95% confidence intervals of estimated parameters
conf = confint(LS, 0.95);
tol  = conf(2, :) - conf(1, :);



% Cramer-Rao lower bound calculation
a     = 150e-9;                                  % [m] Pixel size
N     = sum(ROI, 'all');%97220.38;               % [-] Number of signal photons (here per molecule, not sure  if this is correct)
sigg  = param(3)*a;                              % Converting std. of PSF to m from pixels
sige2 = (sigg^2) + (a^2/12);         
tau   = 2*pi*(sigg^2)*param(5)/(N*(a^2));    % [-] Dimensionless background parameter

% Cramer rao lower bound
dx   = sqrt(sige2*(1 + (4*tau) + sqrt(2*tau/(1 + (4*tau))))/N);

tolx = tol(1)*a;    %[nm] x-tolerance 
toly = tol(2)*a;    %[nm] y-tolerance

tot  = sqrt(tolx^2 + toly^2);

% Calculate how close the estimation is to the cramer-rao bound
prec = tot/dx;

% Plot the ROI and estimated location
imshow(uint16(ROI)*100, 'InitialMagnification', 'fit');
axis on
hold on
plot(param(1) , param(2) , 'r+', 'MarkerSize', 8, 'LineWidth', 1);


close(t);