% Constants
clear all;
tic
a           = 150e-9;   % [m] Pixel size
xsize       = 4;
ysize       = 4;
sensitivity = 0.7;  %threshold value in between the average and maximum intensity (sens = [0,1])
dataLoc     = 'C:\Users\Egon Beyne\MATLAB\Projects\Super-resolution\Data\sequence\';
GtLoc       = 'C:\Users\Egon Beyne\Desktop\Super-Resolution Microscopy\Data\Ground-truth\';
Nfiles      = 41;   %Number of datafiles

% Making a point spread function
psf = @(xs, ys, sg, int, b, x, y)((int/(2*pi*(sg^2)))*exp(-((x-xs).^2 + (y-ys).^2)/(2*(sg^2)))+b); 


% Convert derivatives to normal functions for speed
dfdxs = @(xs, ys, sg, int, b, x, y)(int*exp(-((x - xs).^2 + (y - ys).^2)/(2*sg^2)).*(2*x - 2*xs))/(4*sg^4*pi);
dfdys = @(xs, ys, sg, int, b, x, y)(int*exp(-((x - xs).^2 + (y - ys).^2)/(2*sg^2)).*(2*y - 2*ys))/(4*sg^4*pi);
dfdsg = @(xs, ys, sg, int, b, x, y)((int*exp(-((x - xs).^2 + (y - ys).^2)/(2*sg^2)).*((x - xs).^2 + (y - ys).^2))/(2*sg^5*pi) - (int*exp(-((x - xs).^2 + (y - ys).^2)/(2*sg^2))./(sg^3*pi)));
dfdin = @(xs, ys, sg, int, b, x, y)exp(-((x - xs).^2 + (y - ys).^2)/(2*sg^2))./(2*sg^2*pi);
dfdb  = @(xs, ys, sg, int, b, x, y)(1);

% Empty array to store all localizations
loc_molecules = []; 

% Array to store accuracy of localization in each frame
acc     = zeros(Nfiles, 5);

for iter = 1:Nfiles

% Open the image
imageData = OpenIm(dataLoc, iter);

% Segmenting the data
centroids = segment_frame(imageData,sensitivity);

% Localization
localizations = Fit_Gaussian(centroids, imageData, xsize, ysize, iter, a,dfdxs, dfdys, dfdsg, dfdin, psf);

xc = localizations(:, 1)*a;
yc = localizations(:, 2)*a;

% Store localization and frame number in one big array
loc_molecules = [loc_molecules; [localizations]];

% Cramer-Rao lower bound calculation
N     = sum(imageData, 'all');                            % [-] Number of signal photons 
sigg  = mean(localizations(3))*a;                         % Converting std. of PSF to m from pixels
sige2 = (sigg^2) + ((a^2)/12);         
tau   = 2*pi*(sige2)*mean(localizations(4))/(N*(a^2));    % [-] Dimensionless background parameter

% Cramer rao lower bound
dx   = sqrt(sige2*(1 + (4*tau) + sqrt(2*tau/(1 + (4*tau))))/N)/1e-9;

% Mortensen lower bound
dx_ls = sqrt((sigg^2 + ((a^2)/12))*((16/9) + (4*tau))/N)/1e-9;

% Compare the results to the ground truth
[comp, missed, Nmol] = groundtruth(localizations, GtLoc, iter);


% Store the CRLB and comparison to the ground truth
acc(iter, :) = [dx, comp, missed, Nmol, iter];

% Print progress
progress = string(iter/Nfiles*100) + '% done';
disp(progress)

end

% Average accuracy relative to CRLB
rel_acc = mean((acc(:, 2)./acc(:, 1))*Nmol);


% Total number of missed detections
tot_miss = sum(acc(:, 3));

% Relative accuracy for each frame
bar(acc(:, 2)./acc(:, 1))

%%%%%%%%
%plotting final result
npixels = 256;
nanometer = 5;
%reconstruct(npixels,nanometer,loc_molecules)
% Faster plotting for testing purpose
plot(loc_molecules(:, 1), loc_molecules(:, 2), 'r.')
set(gca,'Color','k')

toc

function [localizations] = Fit_Gaussian(centroids, imageData, xsize, ysize, iter, a,...
    dfdxs, dfdys, dfdsg, dfdin, psf)

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

    fi = Ii(:);


% Initial guess for function parameters
B     = [2;2;1;5000; 100];





step = [4, 4];
maxIt = 30;

it = 0;

% Damping for levenberg marquardt (Has to be tuned)
L_0 = 4.5;
v   = 1.9;

while (step(1)^2 + step(2)^2)>0.0001 && it<maxIt

    % Counting iterations
    it = it + 1;
    
    % Unpack values from b
    xs_0  = B(1);
    ys_0  = B(2);
    sg_0  = B(3);
    int_0 = B(4);
    b_0   = B(5);
    
    % Residuals at the initial guess
    Si    = sum(((fi - psf(xs_0, ys_0, sg_0, int_0, b_0, xi, yi)).^2));

    % Filling in initial guesses and datapoints 
    dfdxs_0  = dfdxs(xs_0, ys_0, sg_0, int_0, b_0, xi, yi);
    dfdys_0  = dfdys(xs_0, ys_0, sg_0, int_0, b_0, xi, yi);
    dfdsg_0  = dfdsg(xs_0, ys_0, sg_0, int_0, b_0, xi, yi);
    dfdin_0  = dfdin(xs_0, ys_0, sg_0, int_0, b_0, xi, yi);
    dfdb_0   = ones(xsize*ysize, 1);

    % Jacobian Matrix
    J = [dfdxs_0, dfdys_0, dfdsg_0, dfdin_0, dfdb_0];
    JT = transpose(J);
    
    % Residuals
    res = (JT*((fi - psf(xs_0, ys_0, sg_0, int_0, b_0, xi, yi))));
    jtj = JT*J;
    
    % Damping matrix according to 'Improvements to the Levenberg-Marquardt algorithm for nonlinear
    %least-squares minimization'
    dtd = jtj;
    
    % Calculate the step size
    step_0 = (jtj + L_0*dtd)\res;
    step_v = (jtj + (L_0/v)*dtd)\res;
    
    % Update the function parameters
    B_0  = B + step_0;
    B_v  = B + step_v;
    
    S_0 = sum(((fi - psf(B_0(1), B_0(2), B_0(3), B_0(4), B_0(5), xi, yi)).^2));
    S_v = sum(((fi - psf(B_v(1), B_v(2), B_v(3), B_v(4), B_v(5), xi, yi)).^2));
    
    % If both choices for damping parameter are worse, multiply by v and
    % guess again
    if S_0>Si && S_v>Si
        L_0 = v*L_0;
        B   = B + ((jtj + L_0*dtd)\res);
    else if S_0>S_v
            B = B_0;
        else
            B = B_v;
            L_0 = L_0/v;
        end
    end
            
end
    
    xs  = B(1) + roi(1) - 1.5; % +0.5 is necessary because pixel count starts at one in the local frame
    ys  = B(2) + roi(2) - 1.5; % and because the center of a pixel is where the count starts
    sg  = B(3);
    b   = B(5);
    
    localizations(i, :) = [xs, ys, sg, b, iter];
    
end


end

function [comp, missed, Nmol] = groundtruth(localizations, groundtruth_loc, iter)

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

% Number of detections
Nmol = N_fd(1);

% Missed detections
missed = N_pr(1)- N_fd(1);

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
comp = norm(sum(sqrt(rsq))/L_rsq(1));


end

function imageData = OpenIm(location, iter)

% Format the iteration number to match the name of the data
num  = '00000';
num(5-strlength(string(iter)) + 1:5) = string(iter);
dir = append(location, num ,'.tif');

% Opening the raw .tif file
t = Tiff(dir, 'r');

% Reading the image
imageData = double(read(t));

close(t);

end


