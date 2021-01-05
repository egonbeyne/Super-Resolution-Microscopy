
clear all;
tic

% Constants
im_px       = 100;      %[nm] Pixel size
rec_px      = 10;       %[nm] accuracy of reconstructed image
xsize       = 4;        %[px] size of fitting region in x
ysize       = 4;        %[px] size of fitting region in y
dataLoc     = 'C:\Users\Egon Beyne\Desktop\Super-Resolution Microscopy\Data\sequence_3\';
GtLoc       = 'C:\Users\Egon Beyne\Desktop\Super-Resolution Microscopy\Data\Ground-truth\';
Nfiles      = 19996;   %Number of datafiles
resolution  = size(OpenIm('C:\Users\Egon Beyne\Desktop\Super-Resolution Microscopy\Data\sequence_3\', 1));
nxpixels    = resolution(2);   % number of pixels in x direction
nypixels    = resolution(1);   % number of pixels in y direction

% Making a point spread function
psf = @(xs, ys, sg, int, b, x, y)((int/(2*pi*(sg^2)))*exp(-((x-xs).^2 + (y-ys).^2)/(2*(sg^2)))+b); 


% Convert derivatives to normal functions for speed
dfdxs = @(xs, ys, sg, int, b, x, y)(int*exp(-((x - xs).^2 + (y - ys).^2)/(2*sg^2)).*(2*x - 2*xs))/(4*sg^4*pi);
dfdys = @(xs, ys, sg, int, b, x, y)(int*exp(-((x - xs).^2 + (y - ys).^2)/(2*sg^2)).*(2*y - 2*ys))/(4*sg^4*pi);
dfdsg = @(xs, ys, sg, int, b, x, y)((int*exp(-((x - xs).^2 + (y - ys).^2)/(2*sg^2)).*((x - xs).^2 + (y - ys).^2))/(2*sg^5*pi) - (int*exp(-((x - xs).^2 + (y - ys).^2)/(2*sg^2))./(sg^3*pi)));
dfdin = @(xs, ys, sg, int, b, x, y)exp(-((x - xs).^2 + (y - ys).^2)/(2*sg^2))./(2*sg^2*pi);
dfdb  = @(xs, ys, sg, int, b, x, y)(1);

% Empty array for reconstruction
tot_im  = zeros((nypixels+9)*im_px/rec_px, (nxpixels+9)*im_px/rec_px);

% Array to store accuracy information
acc     = zeros(Nfiles, 3);

for iter = 1:Nfiles

    % Open the image
    imageData = OpenIm('C:\Users\Egon Beyne\Desktop\Super-Resolution Microscopy\Data\sequence_3\', iter);

    % Segmenting the data
    centroids = segment_frame(imageData);

    % Skip current loop iteration if no molecules are present (lots of
    % detections) (The following lines can probably be removed)
    if sum(isnan(centroids), 'all')>0
        continue
    end

    % Localization
    localizations = Fit_Gaussian(centroids, imageData, xsize, ysize, iter, dfdxs, dfdys, dfdsg, dfdin, psf, resolution);
    
    % Check if any warnings are given
    [warnmsg, msgid] = lastwarn;
    
    % If the warning is about an ill conditioned matrix, skip the frame
    % (Gives NaN as fitting parameters, which makes the reconstruction a black screen)
    if strcmp(msgid,'MATLAB:illConditionedMatrix')
        
        % Reset the warning
        lastwarn(['','']);
        continue
    end
    
    % Stop current iteration if no localizations are present
    if size(localizations, 1)==0
        continue
    end
    
    xc = localizations(:, 1)*im_px;
    yc = localizations(:, 2)*im_px;

    % Cramer-Rao lower bound calculation
    N     = mean(localizations(:, 6)):                                  % [-] Number of signal photons 
    sigg  = mean(localizations(3))*im_px*10^-9;                         % nm] width of blob converted to nm
    sige2 = (sigg^2) + (((im_px*10^-9)^2)/12);                          
    tau   = 2*pi*(sige2)*mean(localizations(4))/(N*((im_px*10^-9)^2));  % [-] Dimensionless background parameter
    
    % Cramer rao lower bound
    dx   = sqrt(sige2*(1 + (4*tau) + sqrt(2*tau/(1 + (4*tau))))/N)/1e-9;

    % Store some statistics
    acc(iter, :) = [iter, dx, mean(localizations(4))];
    
    % Add the localization data to the reconstructed image
    tot_im = reconstruct(tot_im, localizations, 100, rec_px, nxpixels, nypixels);

    % Print progress
    progress = string(iter/Nfiles*100) + '% done';
    disp(progress)
end

% Test statistics
figure(1)
histogram(acc(:, 2))

%%%%%%%%
axis equal
figure(2)
imshow(tot_im)

% Execution time
toc

function [localizations] = Fit_Gaussian(centroids, imageData, xsize, ysize, iter,...
    dfdxs, dfdys, dfdsg, dfdin, psf, resolution)

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

localizations = zeros(L(1), 7);

for i = 1:L(1)
    
    % Extract coordinates of the origin of one region of interest
    roi = origin(i, :);
    
    
    % Prevent the search domain being outside the image
    roi(1) = min(max(roi(1), 1), resolution(1)- ysize+1);
    roi(2) = min(max(roi(2), 1), resolution(2) - xsize +1);
       
    % Select pixel data in the region of interest
    localData = imageData(roi(2):(roi(2)+xsize-1), roi(1):(roi(1)+ysize-1));

    % rearranging ROI
    Ii  = transpose(localData);
    I = Ii(:);

    % Calculating the weights (assuming poissoning noise, based on
    % (Jiaqing,2016))
    w = max(I.^-1, .005);

    fi = Ii(:);

% Initial guesses for x and y
xin = centroids(i, 2) - roi(2);
yin = centroids(i, 1) - roi(1);

% Initial guess for intensity (based on the sum of the intensity in the
% local frame
Iin = sum(localData, 'all');

% Initial guess for b, based on average intensity of the entire image,
% scaled with a factor 2 to account for fluorescent molecules
bin = sum(imageData, 'all')/(resolution(1)*resolution(2)*2);

% Initial guess for function parameters
B     = [xin; yin; 1.5; Iin; bin];

% x and y step size to start up iterations (exaggerated)
step = [4, 4];

% Maximum number of iterations
maxIt = 30;

it = 0;

% Damping for levenberg marquardt (Has to be tuned)
L_0 = 1000;
v   = 1.5;

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
    dtd = eye(5)*diag(jtj);
    
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
    Int = B(4);
    
    
    localizations(i, :) = [xs, ys, sg, b, iter, Int, Iin];
    
    
end


% Filter out localizations using a molecule that is too dim
localizations(localizations(:, 6)<5e3,:) = [];

% Filter out double localizations
localizations(localizations(:, 3)<0.5,:) = [];

localizations(:, 6) = [];

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


