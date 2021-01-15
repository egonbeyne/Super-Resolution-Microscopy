clear all;
tic

% Constants
im_px       = 100;      %[nm] Pixel size of frames
rec_px      = 5;       %[nm] pixel size of reconstructed image
xsize       = 5;        %[px] size of fitting region in x
ysize       = 5;        %[px] size of fitting region in y
%dataLoc     = 'C:\Users\Egon Beyne\Desktop\Super-Resolution Microscopy\Data\sequence_3\';
dataLoc     = 'C:\Users\kaan_\EE\minor\Final project\matlab code\data\ER2.N3.HD\';
Nfiles      = length( dir(dataLoc)) - 2;
%GtLoc       = 'C:\Users\Egon Beyne\Downloads\poitions.csv';%'C:\Users\Egon Beyne\Desktop\Super-Resolution Microscopy\Data\ground-truth\';
resolution  = size(OpenIm(dataLoc, 1));
nxpixels    = resolution(2);   % number of pixels in x direction
nypixels    = resolution(1);   % number of pixels in y direction
psf_pixels  = 9;               %size of psf_block
last_prog   = -1;              % Progress indicator

%init reconstruction, makes psf block                                        
% psf_block = init_reconstruct(rec_px);

warning('off', 'all');

% Making a point spread function
psf = @(xs, ys, sg, int, b, x, y)((int/(2*pi*(sg^2)))*exp(-((x-xs).^2 + (y-ys).^2)/(2*(sg^2)))+b); 


% Convert derivatives to normal functions for speed
dfdxs = @(xs, ys, sg, int, b, x, y)(int*exp(-((x - xs).^2 + (y - ys).^2)/(2*sg^2)).*(2*x - 2*xs))/(4*sg^4*pi);
dfdys = @(xs, ys, sg, int, b, x, y)(int*exp(-((x - xs).^2 + (y - ys).^2)/(2*sg^2)).*(2*y - 2*ys))/(4*sg^4*pi);
dfdsg = @(xs, ys, sg, int, b, x, y)((int*exp(-((x - xs).^2 + (y - ys).^2)/(2*sg^2)).*((x - xs).^2 + (y - ys).^2))/(2*sg^5*pi) - (int*exp(-((x - xs).^2 + (y - ys).^2)/(2*sg^2))./(sg^3*pi)));
dfdin = @(xs, ys, sg, int, b, x, y)exp(-((x - xs).^2 + (y - ys).^2)/(2*sg^2))./(2*sg^2*pi);
dfdb  = @(xs, ys, sg, int, b, x, y)(1);

% Empty array for reconstruction
tot_im  = zeros((nypixels+psf_pixels+2)*im_px/rec_px, (nxpixels+psf_pixels+2)*im_px/rec_px);

% Array to store accuracy information
acc     = zeros(Nfiles, 5);

% Array to store localizations
loc_mol = zeros(20000, 3);

for iter = 1:Nfiles

    % Open the image
    imageData = OpenIm(dataLoc, iter);

    % Segmenting the data
    centroids = segment_frame(imageData, im_px);

    % Skip current loop iteration if no molecules are present
    if sum(isnan(centroids), 'all')>0 || size(centroids, 1) == 0
        %imshow(10*uint16(imageData))
        %error('stop')
        continue
    end

    % Localization
    localizations = Fit_Gaussian(centroids, imageData, xsize, ysize, iter, dfdxs, dfdys, dfdsg, dfdin, psf, resolution);
    
    % Stop current iteration if no localizations are present
    if size(localizations, 1)==0
        continue
    end
    
    % Total number of localizations
    Nloc = nnz(loc_mol(:, 3));
    
    
    % Remove zeros from localization
    localizations(~any(localizations,2), :) = [];
   
    % Number of localizations in the current frame
    Nfr  = size(localizations, 1);
        
    % Store localizations, also add the frame number
    loc_mol((Nloc+1):(Nloc+Nfr), 1:3) = [localizations(:, 1:2)*im_px, iter*ones(Nfr, 1)];
    
    % Cramer-Rao lower bound calculation
    N     = mean(localizations(:, 6)) - (mean(localizations(:, 4))*xsize*ysize);      % [-] Number of signal photons
    sigg  = mean(localizations(:, 3))*im_px*10^-9;                                    % nm] width of blob converted to nm
    sige2 = (sigg^2) + (((im_px*10^-9)^2)/12);                          
    tau   = 2*pi*(sige2)*mean(localizations(:, 4))/(N*((im_px*10^-9)^2));             % [-] Dimensionless background parameter
    
    % Cramer rao lower bound
    dx   = sqrt(sige2*(1 + (4*tau) + sqrt(2*tau/(1 + (4*tau))))/N)/1e-9;
    
    %[comp1, missed, Nmol] = groundtruth(localizations, GtLoc, iter);
    comp1 = 1;
    
    % Store some statistics
    acc(iter, :) = [iter, dx, mean(localizations(:, 3)), N, comp1];
    
    % Add the localization data to the reconstructed image
    tot_im = reconstruct(tot_im, localizations, im_px, rec_px, nxpixels,dx);

    % Print progress
    prog = int16(iter/Nfiles*100);
    if prog ~= last_prog
        disp(string(prog) + '% done ')
    end
    last_prog = prog;
end

%%%%%%%%
figure(1)
axis equal
imshow(tot_im*100,hot(60))
hold on
scalebar(tot_im,rec_px,1000,'nm')

%%%%%%%%
% remove zero elements from the molecule locations
loc_mol(loc_mol(:, 1) == 0, :) = [];


% Progress report
disp("Checking for double localizations... ")

% Checking for molecules that are on in consecutive frames
% guess for number of consecutive frames (best to take high enough)
Ndoubleframes = 10;

% Array to store how many frames a localization is on
N_on          = zeros(Nloc, 1);

% Threshold to consider 2 localization from the same molecule (in nm)
thr           = 5;

% Loop through all frames
for i = 1:Nfiles
    % Select localization of 1 frame, to check with 10 consecutive frames
    local = loc_mol(loc_mol(:, 3) == i, :);
    %i
    % Loop through all localizations in the selected frame
    for j = 1:size(local, 1)
        %j
        % array to attach doubles to
        doubles = local(j, 1:2);  
        
        % Loop through frames to be checked (10 consecutive frames)
        for k = (1+i):(Ndoubleframes+i)
            
            % Select frame to be checked
            fr = loc_mol(loc_mol(:, 3) == k, :);
            
            % Calculate distance to localization
            [diff, index] = min(sqrt(sum((fr(:, 1:2) - local(j, 1:2)).^2   , 2)));
            
            if diff<thr
                
                % Store doubles together
                doubles = [doubles; fr(index, 1:2)];
                
                % Remove localization from other frames ( Will be replaced
                % by average later)
                fr(index, :) = [0, 0, k];
         
            end
            % Replace checked frame by frame with removed doubles
            loc_mol(loc_mol(:, 3)==k, :) = fr;
        end
        
        % Calculate the average from all doubles
        avg_loc = mean(doubles, 1);
        
        % on-time of the fluorophore
        t_on    = size(doubles, 1);
        N_on(i) = t_on;
        
        % replace with average
        local(j, 1:2) = avg_loc;
    end
end

% % Compute the accuracy
% comp = mean(acc(:, 5));            % When using eye dataset
% comp = groundtruth_combined(loc_mol, GtLoc);    % For other dataset
% 
% % Cramer-Rao Lower bound, without zeros
 CRLB = nonzeros(acc(:, 2));
% 
% % number of signal photons 
 N_ph = nonzeros(acc(:, 4));
% 
% figure(2)
% subplot(221)
% histogram(nonzeros(N_on))
% ylabel("Frequency")
% xlabel("# of frames on")
% title("Number of frames fluorophore is on")
% subplot(222)
% histogram(CRLB)
% ylabel("Frequency")
% xlabel("Cramer-Rao bound [nm]")
% title("Cramer-Rao lower bound")
% xlim([0, 15])
% subplot(223)
% histogram((CRLB/comp).^-1)
% title("Cramer-Rao Bound scaled with accuracy")
% subplot(224)
% histogram(N_ph, 'NumBins', 10)
% ylabel("Frequency")
% xlabel("Number of signal photons")
% title("Signal photons per fluorophore")
  
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
origin = round(centroids) - [(xsize+1)/2, (ysize+1)/2];

L = size(centroids);

localizations = zeros(L(1), 7);

for i = 1:L(1)
    
    % Extract coordinates of the origin of one region of interest
    roi = origin(i, :);
    
    % Prevent the search domain being outside the image
    roi(1) = min(max(roi(1), 1), resolution(2)- xsize);
    roi(2) = min(max(roi(2), 1), resolution(1) - ysize);
    
    % Select pixel data in the region of interest
    localData = imageData((roi(2)+1):(roi(2)+xsize), (roi(1)+1):(roi(1)+ysize));
    
    
    % rearranging ROI
    Ii  = transpose(localData);
    I = Ii(:);

    % Calculating the weights (assuming poissoning noise, based on
    % (Jiaqing,2016))
    w = max(I.^-1, .0001);

    fi = Ii(:);

% Initial guesses for x and y
xin = centroids(i, 1) - roi(1);
yin = centroids(i, 2) - roi(2);

% Initial guess for intensity (based on the sum of the intensity in the
% local frame
Iin = sum(localData, 'all');

% Initial guess for b, based on average intensity of the entire image,
% scaled with a factor 2 to account for fluorescent molecules
bin = sum(imageData, 'all')/(resolution(1)*resolution(2)*1.5);

% Initial guess for function parameters
B     = [xin; yin; 1; Iin; bin];

% x and y step size to start up iterations (exaggerated)
step = [4, 4];

% Maximum number of iterations
maxIt = 500;

it = 0;

% Damping for levenberg marquardt (Has to be tuned)
L_o     = 2;
v       = 1.5;
h       = 0.1;
up      = 20;
down    = 5;
L_0 = 4;
alpha = 0.1;
while (step(1)^2 + step(2)^2)>1e-8 && it<maxIt

    % Counting iterations
    it = it + 1;
    
    % Unpack values from b
    xs_0  = B(1);
    ys_0  = B(2);
    sg_0  = B(3);
    int_0 = B(4);
    b_0   = B(5);
    
    % Function value at current iteration
    f   = psf(xs_0, ys_0, sg_0, int_0, b_0, xi, yi);
    
    % Filling in initial guesses and datapoints 
    dfdxs_0  = dfdxs(xs_0, ys_0, sg_0, int_0, b_0, xi, yi);
    dfdys_0  = dfdys(xs_0, ys_0, sg_0, int_0, b_0, xi, yi);
    dfdsg_0  = dfdsg(xs_0, ys_0, sg_0, int_0, b_0, xi, yi);
    dfdin_0  = dfdin(xs_0, ys_0, sg_0, int_0, b_0, xi, yi);
    dfdb_0   = ones(xsize*ysize, 1);

    % Jacobian Matrix
    J = [dfdxs_0, dfdys_0, dfdsg_0, dfdin_0, dfdb_0];
    JT = transpose(J);
    
    % Weights
    W = eye(length(xi)).*w;
    
    % Residuals
    res = JT*W*(fi - f);
    jtj = JT*W*J;
    
    % Damping matrix according to 'Improvements to the Levenberg-Marquardt algorithm for nonlinear
    %least-squares minimization'
    dtd = eye(5).*diag(jtj);

    % Calculate the step size
    step = linsolve((jtj + L_0*dtd),res);
    
    %%%% Geodesic %%%%
    Bs   = B +  h*step;
    fs   = psf(Bs(1), Bs(2), Bs(3), Bs(4), Bs(5), xi, yi);
    
    fvv  = 2*((fs-f)/h - J*step)/h;
    
    % Solve for acceleration
    ag   = linsolve(J, -fvv);
    
    % Add accelerator to step size if convergence is nt slowed down
    if 2*norm(ag)/norm(step)<= alpha
        step = step + 0.5*ag;
    end
    %%%%%%%%%%%%%%%%%%
    
    % Update function parameters
    Bn  = B + step;
    
    % Function value at new point
    fh = psf(Bn(1), Bn(2), Bn(3), Bn(4), Bn(5), xi, yi);
    
    chi     = fi.'*W*fi  -  2*fi.'*W*f  +  f.'*W*f;
    chid    = fi.'*W*fi  -  2*fi.'*W*fh  +  fh.'*W*fh;
    rho     = (chi - chid)/(step.'*((L_0*dtd*step) + (JT*W*(fi - f))));
    
    % Use new parameters if the least squares improved
    if rho>=0.5
        B = Bn;
        L_0 = L_0*max(0.333, 1 - (2*rho - 1)^3);%max(L_0/down, 1e-7);
        v   = 2;
        %disp("better")
    else
        
        L_0 = v*L_0;%min(L_0*up, 1e7);
        v   = 2*v;
        %disp("worse")
    end
    
    
    % Check if any warnings are given
    [~, msgid] = lastwarn;
    
    % If the fitting produces an ill conditioned matrix, or the fluorophore
    % is not bright enough, stop iteration
    if strcmp(msgid,'MATLAB:illConditionedMatrix') || Iin<numel(localData)*8
        
        if strcmp(msgid,'MATLAB:illConditionedMatrix')
            % Reset the warning
            lastwarn(['','']);
            disp('ill conditioned matrix')
        else
            disp('too dim')
        end
        B = [NaN;NaN;NaN;NaN;NaN];
        break
    end
    
    
    
       
end
    if it==maxIt
        disp('itmax reached')
    end
    
    xs  = B(1) + roi(1) - 0.5; % -0.5 is necessary because pixel count starts at one in the local frame
    ys  = B(2) + roi(2) - 0.5; % and because the center of a pixel is where the count starts
    sg  = B(3);
    b   = B(5);
    Int = B(4);
    %{
    if sqrt((B(1) - xin)^2 + (B(2) - yin)^2)>1
        figure(1)
        imshow(uint16(localData)*30, 'InitialMagnification', 800);
        hold on
        axis on
        plot(xin, yin, 'r+')
        plot(B(1), B(2), 'bx')
        Iin
        numel(localData)*150
        figure(2)
        imshow(uint16(imageData)*30);
        hold on
        axis on
        plot(centroids(:,1), centroids(:,2), 'r+')
        pause(5)
    end
    %}
    localizations(i, :) = [xs, ys, sg, b, iter, Int, Iin];
    
    
end

% Remove diverged fittings
localizations((localizations(:, 1) <= -1)|(localizations(:, 2) <=-1), :) = [];
localizations((localizations(:, 1) >= resolution(2)+1)|(localizations(:, 2) >=resolution(1)+1), :) = [];

% Remove invalid fitting (b<0, I<0, NaN)
localizations(localizations(:, 4)<0, :) = [];
localizations(localizations(:, 6)<0, :) = [];
localizations(sum(isnan(localizations), 2)>0, :) = [];


end


function comp = groundtruth_combined(loc_mol, groundtruth_loc)

% Open and read the csv file containing ground truth data
gt = readmatrix(groundtruth_loc);

% Array to store localization errors
C      = zeros(size(gt, 1), 1);

% Iterate through ground truth and find the closest localization to each
% ground truth molecule.

for k=1:size(gt, 1)
  
  % value and index of minimum error
  [val,idx] = min(sqrt(sum((gt(:, 1:2)-loc_mol(k, 1:2)).^2, 2)));
  
  % Store the error
  C(k)=val;
  
  % Remove entry from ground truth
  gt(idx, :)=[];
end

% Trim zero values
C(~any(C), :) = [];

% Compute the average accuracy, excluding outliers
comp = trimmean(C, 99.5);
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