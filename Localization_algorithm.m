
iter = 2;
location = 'C:\Users\Egon Beyne\Desktop\Super-Resolution Microscopy\Data\sequence\';
% Format the iteration number to match the name of the data
num  = '00000';
num(5-strlength(string(iter)) + 1:5) = string(iter);
dir = append(location, num ,'.tif');

%'C:\Users\Egon Beyne\Desktop\Super-Resolution Microscopy\Data\sequence\00002.tif';
% Opening the raw .tif file
t = Tiff(dir, 'r');

% Reading the image
imageData = double(read(t));

sel = imageData(129:132, 31:34);

close(t);


% Empty array to append localizations to
localizations = [];
xsize = 4;
ysize = 4;

X = linspace(1, xsize, xsize);
Y = linspace(1, ysize, ysize);

% Grid coordinates for every point in the ROI
[Xi, Yi] = meshgrid(X, Y);

% Making a 2 dimensional array for all points in the grid
c = cat(2,Xi',Yi');
d = reshape(c,[],2);

xi = d(:, 1);
yi = d(:, 2);

% Function values at datapoints
% rearranging ROI
Ii  = transpose(sel);
fi = Ii(:);


% Initial guess for function parameters
B     = [2;2;1;5000; 100];



% Making a symbolic point spread function
syms psf(xs, ys, sg, int, b, x, y)
psf(xs, ys, sg, int, b, x, y) = ((int/(2*pi*(sg^2)))*exp(-((x-xs).^2 + (y-ys).^2)/(2*(sg^2)))+b); 

% Derivative of the psf wrt to the function parameters
dfdxs = diff(psf, xs);
dfdys = diff(psf, ys);
dfdsg = diff(psf, sg);
dfdin = diff(psf, int);
dfdb  = diff(psf, b);

step = [4, 4];
maxIt = 10;

it = 0;

while (step(1)^2 + step(2)^2)>0.001 && it<maxIt

    % Counting iterations
    it = it + 1;
    
    % Unpack values from b
    xs_0  = B(1);
    ys_0  = B(2);
    sg_0  = B(3);
    int_0 = B(4);
    b_0   = B(5);


    % Filling in initial guesses and datapoints 
    dfdxs_0  = double(dfdxs(xs_0, ys_0, sg_0, int_0, b_0, xi, yi));
    dfdys_0  = double(dfdys(xs_0, ys_0, sg_0, int_0, b_0, xi, yi));
    dfdsg_0  = double(dfdsg(xs_0, ys_0, sg_0, int_0, b_0, xi, yi));
    dfdin_0  = double(dfdin(xs_0, ys_0, sg_0, int_0, b_0, xi, yi));
    dfdb_0   = double(dfdb(xs_0, ys_0, sg_0, int_0, b_0, xi, yi));

    % Jacobian Matrix
    J = [dfdxs_0, dfdys_0, dfdsg_0, dfdin_0, dfdb_0];

    JT = transpose(J);
    
    sum(fi - double(psf(xs_0, ys_0, sg_0, int_0, b_0, xi, yi)))

    % Damping for levenberg marquardt (Has to be tuned)
    damp = 1;
    
    % Calculate the step size
    step = inv(JT*J + damp*eye(5))*(JT*(fi - double(psf(xs_0, ys_0, sg_0, int_0, b_0, xi, yi))));


    % Update the function parameters
    B = B + step;
    
end
imshow(uint16(sel)*100,'InitialMagnification', 'fit');
axis on
hold on
plot(B(1), B(2), 'r+')
