function [] = scalebar(imageName,pixelSize,barScale,units)

%sets a scale bar in the image: imageName = tot_im
%pixelSize = rec_px;
%barScale = length of bar e.g. 1000 for 1000 nm
%units: 'nm' or 'um'

image = imageName;

%plot a scale bar
scaleBarWidth = floor( 1/(pixelSize) * barScale);
scaleBarHeight = floor( 1/(pixelSize) * barScale)/30;

xPos = size(image,2)*0.90 - scaleBarWidth;
yPos = size(image,1)*0.90 - scaleBarHeight;
textCenterX = xPos + floor(scaleBarWidth/2);
textCenterY = yPos + scaleBarHeight*5;
rectPosition = [xPos, yPos, scaleBarWidth, scaleBarHeight];
hRect = rectangle('Position', rectPosition);

%label the scale bar
str = sprintf(['%4d ' units], barScale);
hText = text(textCenterX,textCenterY,str);
set(hText,'HorizontalAlignment','center');

%set white tone for scale bar
set(hRect,'EdgeColor','w');
set(hRect,'FaceColor','w');
set(hText,'Color','w');

end