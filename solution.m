% CT head test

%clc; clear all ; close all;
%% Load images
display('loading images');

close all;

files = dir('./tiff/*.tif');

N = length(files);

I={};
for i=1:N
    im =imread(strcat('./tiff/',files(i).name));
    I{i} = im;
end

%% Normalize images
minVec=zeros(1,N);
maxVec=zeros(1,N);
for i=1:N
    im = I{i};
    minVec(i) = min(im(:));
    maxVec(i) = max(im(:));
end

minVal = min(minVec);
maxVal = max(maxVec);

for i=1:N
    im = I{i}-minVal;
    im = double(im)/(maxVal-minVal);
    I{i} = im;    
end


for i=1:N
    im = I{i};
    minVec(i) = min(im(:));
    maxVec(i) = max(im(:));
end

%% Generate scale-normalized data volume

sliceSep = 0.625;
totalDepth = sliceSep*(N-1);
pixelSize = 1/2.048;
Nnew = round(totalDepth/pixelSize); 

x=[0:N-1]*sliceSep+sliceSep/2;
xq=[0:Nnew-1]*pixelSize+pixelSize/2;

[m,n] = size(I{1});

Ivol = zeros(m,n,Nnew);

for r=1:m
    for c=1:n
        int = zeros(N,1);
        for i=1:N
            int(i) = I{i}(r,c);
        end

        newInt = interp1(x,int,xq);

        for i=1:Nnew
            Ivol(r,c,i) = newInt(i);
        end
    end
end

%% Try crazy plotting

% figure;
% for r=1:m
%     for c=1:n
%         for d=1:Nnew  %d: depth
%             pint = Ivol(r,c,d);
%             if pint > 0.7
%                 plot3(r,c,d,'b.');
%                 hold on;
%             end
%         end
%     end
% end

%% Get top DRR
[m,n]=size(I{1});

DRRtop = zeros(m,n);
for i=1:N
    im = I{i};
    DRRtop = DRRtop + im;
end

DRRtop = DRRtop/N;

%% Get lateral DRR
sliceSep = 0.625;
totalHeight = sliceSep*(N-1);
pixelSize = 1/2.048;
DRRrows = round(totalHeight/pixelSize); % we'll get the same resolution vertically and horizontally, and it will be defined by the pixelSize

x=[0:N-1]*sliceSep+sliceSep/2;
xq=[0:DRRrows-1]*pixelSize+pixelSize/2;

DRRlateralUS = zeros(N,m);
for i=1:N
    im = I{N-(i-1)};
    row = sum(im,2)';
    DRRlateralUS(i,:) = row;
end

DRRlateralUS = DRRlateralUS/n;

DRRlateral = zeros(DRRrows,m);
for i=1:m
    col = DRRlateralUS(:,i);
    newCol = interp1(x,col,xq);
    DRRlateral(:,i) = newCol;
end

%% Get frontal DRR

DRRfrontalUS = zeros(N,m);
for i=1:N
    im = I{N-(i-1)};
    row = sum(im,1);
    DRRfrontalUS(i,:) = row;
end

DRRfrontalUS = DRRfrontalUS/n;

DRRfrontal = zeros(DRRrows,m);
for i=1:m
    col = DRRfrontalUS(:,i);
    newCol = interp1(x,col,xq);
    DRRfrontal(:,i) = newCol;
end

%% Show DDRs
figure;
subplot(131);
imshow(DRRtop,[]);
title('top DRR');
subplot(132);
imshow(DRRfrontal,[]);
title('frontal DRR');
subplot(133);
imshow(DRRlateral,[]);
title('lateral DRR');

%% Prompt for lateral points

[xlateral, ylateral] = extractContour(DRRlateral, 10);

%% Fit ellipse
[ag, bg, x0, y0, phi] = ellipsefit(xlateral,ylateral);
hold on
axis equal
plotellipse([x0;y0], ag, bg, phi, 'k')
hold off

%% Prompt for frontal points
[xfrontal, yfrontal] = extractContour(DRRfrontal, 10);

%% Fit circle
[xc,yc,R] = circfit(xfrontal,yfrontal)
hold on
axis equal
plotcircle(xc,yc,R);
hold off





