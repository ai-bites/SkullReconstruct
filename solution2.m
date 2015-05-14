%% CT head test

clc; clear all ; close all;

%% Load images
display('loading images');

close all;

files = dir('./tiff/*.tif');

N = length(files);

im =imread(strcat('./tiff/',files(1).name));

[dim1,dim2] = size(im);

I=zeros(dim1,dim2,N);
H = padarray(2,[2 2]) - fspecial('gaussian' ,[5 5],2);
for i=1:N
    im = imread(strcat('./tiff/',files(i).name));
    I(:,:,i) = im;
end

%% Normalize images

disp('normalizing image intensitites');
minVal = min(I(:));
maxVal = max(I(:));

I = (I - minVal)/(maxVal-minVal);    

%% Generate scale-normalized data volume

disp('scale normalized data volume generation...');

sliceSep = 0.625;
totalDepth = sliceSep*(N-1);
pixelSize = 1/2.048;
dim3 = round(totalDepth/pixelSize); 

x=[0:N-1]*sliceSep+sliceSep/2;
xq=[0:dim3-1]*pixelSize+pixelSize/2;

Ivol = zeros(dim1,dim2,dim3);

for r=1:dim1
    for c=1:dim2
        Ivol(r,c,:) = interp1(x,squeeze(I(r,c,:)),xq);
    end
end

Ivol = Ivol(:,:,end:-1:2);
dim3 = dim3 - 1;
Ivol2 = zeros(dim1,dim2, dim3+50);
Ivol2(:,:,51:end) = Ivol;
Ivol = Ivol2;
dim3 = size(Ivol,3);



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

disp('generating top DRR');
DRRtop = zeros(dim1,dim2);
for i=2:dim3
    im = squeeze(Ivol(:,:,i));
    DRRtop = DRRtop + im;
end

DRRtop = DRRtop/dim3;

%% Get lateral DRR

disp('generating lateral DRR');
DRRlateral = zeros(dim3,dim1);
for i=1:dim2
    im = squeeze(Ivol(:,i,:));
    DRRlateral = DRRlateral + im';
end

DRRlateral = DRRlateral/dim2;

%% Get frontal DRR

disp('generating frontal DRR');
DRRfrontal = zeros(dim3,dim2);
for i=1:dim1
    im = squeeze(Ivol(i,:,:));
    DRRfrontal = DRRfrontal + im';
end

DRRfrontal = DRRfrontal/dim1;


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
plotellipse([x0;y0], ag, bg, phi, 'k');
hold off

%% Prompt for frontal points
[xfrontal, yfrontal] = extractContour(DRRfrontal, 10);

%% Fit circle
[xc,yc,R] = circfit(xfrontal,yfrontal)
hold on
axis equal
plotcircle(xc,yc,R);
hold off

%% Generate ellipsoid
meshres = 50;

ellCenterX = x0;
ellCenterY = xc;
ellCenterZ = (y0 + yc) /2;
%ellCenterZ = yc;
ellRadX = ag;
ellRadY = R;
ellRadZ = ((y0+bg+yc+R)/2) - ellCenterZ;
%ellRadZ = R;
[XmeshIn,YmeshIn,ZmeshIn] = ellipsoid(ellCenterX,ellCenterY,ellCenterZ, ...
                                        ellRadX,ellRadY,ellRadZ,meshres);
XmeshIn = XmeshIn(1:ceil(meshres/2)+1,:);
YmeshIn = YmeshIn(1:ceil(meshres/2)+1,:);
ZmeshIn = ZmeshIn(1:ceil(meshres/2)+1,:);
XmeshOut = XmeshIn;
YmeshOut = YmeshIn;
ZmeshOut = ZmeshIn;


%% Validate ellipsoid

fh = figure

surf(XmeshIn,YmeshIn,ZmeshIn);
axis equal
axis([1 dim2 1 dim1 1 dim3]);
hold on
xlabel('x');
ylabel('y');
zlabel('z');

view([0 0]);
h = gca;
F = getframe(h);
projLat = F.cdata;


scale = size(DRRlateral, 1)/size(projLat,1);
projLat=imresize(projLat, scale);

projLat=projLat(:,1:dim1,:);
projLat = projLat(end:-1:1,:,:);
projLat = mat2gray(rgb2gray(projLat));
checkLat = (1-projLat)*0.05+mat2gray(DRRlateral)*2;
figure; imshow(checkLat,[])

%surfnorm(Xmesh,Ymesh,Zmesh);

%% Debugging stuff

%figure;
%topX = XmeshIn(51,1:26);
%topY = YmeshIn(51,1:26);
%topZ = ZmeshIn(51,1:26);

%plot3(topX,topY,topZ);
%axis equal

%slice = squeeze(Ivol(round(XmeshIn(parallel,meridian)),:,:));
%slice = slice';
%figure;
%imshow(slice,[]);


%% Compute intensities along normals

[Nx,Ny,Nz] = surfnorm(XmeshIn,YmeshIn,ZmeshIn);

[Xgrid,Ygrid,Zgrid] = meshgrid(1:dim1,1:dim2,1:dim3);

[parallelN, meridianN] = size(XmeshIn);

%parallel = 10;

threshPercent = 0.8;
%figure
for iters=1:20
for parallel=1:parallelN
for meridian=1:meridianN
    
    p = [XmeshIn(parallel,meridian),YmeshIn(parallel,meridian), ...
         ZmeshIn(parallel,meridian)];
    n = [Nx(parallel,meridian),Ny(parallel,meridian),Nz(parallel,meridian)];

    idx = 1;
    pline = [];
    Icurr = [];
    soiLength = 30;
    for i=-soiLength:soiLength
        pline(idx,:) = p+i*n;
        Icurr(idx) = Ivol(round(pline(idx,1)),round(pline(idx,2)),round(pline(idx,3)));
        idx = idx+1;
    end
    % edge detection
    [intMax intMaxIdx] = max(Icurr);
    thresh = intMax*threshPercent;
    
    k = 25;
    cropStart = max(intMaxIdx-k,1);
    cropEnd = min(intMaxIdx+k,length(Icurr));
    IcurrCropped = Icurr(cropStart:cropEnd);
    
    
    [vals idxs] = find(IcurrCropped > thresh);
    minIdx = min(idxs)+cropStart-1;
    maxIdx = max(idxs)+cropStart-1;

    pLeft = pline(minIdx,:);
    pRight = pline(maxIdx,:);
    
    % update in and out ellipsoids
    XmeshIn(parallel,meridian) = pLeft(1);
    XmeshOut(parallel,meridian) = pRight(1);
    YmeshIn(parallel,meridian) = pLeft(2);
    YmeshOut(parallel,meridian) = pRight(2);
    ZmeshIn(parallel,meridian) = pLeft(3);
    ZmeshOut(parallel,meridian) = pRight(3);
    
    % plot profile for debugging
    %figure
    %plot(IcurrCropped);
    
    %plot(Icurr);
    %hold on
    %plot(minIdx,Icurr(minIdx),'rx');
    %plot(maxIdx,Icurr(maxIdx),'rx');
    %title(strcat(num2str(parallel),',',num2str(meridian)))
    %waitforbuttonpress
    %clf
end
end
end

%Icurr = interp3(Xgrid,Ygrid,Zgrid,Ivol,pline(:,1),pline(:,2),pline(:,3));
figure
plot(XmeshIn(parallel, :),YmeshIn(parallel, :));
hold on;
plot(XmeshOut(parallel, :),YmeshOut(parallel, :));
labels = cellstr( num2str([1:meridianN]'));
text(XmeshOut(parallel, :),YmeshOut(parallel, :), labels);


figure
surf(XmeshIn,YmeshIn,dim3-ZmeshIn);
axis equal
axis([1 dim2 1 dim1 1 dim3]);
hold on
xlabel('x');
ylabel('y');
zlabel('z');

hold on
surf(XmeshOut,YmeshOut,dim3-ZmeshOut);

%% Find edges

