% CT head test

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
    im = imfilter(imread(strcat('./tiff/',files(i).name)),H);
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
[Xmesh,Ymesh,Zmesh] = ellipsoid(x0,xc,yc,ag,R,bg,meshres);

[XmeshIn,YmeshIn,ZmeshIn] = ellipsoid(x0,xc,yc,ag,R,bg,meshres);
XmeshIn = XmeshIn(1:ceil(meshres/2)+1,:);
YmeshIn = YmeshIn(1:ceil(meshres/2)+1,:);
ZmeshIn = ZmeshIn(1:ceil(meshres/2)+1,:);
XmeshOut = XmeshIn;
YmeshOut = YmeshIn;
ZmeshOut = ZmeshIn;
figure

surf(XmeshIn,YmeshIn,ZmeshIn);
axis equal
axis([1 dim2 1 dim1 1 dim3]);
hold on
xlabel('x');
ylabel('y');
zlabel('z');
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

parallel = 20;
threshPercent = 0.8;

for meridian=10:30
    
    p = [XmeshIn(parallel,meridian),YmeshIn(parallel,meridian),ZmeshIn(parallel,meridian)];
    n = [Nx(parallel,meridian),Ny(parallel,meridian),Nz(parallel,meridian)];

    idx = 1;
    pline = [];
    Icurr = [];
    soiLength = 70;
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
    figure
    plot(Icurr);
    hold on
    plot(minIdx,Icurr(minIdx),'rx');
    plot(maxIdx,Icurr(maxIdx),'rx');
    title(strcat(num2str(parallel),',',num2str(meridian)))
end

%Icurr = interp3(Xgrid,Ygrid,Zgrid,Ivol,pline(:,1),pline(:,2),pline(:,3));
figure

surf(XmeshIn,YmeshIn,ZmeshIn);
axis equal
axis([1 dim2 1 dim1 1 dim3]);
hold on
xlabel('x');
ylabel('y');
zlabel('z');
figure

surf(XmeshOut,YmeshOut,ZmeshOut);

%% Find edges










