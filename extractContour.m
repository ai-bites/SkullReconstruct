function [x,y] = extractContour(I, N)

figure;
imshow(I,[]);

msg = strcat('Select ', num2str(N), ' points from the external contour of the cranial vault');
title(msg);

disp(msg);

x= [];y = [];
 hold on;
for count = 1:N,
    [xi,yi] = ginput(1);
    
    plot(xi,yi,'+','color',[ 1.000 0.314 0.510 ],'linewidth',2);        
    x = [x;xi];
    y = [y;yi];
    plot(x,y,'-','color',[ 1.000 0.314 0.510 ],'linewidth',2);
    drawnow;
end;
hold off
  