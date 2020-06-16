function im = atImage(vtk, at, limits, numSteps, angles, txt)

res = 500; % resolution
s = 2;     % scale factor, influences the thickness of contour lines
h = figure;
set(h, 'Visible','off', 'Position',[0 0 s*res s*res], 'color','w');
visualizeDataOnMesh(vtk, at, limits, numSteps, angles);
im = getframe(h);
im = cropWhiteBackground(im.cdata);
g = uint8(~(sum(abs(gradient(double(im))),3)~=0));
im = im.*g;
if nargin > 5
    im = [repmat(uint8(255), round(s*res/8), size(im,2), size(im,3)); im];
    im = insertText(im, [size(im,2)/2 s*res/7.5], txt, 'FontSize',round(s*res/15.5), 'BoxOpacity',0.0, 'AnchorPoint','CenterBottom');
end
im = imresize(im, 1/s);
close(h);

end

function im = cropWhiteBackground(im)

x1 = find(sum(sum(im,3),1) < 3*255*size(im,1), 1) - 1;
x2 = find(sum(sum(im,3),1) < 3*255*size(im,1), 1, 'last') + 1;
x2 = x2 + mod(x2-x1+1,2);
y1 = find(sum(sum(im,3),2) < 3*255*size(im,2), 1) - 1;
y2 = find(sum(sum(im,3),2) < 3*255*size(im,2), 1, 'last') + 1;
y2 = y2 + mod(y2-y1+1,2);

xpad = round((x2-x1)/40);
ypad = round((y2-y1)/40);
x1 = max(x1-xpad, 1);
x2 = min(x2+xpad, size(im,2));
y1 = max(y1-ypad, 1);
y2 = min(y2+ypad, size(im,1));

im = im(y1:y2,x1:x2,:);

end