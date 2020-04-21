function im = atImage(vtk, at, limits, numSteps, angles, txt)

s = 4;
h = figure;
set(h, 'Visible','off', 'Position',[0 0 s*500 s*500], 'color','w');
visualizeDataOnMesh(vtk, at, limits, numSteps, angles);
% set(gcf,'color','w');
im = getframe(h);
im = im.cdata(s*125:end-s*110,s*125:end-s*110,:);
g = uint8(~(sum(abs(gradient(double(im))),3)~=0));
im = insertText(im, [s*130 s*40], txt, 'FontSize',s*16, 'BoxOpacity',0.0, 'AnchorPoint','CenterBottom');
im = imresize(im.*g, 1/s);
close(h);

end