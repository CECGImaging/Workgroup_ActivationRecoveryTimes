%{
This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at https://mozilla.org/MPL/2.0/.

Copyright 2019 Steffen Schuler
Institute of Biomedical Engineering
Karlsruhe Institute of Technology
www.ibt.kit.edu
%}

function visualizeDataOnMesh(vtk, data, limits, numSteps, viewAngles)

if nargin < 3
    limits = [min(data) max(data)];
end
if nargin < 4
    numSteps = 20;
end
if nargin < 5
    viewAngles = [0 0];
end

trisurf(vtk.cells, vtk.points(:,1), vtk.points(:,2), vtk.points(:,3), data,  'facecolor','interp', 'edgecolor','none');
axis equal;
camproj('persp');
set(gca,'visible','off');
caxis(limits);
colormap(jet(numSteps));
view(viewAngles(1), viewAngles(2));

end