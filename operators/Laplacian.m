function L = Laplacian(vtk)
% Relies on the gptoolbox by Alec Jacobson

P = double(vtk.points);
C = double(vtk.cells);
L = massmatrix(P,C,'voronoi') \ cotmatrix(P,C);

end