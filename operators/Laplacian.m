function L = Laplacian(vtk)
% Relies on the gptoolbox by Alec Jacobson

P = double(vtk.points);
C = double(vtk.cells);
L = cotmatrix(P,C);

end