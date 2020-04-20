function L = Laplacian(vtk)
% Based on the gptoolbox by Alec Jacobson

P = double(vtk.points);
C = double(vtk.cells);

switch size(C,2)
    case 3
        v1 = P(C(:,2),:)-P(C(:,1),:);
        v2 = P(C(:,3),:)-P(C(:,2),:); 
        v3 = P(C(:,1),:)-P(C(:,3),:);

        l = sqrt([sum(v1.^2,2) sum(v2.^2,2) sum(v3.^2,2)]);
        d = heron(l)/12;
        
        v = -[dot(v2,v3,2)./d dot(v3,v1,2)./d dot(v1,v2,2)./d];
        
        L = sparse(C(:,[1 2 3]), C(:,[2 3 1]), v, size(P,1), size(P,1));
        L = L+L';
        L = L-diag(sum(L,2));
        
        % EXPLANATION:
        % v(:,1) = -cot(v2,v3)/(2*vertexArea)
        % cot(v2,v3) = dot(v2,v3)/(2*triangleArea)
        % vertexArea = triangleArea/3;
        % triangleArea = sqrt(12*d)/4;
        
    case 4
        L = cotmatrix(P,C);
        
    otherwise
        error('Invalid number of columns in cell list.');
end

end

function p = heron(l)

% Heron's formula with Kahan summation:
% triangleArea = sqrt(p)/4

l = sort(l, 2, 'descend');

p = ( l(:,1) + (l(:,2)+l(:,3)) ) .* ...
    ( l(:,3) - (l(:,1)-l(:,2)) ) .* ...
    ( l(:,3) + (l(:,1)-l(:,2)) ) .* ...
    ( l(:,1) + (l(:,2)-l(:,3)) );

end