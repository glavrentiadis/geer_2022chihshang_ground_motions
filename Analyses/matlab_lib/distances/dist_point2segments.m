function [d_norm, d_par] = dist_point2segments(P, V)
%DIST_POINT2SEGMENTS  Normal and parallel distance from a point to 2D segments
%
%   [d_norm, d_par] = dist_point2segments(P, V)
%
% Inputs:
%   P  - 1×2 array [x, y] of the query point
%   V  - N×2 array of vertex coordinates, so segments are V(i,:) → V(i+1,:)
%
% Outputs:
%   d_norm - (N-1)×1 array of signed perpendicular distances:
%            positive if P is to the right of the oriented segment A→B,
%            negative if to the left.
%   d_par  - (N-1)×1 array of “parallel” distances beyond the segment
%            (zero when projection lies within [A, B])
%
% Example:
%   P = [1,2];
%   V = [0,0; 3,0; 3,4];
%   [d_norm, d_par] = dist_point2segments(P, V);

    nV = size(V,1);
    if nV < 2
        error('Need at least two vertices to form segments.');
    end
    nSeg = nV - 1;
    d_norm = zeros(nSeg,1);
    d_par  = zeros(nSeg,1);

    for i = 1:nSeg
        A  = V(i,   :);
        B  = V(i+1, :);
        AB = B - A;
        L2 = AB * AB';
        L  = sqrt(L2);
        
        AP = P - A;
        t  = dot(AP, AB) / L2;
        
        % signed normal distance: positive if P is to the right of A→B
        d_norm(i) = (AP(1)*AB(2) - AP(2)*AB(1)) / L;
        
        % parallel distance beyond segment ends
        if t < 0
            d_par(i) = -t * L;
        elseif t > 1
            d_par(i) =  (t - 1) * L;
        else
            d_par(i) = 0;
        end
    end
end
