function [d_norm_closest, d_par_closest] = clstdist_point2segments(P, V)
%CLSTDIST_POINT2SEGMENTSNormal & parallel distance for the closest 2D segment
%
%   [d_norm_closest, d_par_closest] = clstdist_point2segments(P, V)
%
% Inputs:
%   P  – 1×2 array [x, y] of the query point
%   V  – N×2 array of vertex coords; segments are V(i,:)→V(i+1,:)
%
% Outputs:
%   d_norm_closest – scalar perpendicular distance to the closest segment
%   d_par_closest  – scalar parallel distance (zero if projection falls on it)

    % first compute all distances
    [d_norm_all, d_par_all] = dist_point2segments(P, V);

    % Euclidean distance to each segment
    d_total = hypot(d_norm_all, d_par_all);

    % pick the index of the minimum
    [~, idx] = min(d_total);

    % return only that segment’s distances
    d_norm_closest = d_norm_all(idx);
    d_par_closest  = d_par_all(idx);
end
