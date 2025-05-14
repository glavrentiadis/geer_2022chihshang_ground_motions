function [rot_mat] = axis_rot(theta)
%axis_rot returns the reference system rotation matrix around z axis
%reference: https://en.wikipedia.org/wiki/Rotation_of_axes_in_two_dimensions
%components: E, N, Z

%convert rotation angle to radians;
theta = deg2rad(theta);

%rotation matrix
rot_mat = [ cos(theta), sin(theta), 0;
           -sin(theta), cos(theta), 0;
            0,          0,          1];

end