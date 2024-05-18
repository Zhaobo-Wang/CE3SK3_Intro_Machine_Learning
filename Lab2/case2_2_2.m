% first load the quadratic surface data
data = load('/MATLAB Drive/SK_lab2/implicit_surface_2.mat');

% extract noisy from noisy data
noisy = data.noisy_observations;

% Get each array of x, y, z
x = noisy(:, 1); % every rows from the first column
y = noisy(:, 2); % every rows from the second column
z = noisy(:, 3); % every rows from the third column

% sum 4
% combined x,y,z power of 4
totalSum_x4 = 0;
totalSum_y4 = 0;
totalSum_z4 = 0;
% combined x,y,z power of 3
totalSum_x3y = 0;
totalSum_xy3 = 0;

totalSum_x3z = 0;
totalSum_xz3 = 0;

totalSum_y3z = 0;
totalSum_yz3 = 0;

% combined x,y,z power of 2,2
totalSum_x2y2 = 0;
totalSum_x2z2 = 0;
totalSum_y2z2 = 0;

% combined x,y,z power of 2,1,1
totalSum_x2yz = 0;
totalSum_xyz2 = 0;
totalSum_xy2z = 0;

% sum 3
% combined x,y sum power of 3
totalSum_x3 = 0;
totalSum_y3 = 0;
totalSum_z3 = 0;
% combined x,y,z power of 2,1
totalSum_xy2 = 0;
totalSum_x2y = 0;
totalSum_xz2 = 0;
totalSum_x2z = 0;
totalSum_y2z = 0;
totalSum_yz2 = 0;
% combined x,y,z power of 1,1,1
totalSum_xyz = 0;

%sum2
% combined x,y sum power of 2
totalSum_x2 = 0;
totalSum_y2 = 0;
totalSum_z2 = 0;
% combined x,y sum power of 1,1
totalSum_xz = 0;
totalSum_yz = 0;
totalSum_xy = 0;

%sum1
totalSum_x = 0;
totalSum_y = 0;
totalSum_z = 0;

for i = 1:length(x) 

    % combined x,y sum power of 1
    totalSum_x = totalSum_x + x(i);
    totalSum_y = totalSum_y + y(i);
    totalSum_z = totalSum_z + z(i);

    % combined x,y sum power of 2
    totalSum_x2 = totalSum_x2 + x(i).^2;
    totalSum_y2 = totalSum_y2 + y(i).^2;
    totalSum_z2 = totalSum_z2 + z(i).^2;
    totalSum_xy = totalSum_xy + x(i)*y(i); 
    totalSum_yz = totalSum_yz + y(i)*z(i); 
    totalSum_xz = totalSum_xz + x(i)*z(i); 


    % combined x,y sum power of 3
    totalSum_x3 = totalSum_x3 + x(i).^3;
    totalSum_y3 = totalSum_y3 + y(i).^3;  
    totalSum_z3 = totalSum_z3 + z(i).^3;

    totalSum_xy2 = totalSum_xy2 + x(i)*y(i).^2;
    totalSum_x2y = totalSum_x2y + x(i).^2*y(i);
    totalSum_x2z = totalSum_x2z + x(i).^2*z(i);
    totalSum_xz2 = totalSum_xz2 + x(i)*z(i).^2;
    totalSum_y2z = totalSum_y2z + y(i).^2*z(i);
    totalSum_yz2 = totalSum_yz2 + y(i)*z(i).^2;

    totalSum_xyz = totalSum_xyz + x(i)*y(i)*z(i);


    % combined x,y sum power of 4
    totalSum_x4 = totalSum_x4 + x(i).^4;
    totalSum_y4 = totalSum_y4 + y(i).^4;
    totalSum_z4 = totalSum_z4 + z(i).^4;

    totalSum_x3y = totalSum_x3y + x(i).^3*y(i);
    totalSum_xy3 = totalSum_xy3 + x(i)*y(i).^3;
    totalSum_y3z = totalSum_y3z + y(i).^3*z(i);
    totalSum_yz3 = totalSum_yz3 + y(i)*z(i).^3;
    totalSum_x3z = totalSum_x3z + x(i).^3*z(i);
    totalSum_xz3 = totalSum_xz3 + x(i)*z(i).^3;

    totalSum_x2y2 = totalSum_x2y2 + x(i).^2*y(i).^2;
    totalSum_y2z2 = totalSum_y2z2 + y(i).^2*z(i).^2;
    totalSum_x2z2 = totalSum_x2z2 + x(i).^2*z(i).^2;

    totalSum_x2yz = totalSum_x2yz + x(i).^2*y(i)*z(i);
    totalSum_xyz2 = totalSum_xyz2 + x(i)*y(i)*z(i).^2;
    totalSum_xy2z = totalSum_xy2z + x(i)*y(i).^2*z(i);

end


left_matrix = [totalSum_x4,totalSum_x2y2,totalSum_x2z2,totalSum_x3y,totalSum_x2yz,totalSum_x3z,totalSum_x3,totalSum_x2y,totalSum_x2z,totalSum_x2;
    totalSum_x2y2, totalSum_y4,totalSum_y2z2,totalSum_xy3,totalSum_y3z,totalSum_xy2z,totalSum_xy2,totalSum_y3,totalSum_y2z,totalSum_y2;
    totalSum_x2z2,totalSum_y2z2,totalSum_z4,totalSum_xyz2,totalSum_yz3,totalSum_xz3,totalSum_xz2,totalSum_yz2,totalSum_z3,totalSum_z2;
    totalSum_x3y,totalSum_xy3,totalSum_xyz2,totalSum_x2y2,totalSum_xy2z,totalSum_x2yz,totalSum_x2y,totalSum_xy2,totalSum_xyz,totalSum_xy;
    totalSum_x2yz,totalSum_y3z,totalSum_yz3,totalSum_xy2z,totalSum_y2z2,totalSum_xyz2,totalSum_xyz,totalSum_y2z,totalSum_yz2,totalSum_yz;
    totalSum_x3z,totalSum_xy2z,totalSum_xz3,totalSum_x2yz,totalSum_xyz2,totalSum_x2z2,totalSum_x2z,totalSum_xyz,totalSum_xz2,totalSum_xz;
    totalSum_x3,totalSum_xy2,totalSum_xz2,totalSum_x2y,totalSum_xyz,totalSum_x2z,totalSum_x2,totalSum_xy,totalSum_xz,totalSum_x;
    totalSum_x2y,totalSum_y3,totalSum_yz2,totalSum_xy2,totalSum_y2z,totalSum_xyz,totalSum_xy,totalSum_y2,totalSum_yz,totalSum_y;
    totalSum_x2z,totalSum_y2z,totalSum_z3,totalSum_xyz,totalSum_yz2,totalSum_xz2,totalSum_xz,totalSum_yz,totalSum_z2,totalSum_z;
    totalSum_x2,totalSum_y2,totalSum_z2,totalSum_xy,totalSum_yz,totalSum_xz,totalSum_x,totalSum_y,totalSum_z,500];


right_matrix = [0;0;0;0;0;0;0;0;0;0];
[v,D] = eig(left_matrix);

[~, smallest_index] = min(diag(D));
estimated_value = v(:,smallest_index);
