% Load the .mat file
data = load('cubic_surface.mat');

% Extract the noisy_observations array from the loaded data
noisy_observations = data.noisy_observations;

% Separate the noisy_observations array into three arrays: x, y, and z
x = noisy_observations(:, 1); % All rows, first column
y = noisy_observations(:, 2); % All rows, second column
z = noisy_observations(:, 3); % All rows, third column
totalSum_x = 0;
totalSum_y = 0;
totalSum_z = 0;

%power of 6
totalSum_x6 = 0;
totalSum_y6 = 0;
totalSum_x5y = 0;
totalSum_xy5 = 0;
totalSum_x4y2 = 0;
totalSum_x2y4 = 0;
totalSum_x3y3 = 0;

%power of 5
totalSum_x5 = 0;
totalSum_y5 = 0;
totalSum_x4y = 0;
totalSum_xy4 = 0;
totalSum_x3y2 = 0;
totalSum_x2y3 = 0;

%power of 4
totalSum_x4 = 0;
totalSum_y4 = 0;
totalSum_x3y = 0;
totalSum_xy3 = 0;
totalSum_x2y2 = 0;

%power of 3
totalSum_x3 = 0;
totalSum_y3 = 0;
totalSum_xy2 = 0;
totalSum_x2y = 0;

%power of 2
totalSum_x2 = 0;
totalSum_y2 = 0;
totalSum_xy = 0;

totalSum_xz = 0;
totalSum_yz = 0;
totalSum_xyz = 0;
totalSum_x2z = 0;
totalSum_y2z = 0;
totalSum_x2yz = 0;
totalSum_xy2z = 0;
totalSum_x3z = 0;
totalSum_y3z = 0;

for i = 1:length(x) 
    %power of 1
    totalSum_x = totalSum_x + x(i);
    totalSum_y = totalSum_y + y(i);
    totalSum_z = totalSum_z + z(i);

    %power of 6
    totalSum_x6 = totalSum_x6 + x(i).^6;
    totalSum_y6 = totalSum_y6 + y(i).^6;
    totalSum_x5y = totalSum_x5y + x(i).^5*y(i);
    totalSum_xy5 = totalSum_xy3 + x(i)*y(i).^5;
    totalSum_x4y2 = totalSum_x4y2 + x(i).^4*y(i).^2;
    totalSum_x2y4 = totalSum_x2y4 + x(i).^2*y(i).^4;
    totalSum_x3y3 = totalSum_x3y3 + x(i).^3*y(i).^3;

    %power of 5
    totalSum_x5 = totalSum_x5 + x(i).^5;
    totalSum_y5 = totalSum_y5 + y(i).^5;
    totalSum_x4y = totalSum_x4y + x(i).^4*y(i);
    totalSum_xy4 = totalSum_xy4 + x(i)*y(i).^4;
    totalSum_x3y2 = totalSum_x3y2 + x(i).^3*y(i).^2;
    totalSum_x2y3 = totalSum_x2y3 + x(i).^2*y(i).^3;
    
    %power of 4
    totalSum_x4 = totalSum_x4 + x(i).^4;
    totalSum_x3y = totalSum_x3y + x(i).^3*y(i);
    totalSum_xy3 = totalSum_xy3 + x(i)*y(i).^3;
    totalSum_x2y2 = totalSum_x2y2 + x(i).^2*y(i).^2;
    totalSum_y4 = totalSum_y4 + y(i).^4;

    %power of 3
    totalSum_x3 = totalSum_x3 + x(i).^3;
    totalSum_xy2 = totalSum_xy2 + x(i)*y(i).^2;
    totalSum_x2y = totalSum_x2y + x(i).^2*y(i);
    totalSum_y3 = totalSum_y3 + y(i).^3;
    %power of 2
    totalSum_x2 = totalSum_x2 + x(i).^2;
    totalSum_y2 = totalSum_y2 + y(i).^2;
    totalSum_xy = totalSum_xy + x(i)*y(i);   
    %relate to z 
    totalSum_xz = totalSum_xz + x(i)*z(i);
    totalSum_yz = totalSum_yz + y(i)*z(i);
    totalSum_xyz = totalSum_xyz + x(i)*y(i)*z(i);
    totalSum_x2z = totalSum_x2z + x(i).^2*z(i);
    totalSum_y2z = totalSum_y2z + y(i).^2*z(i);
end


left_matrix = [totalSum_x6,totalSum_x3y3,totalSum_x5y,totalSum_x4y2,totalSum_x5,totalSum_x3y2,totalSum_x4y,totalSum_x4,totalSum_x3y,totalSum_x3;
    totalSum_x3y3, totalSum_y6,totalSum_x2y4,totalSum_xy5,totalSum_x2y3,totalSum_y5,totalSum_xy4,totalSum_xy3,totalSum_y4,totalSum_y3;
    totalSum_x5y,totalSum_x2y4,totalSum_x4y2,totalSum_x3y3,totalSum_x4y,totalSum_x2y3,totalSum_x3y2,totalSum_x3y,totalSum_x2y2,totalSum_x2y;
    totalSum_x4y2,totalSum_xy5,totalSum_x3y3,totalSum_x2y4,totalSum_x3y2,totalSum_xy4,totalSum_x2y3,totalSum_x2y2,totalSum_xy3,totalSum_xy2;
    totalSum_x5,totalSum_x2y3,totalSum_x4y,totalSum_x3y2,totalSum_x4,totalSum_x2y2,totalSum_x3y,totalSum_x3,totalSum_x2y,totalSum_x2;
    totalSum_x3y2,totalSum_y5,totalSum_x2y3,totalSum_xy4,totalSum_x2y2,totalSum_y4,totalSum_xy3,totalSum_xy2,totalSum_y3,totalSum_y2;
    totalSum_x4y,totalSum_xy4,totalSum_x3y2,totalSum_x2y3,totalSum_x3y,totalSum_xy3,totalSum_x2y2,totalSum_x2y,totalSum_xy2,totalSum_xy;
    totalSum_x4,totalSum_xy3,totalSum_x3y,totalSum_x2y2,totalSum_x3,totalSum_xy2,totalSum_x2y,totalSum_x2,totalSum_xy,totalSum_x;
    totalSum_x3y,totalSum_y4,totalSum_x2y2,totalSum_xy3,totalSum_x2y,totalSum_y3,totalSum_xy2,totalSum_xy,totalSum_y2,totalSum_y;
    totalSum_x3,totalSum_y3,totalSum_x2y,totalSum_xy2,totalSum_x2,totalSum_y2,totalSum_xy,totalSum_x,totalSum_y,1600];

right_matrix = [totalSum_x3z;totalSum_y3z;totalSum_x2yz;totalSum_xy2z;totalSum_x2z;totalSum_y2z;totalSum_xyz;totalSum_xz;totalSum_yz;totalSum_z];

% Solve Ax = b for x without computing the inverse
% x = A\b;
Results = left_matrix\right_matrix;

% Define the quadratic function with the estimated parameters
A = 1.0001;
B = -0.00002829;
C = -0.0001198;
D = -3.0007;
E = 0.4852;
F = -0.4848;
G = 0.02888;
H = 0.9949;
I = 0.005009;
J = 0.9994;

cubicFunc = @(x, y) A*x.^3 + B*y.^3 + C*x.^2.*y + D*x.*y.^2 + E*x.^2 + F*y.^2 + G*x.*y + H*x + I*y + J;

% Define the bounds of the integration
% Find from given data
xmin = -2.9999; 
xmax = 3.01; 
ymin = -2.9996; 
ymax = 3.0098; 

% Perform the double integral
V = integral2(cubicFunc, xmin, xmax, ymin, ymax);

% Display the estimated volume
disp(['The estimated volume is: ', num2str(V)]);

[x, y] = meshgrid(linspace(xmin, xmax, 50), linspace(ymin, ymax, 50));

% Evaluate the cubic function on the grid
z = A*x.^3 + B*y.^3 + C*x.^2.*y + D*x.*y.^2 + E*x.^2 + F*y.^2 + G*x.*y + H*x + I*y + J;

% Plot the surface
figure;
surf(x, y, z);
