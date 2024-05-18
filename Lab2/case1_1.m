% Load the .mat file
data = load('quadratic_surface.mat');

% Extract the noisy_observations array from the loaded data
noisy_observations = data.noisy_observations;

% Separate the noisy_observations array into three arrays: x, y, and z
x = noisy_observations(:, 1); % All rows, first column
y = noisy_observations(:, 2); % All rows, second column
z = noisy_observations(:, 3); % All rows, third column
totalSum_x = 0;
totalSum_y = 0;
totalSum_z = 0;
totalSum_x4 = 0;
totalSum_x3y = 0;
totalSum_xy3 = 0;
totalSum_x2y2 = 0;
totalSum_y4 = 0;

totalSum_x3 = 0;
totalSum_xy2 = 0;
totalSum_x2y = 0;
totalSum_y3 = 0;

totalSum_x2 = 0;
totalSum_y2 = 0;
totalSum_xy = 0;

totalSum_xz = 0;
totalSum_yz = 0;
totalSum_xyz = 0;
totalSum_x2z = 0;
totalSum_y2z = 0;

for i = 1:length(x) 
    %power of 1
    totalSum_x = totalSum_x + x(i);
    totalSum_y = totalSum_y + y(i);
    totalSum_z = totalSum_z + z(i);
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


left_matrix = [totalSum_x4,totalSum_x2y2,totalSum_x3y,totalSum_x3,totalSum_x2y,totalSum_x2;
    totalSum_x2y2, totalSum_y4,totalSum_xy3,totalSum_xy2,totalSum_y3,totalSum_y2;
    totalSum_x3y,totalSum_xy3,totalSum_x2y2,totalSum_x2y,totalSum_xy2,totalSum_xy;
    totalSum_x3,totalSum_xy2,totalSum_x2y,totalSum_x2,totalSum_xy,totalSum_x;
    totalSum_x2y,totalSum_y3,totalSum_xy2,totalSum_xy,totalSum_y2,totalSum_y;
    totalSum_x2,totalSum_y2,totalSum_xy,totalSum_x,totalSum_y,900];

right_matrix = [totalSum_x2z;totalSum_y2z;totalSum_xyz;totalSum_xz;totalSum_yz;totalSum_z];

% Solve Ax = b for x without computing the inverse
% x = A\b;
Results = left_matrix\right_matrix;


% Define the quadratic function with the estimated parameters

A = 1.0018;
B = 1.9955;
C = 3.0001;
D = -0.0299;
E = -0.0236;
F = 1.0042;

quadFunc = @(x, y) A*x.^2 + B*y.^2 + C*x.*y + D*x + E*y + F;

% Define the bounds of the integration
% Find from given data
xmin = 0.0007; 
xmax = 3.0096; 
ymin = 0.0001; 
ymax = 3.0093; 

% Perform the double integral
V = integral2(quadFunc, xmin, xmax, ymin, ymax);

% Display the estimated volume
disp(['The estimated volume is: ', num2str(V)]);
