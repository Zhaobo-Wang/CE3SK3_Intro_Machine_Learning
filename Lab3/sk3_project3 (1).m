function project3_code
    % load data 
    data = load('castle.mat');
    
    % load kernel and given image
    mykernel = data.kernel_weights;
    given_image = data.blurred_image;

    % given image size m/n
    [m, n] = size(given_image);
    A = construct_Convolution_Matrix(mykernel, m, n);

    given_image = double(given_image);    % double precision

    my_deblurred_image = my_deblurImage(A, given_image(:), m, n);% this is my deblurred image function

    figure;
    subplot(1, 2, 1); 
    imshow(given_image, []); title('Blurred castle');
    subplot(1, 2, 2); 
    imshow(my_deblurred_image, []); title('Deblurred castle');
     % Display two image
end

%construct the convolution matrix
function A = construct_Convolution_Matrix(blur_kernel, m, n)
    A = sparse(m*n, m*n);
    Center_of_Kernel = floor((size(blur_kernel) + 1) / 2); %find the center of the kernel

    for i = 1:m % loop for the row
        for j = 1:n % loop for the column 
            my_index = (j-1)*m + i; % find the index
            for p = 1:size(blur_kernel, 1) % loop kernel row
                for q = 1:size(blur_kernel, 2) % loop kernel column
                    rowOffset = p - Center_of_Kernel(1); % offset
                    colOffset = q - Center_of_Kernel(2);
                    newRow = i + rowOffset;
                    newCol = j + colOffset;
                    if newRow >= 1 && newRow <= m && newCol >= 1 && newCol <= n
                        idx_neighbor = (newCol-1)*m + newRow;
                        A(my_index, idx_neighbor) = blur_kernel(p, q);
                    end
                end
            end
        end
    end
    A = sparse(A); 
end

function deblurred_image = my_deblurImage(A, blurred_image_v, m, n)
    A_inverse = Lu_inverse(A); % define a inverse by my self
    deblurred_image_v = A_inverse * blurred_image_v;
    deblurred_image = reshape(deblurred_image_v, m, n);
end

function A_inverse = Lu_inverse(A)
    [L, U] = myLU(A);
    A_inverse = zeros(size(A));
    n = size(A, 1);
    for i = 1:n
        b = zeros(n, 1);
        b(i) = 1;
        d = sub_forward(L, b); %call sub_forward
        x = sub_backward(U, d); %call sub_backward
        A_inverse(:, i) = x;
    end
end

function x = sub_forward(L, b)
    n = length(b);
    x = zeros(n, 1);
    for i = 1:n
        x(i) = (b(i) - L(i, 1:i-1) * x(1:i-1)) / L(i, i); 
    end
end

function x = sub_backward(U, b)
    n = length(b);
    x = zeros(n, 1);
    for i = n:-1:1
        x(i) = (b(i) - U(i, i+1:end) * x(i+1:end)) / U(i, i);
    end
end

function [L, U] = myLU(A) %define my own lu function
    n = size(A, 1);
    L = eye(n);
    U = zeros(n);
    for j = 1:n
        for i = 1:j
            U(i, j) = A(i, j) - L(i, 1:i-1) * U(1:i-1, j); %my lu get u decomposion
        end
        for i = j+1:n
            L(i, j) = (A(i, j) - L(i, 1:j-1) * U(1:j-1, j)) / U(j, j); %my lu get l decomposion
        end
    end
end

