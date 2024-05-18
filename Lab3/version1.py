import numpy as np
from scipy.signal import convolve2d
from scipy.linalg import lu
from PIL import Image
from scipy.io import loadmat
import cv2
# Load data from .mat file
mat_data = loadmat('butterfly.mat')

# Get 'kernel' variable
kernel = mat_data['kernel_weights']

# Get 64x64 matrix
blurred_image = mat_data['blurred_image']

# Convert to uint8 to ensure correct data type for an image
image_data = np.uint8(blurred_image)
image = Image.fromarray(image_data)

# Convert PIL Image to NumPy array to use its shape
image_array = np.array(image)

# Now use image_array.shape for your operations
image_shape = image_array.shape
#image_jpg = cv2.imread('D:/McMaster/3SK3/Lab3/but_64_64.jpg')
#image_gray = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)

blurred_image = convolve2d(image, kernel, mode='same', boundary='fill', fillvalue=0)

def construct_conv_matrix(image_shape, kernel):
    m, n = image_shape
    print(m)
    print(n)  
    k_m, k_n = kernel.shape  

    center_k_m = k_m // 2
    center_k_n = k_n // 2

    # 创建一个零矩阵A
    A = np.zeros((m*n, m*n))

    # 为每个像素填充卷积矩阵A
    for i in range(m):
        for j in range(n):
            for k in range(k_m):
                for l in range(k_n):
                    # 计算卷积核对应于图像中的位置
                    row = i + (k - center_k_m)
                    col = j + (l - center_k_n)

                    # 检查边界
                    if 0 <= row < m and 0 <= col < n:
                        # 在矩阵A中设置相应的值
                        A[i*n+j, row*n+col] = kernel[k, l]

    return A

def myLUDecomp(A):
    n = A.shape[0]
    U = A.copy().astype(float)
    L = np.eye(n, dtype=float)
    P = np.eye(n, dtype=float)

    for i in range(n-1):
        # 查找最大的主元
        pivotRow = findMaxPivotRow_new(U, i, n)
        
        # 交换行
        if pivotRow != i:
            U[[i, pivotRow]] = U[[pivotRow, i]]
            P[[i, pivotRow]] = P[[pivotRow, i]]
            if i > 0:  # 交换L矩阵的部分行
                L[[i, pivotRow], :i] = L[[pivotRow, i], :i]

        # 检查除数是否为零，如果是，则跳过当前步骤
        if U[i, i] == 0:
            continue

        # 向量化的更新
        L[i+1:n, i] = U[i+1:n, i] / U[i, i]
        U[i+1:n, i:] -= np.outer(L[i+1:n, i], U[i, i:])

    return L, U, P

def findMaxPivotRow_new(U, currentRow, totalRows):
    return np.argmax(np.abs(U[currentRow:totalRows, currentRow])) + currentRow
 

A = construct_conv_matrix(image_shape, kernel)
print("A:",A.shape)
b = blurred_image.ravel()
L, U, P = myLUDecomp(A)
y = np.linalg.solve(L, b)
x = np.linalg.solve(U, y)
deblurred_image = x.reshape(image_shape)

screen_res = 1280,720
# Calculate the scale and the new window dimensions
scale_width = screen_res[0] / deblurred_image.shape[1]
scale_height = screen_res[1] / deblurred_image.shape[0]
scale = min(scale_width, scale_height)
window_width = int(deblurred_image.shape[1] * scale)
window_height = int(deblurred_image.shape[0] * scale)

# Resize images for display
# Convert PIL Image to NumPy array for OpenCV operations
image_cv = np.array(image)
resized_original = cv2.resize(image_cv, (window_width, window_height))
resized_blurred = cv2.resize(blurred_image, (window_width, window_height))
resized_deblurred = cv2.resize(deblurred_image, (window_width, window_height))

resized_blurred = np.clip(resized_blurred, 0, 255).astype('uint8')
resized_deblurred = np.clip(resized_deblurred, 0, 255).astype('uint8')

# Display resized images
cv2.imshow('Original Image', resized_original)
cv2.imshow('Blurred Image', resized_blurred)
cv2.imshow('Deblurred Image', resized_deblurred)
cv2.waitKey(0)
cv2.destroyAllWindows()


