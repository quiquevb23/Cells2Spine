'''
Code to register 2 images to a common space of coordinates, modifying the barcode position too

'''

import cv2
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import json

# Load the histology images
image1 = cv2.imread('/home/quiquevb/Desktop/Data/Yadav/Sample1/GSM6919905_V160-CGND-HRA-02636-C_tissue_hires_image.png')
image2 = cv2.imread('/home/quiquevb/Desktop/Data/Yadav/Sample2/GSM6919906_V160-CGND-HRA-02636-D_tissue_hires_image.png')

#Load the manually selected landmarks
landmarks_sample1 = pd.read_csv('/home/quiquevb/Desktop/Data/Yadav/Sample1_landmarks.csv')
landmarks_sample2 = pd.read_csv('/home/quiquevb/Desktop/Data/Yadav/Sample2_landmarks.csv')

# Load the scale factors from JSON files
with open('/home/quiquevb/Desktop/Data/Yadav/scalefactors/V160-CGND-HRA-02636-C_scalefactors_json.json', 'r') as file:
    scale_factors1 = json.load(file)
with open('/home/quiquevb/Desktop/Data/Yadav/scalefactors/V160-CGND-HRA-02636-D_scalefactors_json.json', 'r') as file:
    scale_factors2 = json.load(file)

# Extract scaling factors
tissue_hires_scalef1 = scale_factors1['tissue_hires_scalef']
tissue_hires_scalef2 = scale_factors2['tissue_hires_scalef']

'''
# Rescale landmark coordinates
landmarks_sample2['X'] *= tissue_hires_scalef2
landmarks_sample2['Y'] *= tissue_hires_scalef2
landmarks_sample1['X'] *= tissue_hires_scalef1
landmarks_sample1['Y'] *= tissue_hires_scalef1
'''

# Convert to NumPy arrays
points1 = landmarks_sample1[['X', 'Y']].values
points2 = landmarks_sample2[['X', 'Y']].values

'''
# Convert to grayscale for better processing
gray1 = cv2.cvtColor(image1, cv2.COLOR_BGR2GRAY)
gray2 = cv2.cvtColor(image2, cv2.COLOR_BGR2GRAY)

# Display the images
plt.figure(figsize=(10, 5))
plt.subplot(1, 2, 1)
plt.title('Sample 1')
plt.imshow(gray1, cmap='gray')
plt.subplot(1, 2, 2)
plt.title('Sample 2')
plt.imshow(gray2, cmap='gray')
plt.show()
'''

# Compute the homography matrix using the manually selected landmarks
H, status = cv2.findHomography(points2, points1, cv2.RANSAC, 5.0)

# Use homography to warp the second image to the first
height, width = image1.shape[:2]
image2_aligned = cv2.warpPerspective(image2, H, (width, height))

# Display the aligned images
plt.figure(figsize=(15, 5))
plt.subplot(1, 2, 1)
plt.title('Original Sample 1')
plt.imshow(cv2.cvtColor(image1, cv2.COLOR_BGR2RGB))

plt.subplot(1, 2, 2)
plt.title('Aligned Sample 2')
plt.imshow(cv2.cvtColor(image2_aligned, cv2.COLOR_BGR2RGB))
plt.show()

# Load barcode positions (assuming CSV format)
barcodes_sample1 = pd.read_csv('/home/quiquevb/Desktop/Data/Yadav/Sample1/GSM6919905_V160-CGND-HRA-02636-C_tissue_positions.csv.gz', header=None)
barcodes_sample2 = pd.read_csv('/home/quiquevb/Desktop/Data/Yadav/Sample2/GSM6919906_V160-CGND-HRA-02636-D_tissue_positions.csv.gz', header=None)

# Scale the barcode positions, 4th and 5th columns are X and Y coordinates
barcodes_sample2[4] *= tissue_hires_scalef2
barcodes_sample2[5] *= tissue_hires_scalef2
barcodes_sample1[4] *= tissue_hires_scalef1
barcodes_sample1[5] *= tissue_hires_scalef1


# Transform barcodes_sample2 using the homography matrix with cv2.perspectiveTransform
barcode_positions_sample2 = np.column_stack([barcodes_sample2[4], barcodes_sample2[5]]).astype(np.float32)
barcode_positions_sample2 = np.expand_dims(barcode_positions_sample2, axis=1)  # Prepare for perspectiveTransform

# Apply the perspective transformation
transformed_barcode_positions = cv2.perspectiveTransform(barcode_positions_sample2, H)
transformed_barcode_positions = transformed_barcode_positions.squeeze()

'''
#scale barcodes of image 2 afterwards
transformed_barcode_positions[4] *= tissue_hires_scalef2
transformed_barcode_positions[5] *= tissue_hires_scalef2
'''

# Create a new DataFrame for the transformed barcodes
transformed_barcodes_df = barcodes_sample2.copy()
transformed_barcodes_df['X'] = transformed_barcode_positions[:, 0]
transformed_barcodes_df['Y'] = transformed_barcode_positions[:, 1]

# Clip barcode positions that fall outside the image boundaries
transformed_barcodes_df = transformed_barcodes_df[
    (transformed_barcodes_df['X'] >= 0) & (transformed_barcodes_df['X'] < width) &
    (transformed_barcodes_df['Y'] >= 0) & (transformed_barcodes_df['Y'] < height)
]

aligned_barcodes_df = transformed_barcodes_df.copy()
aligned_barcodes_df.drop(aligned_barcodes_df.columns[[4, 5]], axis=1, inplace=True)
aligned_barcodes_df['X'] /= tissue_hires_scalef2
aligned_barcodes_df['Y'] /= tissue_hires_scalef2

# Convert 'X' and 'Y' columns to integers
aligned_barcodes_df['X'] = aligned_barcodes_df['X'].astype(int)
aligned_barcodes_df['Y'] = aligned_barcodes_df['Y'].astype(int)

aligned_barcodes_df.insert(4, 'X', aligned_barcodes_df.pop('X'))
aligned_barcodes_df.insert(5, 'Y', aligned_barcodes_df.pop('Y'))

# Save the transformed barcodes
transformed_barcodes_df.to_csv('/home/quiquevb/Desktop/Data/Yadav/Sample2/aligned_scaled_barcodes_sample2.csv', index=False)
# Save the aligned but reverted barcodes for use
aligned_barcodes_df.to_csv('/home/quiquevb/Desktop/Data/Yadav/Sample2/aligned_barcodes_sample2.csv', header=False, index=False)


# Also save original barcodes to csv for comparison of spots
barcodes_sample1.to_csv('/home/quiquevb/Desktop/Data/Yadav/Sample1/scaled_barcodes_sample1.csv', index=False)

# Optionally, visualize the original and transformed barcodes
plt.figure(figsize=(10, 5))
plt.subplot(1, 2, 1)
plt.title('Original Sample 1 Barcodes')
plt.imshow(cv2.cvtColor(image1, cv2.COLOR_BGR2RGB))
plt.scatter(barcodes_sample1[4], barcodes_sample1[5], color='red', s=10)

plt.subplot(1, 2, 2)
plt.title('Transformed Sample 2 Barcodes')
plt.imshow(cv2.cvtColor(image2_aligned, cv2.COLOR_BGR2RGB))
plt.scatter(transformed_barcodes_df['X'], transformed_barcodes_df['Y'], color='blue', s=10)

plt.show()


