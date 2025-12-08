from PIL import Image
import numpy as np

def image_to_binary_txt(image_path, output_txt):
    # Open the image and convert it to grayscale
    img = Image.open(image_path).convert('L')

    # Obtain image size (width, height)
    width, height = img.size
    print(f"Image pixel size:{width}（width）x {height}（height）")
    
    # Convert image data to numpy array
    img_array = np.array(img)
    
    # Threshold processing
    # Set grayscale values greater than 128 as 1 (white) and the rest as 0 (black)
    binary_array = (img_array > 128).astype(int)
    
    # Flatten the matrix into a single column and save as txt
    with open(output_txt, 'w') as f:
        for value in binary_array.flatten():
            f.write(f"{value}\n")

image_to_binary_txt('Ex_comparison3.png', 'output.txt')
