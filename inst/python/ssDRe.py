import rasterio
import numpy as np
import matplotlib.pyplot as plt
from PIL import Image

Image.MAX_IMAGE_PIXELS = None 

def ssDNAReSize(input_path, output_path, dpi=300):
    img = Image.open(input_path)
    original_width, original_height = img.size

    if output_path is not None:
        with rasterio.open(input_path) as src:
            data = src.read(1)  

        plt.imshow(data, cmap='gray') 
        plt.axis('off')  
        plt.savefig(output_path, bbox_inches='tight', pad_inches=0, dpi=dpi) 

        with Image.open(output_path) as png_img:
            png_width, png_height = png_img.size
    else:
        png_width = original_width
        png_height = original_height
        original_width = None
        original_height = None     

    return original_width, original_height, png_width, png_height