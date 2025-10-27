import rasterio
import numpy as np
import matplotlib.pyplot as plt
from PIL import Image

Image.MAX_IMAGE_PIXELS = None  

def SigNSize(input_path):
    img = Image.open(input_path)
    original_width, original_height = img.size

    png_width = original_width
    png_height = original_height
    original_width = None
    original_height = None  

    return original_width, original_height, png_width, png_height
