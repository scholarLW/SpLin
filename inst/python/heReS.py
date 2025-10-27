import rasterio
import numpy as np
import matplotlib.pyplot as plt
from PIL import Image

Image.MAX_IMAGE_PIXELS = None  

def HEReSize(input_path, output_path, dpi=300):
    img = Image.open(input_path)
    original_width, original_height = img.size

    if output_path is not None:
        with rasterio.open(input_path) as src:
            num_bands = src.count
            if num_bands < 3:
                raise ValueError(f"The TIFF image has fewer than 3 bands, with only {num_bands} bands.")


            bands = src.read([1, 2, 3])
            data = np.transpose(bands, (1, 2, 0))

            if data.dtype != np.uint8:
                data = (data.astype(np.float32) / 65535.0 * 255).astype(np.uint8)


            plt.imshow(data)
            plt.axis('off')
            plt.savefig(output_path, dpi=dpi)
            plt.close()

            with Image.open(output_path) as png_img:
                png_width, png_height = png_img.size
    else:
        png_width = original_width
        png_height = original_height
        original_width = None
        original_height = None  

    return original_width, original_height, png_width, png_height
