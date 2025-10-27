def convertCoordinate(original_width, original_height, new_width, new_height, x_resized, y_resized):
    x_original = round(x_resized * original_width / new_width, 13)
    y_original = original_height - round(y_resized * original_height / new_height, 13)

    return x_original, y_original
