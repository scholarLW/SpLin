from cCo import convertCoordinate as conv

def adjust_shape_coordinates(shapes, neww, newh, oldw, oldh):
    adjusted_shapes = []
    for shape in shapes:
        adjusted_points = [conv(oldw, oldh, neww, newh, x, y) for x, y in shape['points']]
        shape['points'] = adjusted_points
        adjusted_shapes.append(shape)
    return adjusted_shapes
