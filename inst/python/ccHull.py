from scipy.spatial import ConvexHull

def create_convex_hull(points):
    hull = ConvexHull(points)
    hull_points = points[hull.vertices]
    return hull_points
