import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Function to normalize a vector
def normalize(vector):
    return vector / np.linalg.norm(vector)

# Function to project points onto a plane perpendicular to the line defined by the gradient
def project_onto_plane(point, gradient):
    gradient_unit = normalize(gradient)
    # Generating a random vector orthogonal to the gradient vector
    v = np.random.randn(len(gradient))
    orthogonal_vector = v - np.dot(v, gradient_unit) * gradient_unit
    orthogonal_vector = normalize(orthogonal_vector)
    
    # Projecting point onto the plane perpendicular to the line
    projection = point - np.dot(point, orthogonal_vector) * orthogonal_vector
    return projection

# Generating random scatter points in 3D
np.random.seed(42)
scatter_points = np.random.randn(100, 3)  # 100 random 3D points
gradient_vector = np.array([1, 2, 1])     # Example gradient vector

# Projecting scatter points onto the plane perpendicular to the line defined by the gradient vector
collapsed_points = np.array([project_onto_plane(point, gradient_vector) for point in scatter_points])

# Plotting the original scatter points in 3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.scatter(scatter_points[:, 0], scatter_points[:, 1], scatter_points[:, 2], c='blue', label='Original Points')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

# Plotting the collapsed points onto the plane perpendicular to the line
ax.scatter(collapsed_points[:, 0], collapsed_points[:, 1], collapsed_points[:, 2], c='red', label='Collapsed Points')
ax.legend()
plt.title('Scatter points and their projections onto the plane perpendicular to the gradient vector')
plt.show()

