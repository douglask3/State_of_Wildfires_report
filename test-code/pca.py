import numpy as np
import matplotlib.pyplot as plt

# Define the 2-dimensional vector 'v'
v = np.array([2, 1])  # Replace this with your 2-dimensional vector

# Define the points as a NumPy array
points = np.array([
    [1, 2],
    [3, 2],
    [5, 4],
    [7, 5]
])  # Replace this with your points

# Calculate the projection onto the line perpendicular to 'v'
projection = points - np.dot(points, v)[:, np.newaxis] * v / np.dot(v, v)

# Create a plot
plt.figure(figsize=(6, 6))

# Plot original points
plt.scatter(points[:, 0], points[:, 1], color='blue', label='Original Points')

# Plot the projected points onto the line perpendicular to 'v'
plt.scatter(projection[:, 0], projection[:, 1], color='red', label='Projected Points')

# Plot the vector 'v'
plt.quiver(0, 0, v[0], v[1], angles='xy', scale_units='xy', scale=1, color='green', label='Vector v')

# Set plot labels and legend
plt.xlabel('X-axis')
plt.ylabel('Y-axis')
plt.axhline(0, color='black',linewidth=0.5)
plt.axvline(0, color='black',linewidth=0.5)
plt.legend()
plt.grid(visible=True)

plt.title('Projection of Points onto Line Perpendicular to Vector')
plt.show()

