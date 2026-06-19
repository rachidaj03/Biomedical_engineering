import numpy as np
import matplotlib.pyplot as plt

def create_maze(dim):
    maze = np.zeros((dim, dim))
    stack = [(0, 0)]
    while stack:
        x, y = stack[-1]
        maze[y, x] = 1
        options = [(x-2, y), (x+2, y), (x, y-2), (x, y+2)]
        options = [(x, y) for x, y in options if 0 <= x < dim and 0 <= y < dim and maze[y, x] == 0]
        if options:
            stack.append(options[np.random.randint(0, len(options))])
        else:
            stack.pop()
    return maze

dim = 21
maze = create_maze(dim)
plt.imshow(maze, cmap='binary')
plt.show()
