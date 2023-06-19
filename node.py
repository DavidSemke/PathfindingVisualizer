from pygame import draw

'''The following dict provides colors signifying node status.'''

colors = {
    # signifies closed node
    'red': (225, 0, 0),
    # signifies open node 
    'green': (0, 255, 0),
    # signifies invis barrier 
    'blue': (0, 0, 255),
    # signifies original start node if start node is not static 
    'yellow': (255, 255, 0),
    # signifies regular empty node
    'white': (255, 255, 255),
    # signifies vis barrier
    'black': (0, 0, 0),
    # signifies node which is in the path taken to reach goal
    'purple': (128, 0, 128),
    # signifies start node
    'orange': (255, 165, 0),
    # for drawing lines in graph
    'grey': (128, 128, 128),
    # signifies goal node
    'turquoise': (64, 224, 208) 
}


'''A class for defining nodes in the graph. f-values, rhs-values, 
and g-values are stored separately in dictionaries.'''
class Node:
    def __init__(self, row, col, width, total_rows):
        self.row = row
        self.col = col
        self.x = row * width
        self.y = col * width
        self.color = colors['white']
        self.neighbors = []
        self.width = width
        self.total_rows = total_rows

    def get_pos(self):
        return self.row, self.col

    def is_closed(self):
        return self.color == colors['red']

    def is_open(self):
        return self.color == colors['green']

    def is_vis_barrier(self):
        return self.color == colors['black']

    def is_invis_barrier(self):
        return self.color == colors['blue']

    def is_start(self):
        return self.color == colors['orange']

    def is_original_start(self):
        return self.color == colors['yellow']

    def is_end(self):
        return self.color == colors['turquoise']

    def is_path(self):
        return self.color == colors['purple']

    def reset(self):
        self.color = colors['white']

    def make_start(self):
        self.color = colors['orange']

    def make_original_start(self):
        self.color = colors['yellow']

    def make_closed(self):
        self.color = colors['red']

    def make_open(self):
        self.color = colors['green']

    def make_vis_barrier(self):
        self.color = colors['black']

    def make_invis_barrier(self):
        self.color = colors['blue']

    def make_end(self):
        self.color = colors['turquoise']

    def make_path(self):
        self.color = colors['purple']

    def draw(self, window):
        draw.rect(
            window, self.color, 
            (self.x, self.y, self.width, self.width)
        )

    def update_neighbors(self, grid):
        self.neighbors = []
        # DOWN
        if ( 
            self.row < self.total_rows - 1 
            and not grid[self.row + 1][self.col].is_vis_barrier()
        ):  
            self.neighbors.append(grid[self.row + 1][self.col])

        # UP
        if (
            self.row > 0 
            and not grid[self.row - 1][self.col].is_vis_barrier()
        ):  
            self.neighbors.append(grid[self.row - 1][self.col])

        # RIGHT
        if (
            self.col < self.total_rows - 1 
            and not grid[self.row][self.col + 1].is_vis_barrier()
        ):
            self.neighbors.append(grid[self.row][self.col + 1])

        # LEFT
        if (
            self.col > 0 
            and not grid[self.row][self.col - 1].is_vis_barrier()
        ):  
            self.neighbors.append(grid[self.row][self.col - 1])

    def __lt__(self, _):
        return False