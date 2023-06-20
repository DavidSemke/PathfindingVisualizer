import pygame as pg
import random as rand
from node import *

# DIMNS = (ROWS, WIDTH)
# set (40, 600) for what was used during testing
# for demonstration, set (20, 400)
DIMNS = (40, 600)
WINDOW = pg.display.set_mode((DIMNS[1], DIMNS[1]))
pg.display.set_caption("A*, LPA*, and D* Lite Pathfinder")

# integer that determines the likelihood of a node being generated 
# as a barrier
# 2 means 1/2 chance, 3 means 1/3 chance...
BARRIER_RAND_CONST = 4


# returns list of lists containing nodes
def make_grid():
    rows, width = DIMNS
    
    # equal to the width of a single tile in the graph
    gap = width // rows
    grid = []

    for i in range(rows):
        grid.append([])

        for j in range(rows):
            node = Node(i, j, gap, rows)
            grid[i].append(node)

    return grid

# make a copy of the current grid state
def duplicate_grid(grid):
    dup_grid = make_grid()
    prev_start = None
    prev_end = None
    
    for row in grid:
        
        for node in row:
            x, y = node.get_pos()
            dup_node = dup_grid[x][y]
            dup_node.color = node.color
            
            if dup_node.is_start():
                prev_start = dup_node
            
            elif dup_node.is_end():
                prev_end = dup_node

    return dup_grid, prev_start, prev_end


# draw grid lines
def draw_grid():
    rows, width = DIMNS
    gap = width // rows

    for i in range(rows):
        pg.draw.line(
            # window, color, start point, end point
            WINDOW, colors['grey'], (0, i * gap), (width, i * gap)
        )  
        
        for j in range(rows):
            pg.draw.line(
                WINDOW, colors['grey'], (j * gap, 0), (j * gap, width)
            )


def draw(grid):
    WINDOW.fill(colors['white'])
    
    [node.draw(WINDOW) for row in grid for node in row]

    draw_grid()
    pg.display.update()


def get_clicked_pos(pos):
    rows, width = DIMNS
    gap = width // rows
    y, x = pos
    row = y // gap
    col = x // gap

    return row, col

# generate barriers that are solely vis
def generate_barriers(grid, barriers_are_vis):
    
    for row in grid:

        for node in row:

            if rand.randint(1, BARRIER_RAND_CONST) != 1: continue
                
            if barriers_are_vis: 
                node.make_vis_barrier()
            
            else:
                node.make_invis_barrier()
        

# generate a mix of vis and invis barriers
def generate_barriers_mixed(grid):
    
    for row in grid:
        
        for node in row:

            if rand.randint(1, BARRIER_RAND_CONST) != 1: continue
            
            if rand.randint(0, 1) == 1:
                node.make_vis_barrier()
            
            else:
                node.make_invis_barrier()