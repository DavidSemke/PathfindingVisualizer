import pygame as pg
import random as rand
from node import *
from a_star import perform_a_star
from lpa_star import perform_lpa_star
from d_star_lite import perform_d_star_lite

'''The following are the available command keys for this pathfinding 
program.'''

# left mouse click -
#   The first click places the start node.
#   The second click places the goal node.
#   Clicks after these nodes are placed will create barriers of the 
# type determined by whatever type has been
#   selected using the 'z' key toggle.

# ** Note that the program will not start until a start and end node 
# have been selected.

# right mouse click - clears a node of any unique status (becomes an empty node).

# Key 'z' - toggle between vis (visible) and invis (invisible) 
# barriers

# Key 'g' - this will generate a graph with barriers; the number of 
# barriers corresponds to the set obstacle density
# of the program, which is by default 25%. Toggle the default 
# obstacle density by changing the
# BARRIER_RAND_CONST in the program's code (below). In the program, 
# you can toggle the type of barrier to be
# generated by pressing the key 'z'. Pressing 'z' toggles between 
# invis barrier generation and vis barrier
# generation.

# Key 'm' - this will generate a graph with barriers like when 
# pressing key 'g', except in this case a mixture of invis and 
# vis barriers will be generated. When a node is
# selected to be a barrier during generation, there is a 50-50 chance 
# of an invis or vis barrier being
# generated.

# Key 'c' - clear the grid, meaning all nodes become empty nodes.

# Key 's' - clear the grid of all nodes that are paths, open, or 
# closed. These node statuses are described below.

# Key 'b' - restore the graph to its state before the last search 
# algorithm was initiated.

# Key '1' - initiate A* without travel.

# Key '2' - initiate A* with travel (compute shortest path and 
# traverse to goal position).

# Key '3' - initiate LPA*.

# Key '4' - initiate D* Lite.


'''The following are descriptions of all node statuses. Corresponding 
colors are defined further down.'''

# Closed node - In the case of A*, a closed node is a node that has 
# been visited due to its f-value being the smallest.
# For LPA* and D* Lite, however, a closed node is a node that is 
# locally consistent, where a node is locally consistent when its 
# rhs-value is equal to its g-value.

# Open node - In all algorithms, the open node is a node that is in 
# the open set, which is a priority queue.

# Path node - A path node signifies that the node is part of the path 
# to the goal node (may not be part of the shortest
# path, in the case of locating and moving to the goal).

# invis Barrier - A node which is not able to be traversed to, 
# but the agent is not yet aware of this.

# vis barrier - A node which is not able to be traversed to, and 
# the agent is aware of this.

# Original start node - This node is simply the starting point for an 
# algorithm which seeks to move the agent
#   from the start node to the goal node.

# Empty node - A node which has no unique status.

# Start node - The node at which the agent starts calculating the 
# shortest path. This node changes when the algorithm
#   demands the agent to traverse from the start node to the end node.

# End (goal) node - The node that the agent is attempting to find a 
# shortest path to, and perhaps traverse to.


'''Global variables.'''

# DIMNS = (ROWS, WIDTH)
# set (40, 600) for what was used during testing
# for demonstration, set (20, 400)
DIMNS = (40, 600)
WINDOW = pg.display.set_mode((DIMNS[1], DIMNS[1]))
pg.display.set_caption("A*, LPA*, and D* Lite Pathfinder")

# determines how likely it is for a node to be generated as a barrier
# 2 means 1/2 chance, 3 means 1/3 chance...
# 1/4 chance per node means the obstacle density is (1/4 * 100) %
# the barrier random constant should be an integer
BARRIER_RAND_CONST = 4


# this defines the grid as a list, which contains inner-lists 
# containing Node instances
# each inner list contains one row's worth of nodes
def make_grid(dimns):
    rows, width = dimns
    
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
def duplicate_grid(grid, dimns):
    dup_grid = make_grid(dimns)
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
def draw_grid(window, dimns):
    rows, width = dimns
    gap = width // rows

    for i in range(rows):
        pg.draw.line(
            # window, color, start point, end point
            window, colors['grey'], (0, i * gap), (width, i * gap)
        )  
        
        for j in range(rows):
            pg.draw.line(
                window, colors['grey'], (j * gap, 0), (j * gap, width)
            )


def draw(window, grid, dimns):
    window.fill(colors['white'])
    
    [node.draw(window) for row in grid for node in row]

    draw_grid(window, dimns)
    pg.display.update()


def get_clicked_pos(pos, dimns):
    rows, width = dimns
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


def handle_left_click(env, node, barriers_are_vis):
    grid, start, end = env
    
    if not start and not end:
        start = node
        start.make_start()

    elif not end and node is not start:
        end = node
        end.make_end()

    elif not start and node is not end:
        start = node
        start.make_start()

    elif node is not end and node is not start:
        
        if barriers_are_vis:
            node.make_vis_barrier()
        
        else:
            node.make_invis_barrier()
    
    return (grid, start, end)


def handle_right_click(env, node):
    grid, start, end = env
    node.reset()
    
    if node is start:
        start = None
    
    elif node is end:
        end = None

    return (grid, start, end)


def handle_search_keys(
        event, env, invis_barriers, dimns, window
):
    grid, start, end = env
    
    [node.update_neighbors(grid) for row in grid for node in row]
    
    start.make_start()
    end.make_end()
    prior_env = None

    # do A* without traversal to goal
    if event.key == pg.K_1:
        prior_env = duplicate_grid(grid, dimns)
        invis_barriers = [
            node for row in grid for node in row 
            if node.is_invis_barrier()
        ]
        
        perform_a_star(
            lambda: draw(window, grid, dimns), 
            env, invis_barriers, False
        )

    # do A* with traversal to goal
    elif event.key == pg.K_2:
        prior_env = duplicate_grid(grid, dimns)
        invis_barriers = [
            node for row in grid for node in row 
            if node.is_invis_barrier()
        ]

        perform_a_star(
            lambda: draw(window, grid, dimns), 
            env, invis_barriers, True
        )

    # do LPA*
    elif event.key == pg.K_3:
        prior_env = duplicate_grid(grid, dimns)
        invis_barriers = [
            node for row in grid for node in row 
            if node.is_invis_barrier()
        ]

        perform_lpa_star(
            lambda: draw(window, grid, dimns),
            env, invis_barriers
        )

    # do D* Lite
    elif event.key == pg.K_4:
        prior_env = duplicate_grid(grid, dimns)

        perform_d_star_lite(
            lambda: draw(window, grid, dimns), env
        )
    
    return prior_env, invis_barriers


def handle_prep_keys(
        event, env, prior_env, barriers_are_vis, dimns
):
    grid, start, end = env

    # clear grid
    if event.key == pg.K_c:
        start = None
        end = None
        grid = make_grid(dimns)

        # "safe clear" grid - remove everything except barriers, 
        # start, and end
    elif event.key == pg.K_s:
        l = lambda n: n.is_open() or n.is_closed() or n.is_path()
        [node.reset() for row in grid for node in row if l(node)]

    # toggle barrier type (vis or invis)
    elif event.key == pg.K_z:
        barriers_are_vis = not barriers_are_vis

    # generate barriers
    elif event.key == pg.K_g:
        start = None
        end = None
        grid = make_grid(dimns)
        generate_barriers(grid, barriers_are_vis)

    # generate a mix of invis and vis barriers
    elif event.key == pg.K_m:
        start = None
        end = None
        grid = make_grid(dimns)
        generate_barriers_mixed(grid)

    # restore grid to original version preceding latest algorithm 
    # visualization
    elif event.key == pg.K_b:
        prev_grid, prev_start, prev_end = prior_env
        grid = prev_grid
        start = prev_start
        end = prev_end
    
    env = (grid, start, end)
    
    return (env, barriers_are_vis)


def handle_interaction(
        event, env, prior_env, barriers_are_vis, 
        invis_barriers, dimns, window
):
    grid, start, end = env

    pos = pg.mouse.get_pos()
    row, col = get_clicked_pos(pos, dimns)
    node = grid[row][col]

    if pg.mouse.get_pressed()[0]: # left click
        env = handle_left_click(env, node, barriers_are_vis)
    
    elif pg.mouse.get_pressed()[2]:  # Right click
        env = handle_right_click(env, node)

    elif event.type == pg.KEYDOWN:

        if start and end:
            results = handle_search_keys(
                event, env, invis_barriers, dimns, window
            )
            triple, barriers = results
            prior_env = triple if triple else prior_env
            invis_barriers = barriers
        
        results = handle_prep_keys(
            event, env, prior_env, barriers_are_vis, dimns
        )
        env, barriers_are_vis = results
    
    return (env, prior_env, barriers_are_vis, invis_barriers)


def main(window, dimns):
    grid = make_grid(dimns)
    start = None
    end = None
    run = True
    barriers_are_vis = True
    invis_barriers = []
    env = (grid, start, end)
    prior_env = env

    while run:
        grid, _, _ = env
        draw(window, grid, dimns)

        for event in pg.event.get():

            if event.type == pg.QUIT:
                run = False

            results = handle_interaction(
                event, env, prior_env, barriers_are_vis, 
                invis_barriers, dimns, window
            )

            env, prior_env, barriers_are_vis, invis_barriers = results

    pg.quit()


main(WINDOW, DIMNS)