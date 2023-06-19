import min_heap
from search_tools import *

def a_star_compute_shortest_path(
        draw_func, env, g_dict, f_dict, open_set
):
    global PERCOLATES
    global EXPANSIONS
    global ACCESSES
    global COUNT

    _, start, end = env

    while open_set:  # heap is not empty
        item, percolates = min_heap.heappop(open_set)
        current = item[2]
        PERCOLATES += percolates

        if current is end: return True

        for neighbor in current.neighbors:
            temp_g_score = g_dict[current] + NODE_TO_NODE_DIST
            ACCESSES += 2
            
            if temp_g_score < g_dict[neighbor]:
                ACCESSES += 2
                g_dict[neighbor] = temp_g_score
                f_dict[neighbor] = temp_g_score + heuristic(neighbor.get_pos(), end.get_pos())

                # check for presence of neighbor in heap
                neighbor_not_in_heap = True
                
                for item in open_set:
                    if neighbor is item[2]:
                        neighbor_not_in_heap = False

                if neighbor_not_in_heap:
                    ACCESSES += 1
                    PERCOLATES += min_heap.heappush(open_set, (f_dict[neighbor], COUNT, neighbor))
                    COUNT += 1
                    if not neighbor.is_invis_barrier() and not neighbor.is_path() and not neighbor.is_start():
                        neighbor.make_open()

        draw_func()

        if not (
            current is start 
            or current.is_invis_barrier() 
            or current.is_path()
        ):
            current.make_closed()

        EXPANSIONS += 1

    return False


def perform_a_star(draw_func, env, invis_barriers, travel):
    global PERCOLATES
    global EXPANSIONS
    global ACCESSES
    global COUNT
    PERCOLATES = 0
    EXPANSIONS = 0
    ACCESSES = 0
    COUNT = 0

    grid, _, _ = env

    g_dict = {node: float("inf") for row in grid for node in row}
    f_dict = {node: float("inf") for row in grid for node in row}

    # if travel true, the agent will attempt to move from start to 
    # goal position; else, the shortest path is computed
    if travel:
        a_star_with_travel(draw_func, env, g_dict, f_dict)
    
    else:
        a_star_without_travel(
            draw_func, env, g_dict, f_dict, invis_barriers
        )

    print("Heap percolates: " + str(PERCOLATES))
    print("Node expansions: " + str(EXPANSIONS))
    print("Node accesses: " + str(ACCESSES))


# compute shortest path with A*, and traverse to goal position
def a_star_with_travel(draw_func, env, g_dict, f_dict):
    global PERCOLATES
    global ACCESSES
    global COUNT

    grid, start, end = env
    # reorder to start at end node when computing shortest path
    env = (grid, end, start)

    origin = start

    ACCESSES += 2
    g_dict[end] = 0
    f_dict[end] = heuristic(start.get_pos(), end.get_pos())

    open_set_heap = [(0, COUNT, end)]
    COUNT += 1

    # make all invis barriers within range of agent's vision at start position vis
    for n in origin.neighbors:
        
        if n.is_invis_barrier():
            n.make_vis_barrier()
            
            for n_neighbor in n.neighbors:
                n_neighbor.update_neighbors(grid)

    if a_star_compute_shortest_path(
        draw_func, env, g_dict, f_dict, open_set_heap
    ):
        
        while start is not end:
            next_start_node = None
            min_distance = float("inf")
            
            for n in start.neighbors:
                ACCESSES += 1
                possible_new_min_distance = g_dict[n]
                
                if possible_new_min_distance < min_distance:
                    min_distance = possible_new_min_distance
                    next_start_node = n

            start.make_path()
            start = next_start_node
            start.make_start()

            draw_func()

            # the next step simulates scanning for changes in edge costs
            # changes to graph can occur one node away from the start node in any direction
            nodes_changed = []
            
            for n in start.neighbors:
                
                if n.is_invis_barrier():
                    n.make_vis_barrier()
                    nodes_changed.append(n)
                    # remove from heap if present
                    index = 0
                    
                    for item in open_set_heap:
                        
                        if n is item[2]:
                            PERCOLATES += min_heap.heapremove(
                                open_set_heap, index
                            )
                            break
                        
                        index += 1

            if nodes_changed:
                
                # note that some code has been omitted here, as the code would not apply to an environment with
                # solely traversable and not traversable edges (edge is either some constant or infinity)
                for n in nodes_changed:
                    n.update_neighbors(grid)
                    
                    for n_neighbor in n.neighbors:
                        n_neighbor.update_neighbors(grid)

                # reset everything and start search fresh
                for node in [node for row in grid for node in row]:

                    if (
                        node.is_closed() 
                        or node.is_open() 
                        or node.is_path() 
                        or node.is_start() 
                        or node.is_end()
                        or node.is_invis_barrier()
                    ):
                        
                        if node.is_open() or node.is_closed():
                            node.reset()
                        
                        ACCESSES += 2
                        g_dict[node] = float("inf")
                        f_dict[node] = float("inf")

                COUNT = 0
                open_set_heap = [(0, COUNT, end)]
                COUNT += 1
                g_dict[end] = 0
                f_dict[end] = heuristic(
                    start.get_pos(), end.get_pos()
                )
                ACCESSES += 2

                if not a_star_compute_shortest_path(
                    draw_func, env, g_dict, f_dict, open_set_heap
                ):
                    # leave while loop if a path to start from goal 
                    # does not exist
                    break  

        if start is end:
            origin.make_original_start()
            print("Journey completed via A* (with travel).")
        
        else:
            print("Journey unsuccessful - A* (with travel) failed to find path to goal.")

    else:
        print("A* (with travel) failed to find path to goal.")


# compute shortest path with A*, and adapt path to changes if they 
# occur (without traversing to goal)
def a_star_without_travel(
        draw_func, env, g_dict, f_dict, invis_barriers
):
    global PERCOLATES
    global ACCESSES
    global COUNT

    ACCESSES += 2

    grid, start, end = env

    g_dict[start] = 0
    f_dict[start] = heuristic(start.get_pos(), end.get_pos())

    open_set_heap = [(0, COUNT, start)]
    COUNT += 1

    invis_barriers_index = 0

    if a_star_compute_shortest_path(
        draw_func, env, g_dict, f_dict, open_set_heap
    ):
        reconstruct_path(start, end, g_dict, draw_func)
        end.make_end()
        start.make_start()

        draw_func()

    while invis_barriers_index < len(invis_barriers):
        b = invis_barriers[invis_barriers_index]
        b.make_vis_barrier()
        invis_barriers_index += 1
        index = 0
        
        for item in open_set_heap:
            
            if b is item[2]:
                PERCOLATES += min_heap.heapremove(
                    open_set_heap, index
                )
                break

            index += 1

        b.update_neighbors(grid)

        for b_neighbor in b.neighbors:
            b_neighbor.update_neighbors(grid)

        # clear previous activity
        for node in [node for row in grid for node in row]:
            if (
                node.is_path() 
                or node.is_closed() 
                or node.is_open() 
                or node.is_start() 
                or node.is_end()
                or node.is_invis_barrier()
            ):
                if (
                    node.is_closed() 
                    or node.is_open() 
                    or node.is_path()
                ):
                    node.reset()

                ACCESSES += 2
                g_dict[node] = float("inf")
                f_dict[node] = float("inf")

        COUNT = 0
        open_set_heap = [(0, COUNT, start)]
        COUNT += 1
        g_dict[start] = 0
        f_dict[start] = heuristic(start.get_pos(), end.get_pos())

        ACCESSES += 2

        if a_star_compute_shortest_path(
            draw_func, env, g_dict, f_dict, open_set_heap
        ):
            reconstruct_path(env, g_dict, draw_func)
            end.make_end()
            start.make_start()
            draw_func()

    else:
        print("A* (without travel) completed trials.")