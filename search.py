import min_heap

__all__ = ['perform_lpa_star', 'perform d_star_lite', 'perform_a_star']

# universal edge cost (edges costing infinity are simply ignored)
NODE_TO_NODE_DISTANCE = 1

'''The following are global counters used in each pathfinding algorithm.'''
# PERCOLATES - Total heap percolates
# EXPANSIONS - Total expansions of nodes (expansion = visiting a node AND updating its g_value)
# ACCESSES - Total accesses of nodes (updating and reading node attributes, including g-values, f_values,
# and rhs-values)
# COUNT - Ensures each entry into the min heap is unique, so that a tie breaker between priorities leads to the oldest
# entry being prioritized

# Each of these variables reset to 0 at the start of a new algorithm
PERCOLATES = 0
EXPANSIONS = 0
ACCESSES = 0
COUNT = 0

# Define the heuristic for a node. Uses Manhattan distance.
def heuristic(p1, p2):  # note that the order of points as parameters does not matter
    x1, y1 = p1
    x2, y2 = p2
    return abs(x1 - x2) + abs(y1 - y2)


# to reconstruct path, find the predecessor which minimizes g_value + edge cost, and
# repeat until the predecessor is the start node
# edge cost is fixed at NODE_TO_NODE_DISTANCE, so can just compare g-values
def reconstruct_path(start, end, g_value_dict, draw_func):
    current = end
    next_node = None
    while current is not start:
        min_distance = float("inf")
        for n in current.neighbors:
            possible_new_min_distance = g_value_dict[n]
            if possible_new_min_distance < min_distance:
                min_distance = possible_new_min_distance
                next_node = n
        current = next_node
        if not current.is_invisible_barrier():
            current.make_path()
        draw_func()

# LPA*______________________________________________________________________________________________________________

# nodes are prioritized by their keys, where smaller key values are of greater priority
def calculate_keys_lpa(node, end, g_value_dict, rhs_value_dict):
    global ACCESSES

    g_value = g_value_dict[node]
    rhs_value = rhs_value_dict[node]
    # f_value correspondence
    key1 = min(g_value, rhs_value) + heuristic(end.get_pos(), node.get_pos())
    # g_value correspondence
    key2 = min(g_value, rhs_value)

    ACCESSES += 2

    return key1, key2


def update_node_lpa(node_to_update, end, g_value_dict, rhs_value_dict, open_set):
    global COUNT
    global PERCOLATES
    global ACCESSES

    index_of_item = None
    index = 0
    for item in open_set:
        if node_to_update is item[2]:
            index_of_item = index
            break
        index += 1

    locally_inconsistent = g_value_dict[node_to_update] != rhs_value_dict[node_to_update]

    ACCESSES += 2

    if locally_inconsistent and index_of_item is not None:
        k1, k2 = calculate_keys_lpa(node_to_update, end, g_value_dict, rhs_value_dict)
        PERCOLATES += min_heap.heapremove(open_set, index_of_item)
        PERCOLATES += min_heap.heappush(open_set, ((k1, k2), COUNT, node_to_update))
        COUNT += 1

    elif locally_inconsistent and index_of_item is None:
        k1, k2 = calculate_keys_lpa(node_to_update, end, g_value_dict, rhs_value_dict)
        PERCOLATES += min_heap.heappush(open_set, ((k1, k2), COUNT, node_to_update))
        COUNT += 1

        if not node_to_update.is_invisible_barrier():
            node_to_update.make_open()

    elif not locally_inconsistent and index_of_item is not None:
        PERCOLATES += min_heap.heapremove(open_set, index_of_item)


#  draw_func is necessary for updating grid
#  start and end are start and goal nodes
def lpa_star_compute_shortest_path(draw_func, g_value_dict, rhs_value_dict, open_set, start, end):
    global PERCOLATES
    global EXPANSIONS
    global ACCESSES

    # get top item in priority queue
    while ((open_set and open_set[0][0] < calculate_keys_lpa(end, end, g_value_dict, rhs_value_dict))
           or rhs_value_dict[end] != g_value_dict[end]):

        current = open_set[0][2]

        ACCESSES += 4  # accesses in while loop and in following if statement
        if g_value_dict[current] > rhs_value_dict[current]:
            g_value_dict[current] = rhs_value_dict[current]
            PERCOLATES += min_heap.heappop(open_set)[1]
            ACCESSES += 2
            for node in current.neighbors:
                if node is not start:
                    rhs_value_dict[node] = min(rhs_value_dict[node], g_value_dict[current] + NODE_TO_NODE_DISTANCE)
                    ACCESSES += 3

                update_node_lpa(node, end, g_value_dict, rhs_value_dict, open_set)

            if current is not start and not current.is_invisible_barrier():
                current.make_closed()

        else:
            old_g_value = g_value_dict[current]
            g_value_dict[current] = float('inf')
            ACCESSES += 2

            neighbors_and_current = current.neighbors + [current]
            for node in neighbors_and_current:
                ACCESSES += 1
                if rhs_value_dict[node] == old_g_value + NODE_TO_NODE_DISTANCE or node is current:
                    if node is not start:
                        min_distance = float("inf")
                        for n in node.neighbors:
                            ACCESSES += 1
                            possible_new_min_distance = g_value_dict[n]
                            if possible_new_min_distance < min_distance:
                                min_distance = possible_new_min_distance

                        rhs_value_dict[node] = min_distance + NODE_TO_NODE_DISTANCE
                        ACCESSES += 1

                update_node_lpa(node, end, g_value_dict, rhs_value_dict, open_set)

        EXPANSIONS += 1

        draw_func()

    ACCESSES += 1
    if g_value_dict[end] == float('inf'):
        return False
    else:
        reconstruct_path(start, end, g_value_dict, draw_func)
        end.make_end()
        start.make_start()
        return True

# main method for LPA*
def perform_lpa_star(draw_func, grid, start, end, invisible_barriers):
    global PERCOLATES
    global EXPANSIONS
    global ACCESSES
    global COUNT
    PERCOLATES = 0
    EXPANSIONS = 0
    ACCESSES = 0
    COUNT = 0

    open_set_heap = []
    # create dictionary of g_values set to infinity
    g_value_dict = {node: float("inf") for row in grid for node in row}
    # create dictionary of rhs_values set to infinity
    rhs_value_dict = {node: float("inf") for row in grid for node in row}

    # set rhs_value of start node to 0
    rhs_value_dict[start] = 0
    ACCESSES += 1

    open_set_heap.append(((heuristic(end.get_pos(), start.get_pos()), 0), COUNT, start))
    COUNT += 1

    invisible_barriers_index = 0

    lpa_star_compute_shortest_path(draw_func, g_value_dict, rhs_value_dict, open_set_heap, start, end)

    while invisible_barriers_index < len(invisible_barriers):
        b = invisible_barriers[invisible_barriers_index]
        b.make_visible_barrier()

        invisible_barriers_index += 1

        index = 0
        for item in open_set_heap:
            if b is item[2]:
                PERCOLATES += min_heap.heapremove(open_set_heap, index)
                break
            index += 1
        b.update_neighbors(grid)

        for n_neighbor in b.neighbors:
            n_neighbor.update_neighbors(grid)
            ACCESSES += 2
            if rhs_value_dict[n_neighbor] == g_value_dict[b] + NODE_TO_NODE_DISTANCE:
                if n_neighbor is not start:
                    min_distance = float("inf")
                    for neighbor_of_n_neighbor in n_neighbor.neighbors:
                        ACCESSES += 1
                        possible_new_min_distance = g_value_dict[neighbor_of_n_neighbor]
                        if possible_new_min_distance < min_distance:
                            min_distance = possible_new_min_distance
                    rhs_value_dict[n_neighbor] = min_distance + NODE_TO_NODE_DISTANCE
                    ACCESSES += 1

            update_node_lpa(n_neighbor, end, g_value_dict, rhs_value_dict, open_set_heap)

        # clear previous activity
        for row in grid:
            for node in row:
                if node.is_path() or node.is_closed() or node.is_open():
                    node.reset()

        draw_func()

        lpa_star_compute_shortest_path(draw_func, g_value_dict, rhs_value_dict, open_set_heap, start, end)

    else:
        print("LPA* completed trials.")

    print("Heap percolates: " + str(PERCOLATES))
    print("Node expansions: " + str(EXPANSIONS))
    print("Node accesses: " + str(ACCESSES))

# _________________________________________________________________________________________________________________




# D* Lite__________________________________________________________________________________________________________

# global k_variable
K_VARIABLE = 0

# k_variable is the variable described in the D* Lite algorithm which adjusts the heuristic key values
# in order to avoid reordering heap
def calculate_keys_d_lite(node, start, g_value_dict, rhs_value_dict):
    global ACCESSES

    g_value = g_value_dict[node]
    rhs_value = rhs_value_dict[node]

    key1 = min(g_value, rhs_value) + heuristic(start.get_pos(), node.get_pos()) + K_VARIABLE
    key2 = min(g_value_dict[node], rhs_value_dict[node])

    ACCESSES += 2

    return key1, key2


def update_node_d_lite(node_to_update, start, g_value_dict, rhs_value_dict, open_set):
    global PERCOLATES
    global ACCESSES
    global COUNT

    index_of_item = None
    index = 0
    for item in open_set:
        if node_to_update is item[2]:
            index_of_item = index
            break
        index += 1

    locally_inconsistent = g_value_dict[node_to_update] != rhs_value_dict[node_to_update]

    ACCESSES += 2

    if locally_inconsistent and index_of_item is not None:
        k1, k2 = calculate_keys_d_lite(node_to_update, start, g_value_dict, rhs_value_dict)
        PERCOLATES += min_heap.heapremove(open_set, index_of_item)
        PERCOLATES += min_heap.heappush(open_set, ((k1, k2), COUNT, node_to_update))
        COUNT += 1

    elif locally_inconsistent and index_of_item is None:
        k1, k2 = calculate_keys_d_lite(node_to_update, start, g_value_dict, rhs_value_dict)
        PERCOLATES += min_heap.heappush(open_set, ((k1, k2), COUNT, node_to_update))
        COUNT += 1

        if not node_to_update.is_invisible_barrier() and not node_to_update.is_path():
            node_to_update.make_open()

    elif not locally_inconsistent and index_of_item is not None:
        PERCOLATES += min_heap.heapremove(open_set, index_of_item)


def d_star_lite_compute_shortest_path(draw_func, g_value_dict, rhs_value_dict, open_set, start, end):
    global PERCOLATES
    global EXPANSIONS
    global ACCESSES
    global COUNT

    while ((open_set and open_set[0][0] < calculate_keys_d_lite(start, start, g_value_dict, rhs_value_dict))
           or rhs_value_dict[start] != g_value_dict[start]):

        current = open_set[0][2]

        k_old = open_set[0][0]
        k_new = calculate_keys_d_lite(current, start, g_value_dict, rhs_value_dict)

        ACCESSES += 2  # referring to those in while loop check

        if k_old < k_new:
            PERCOLATES += min_heap.heapreplace(open_set, (k_new, COUNT, current))[1]
            COUNT += 1
        else:
            ACCESSES += 2
            if g_value_dict[current] > rhs_value_dict[current]:
                ACCESSES += 2
                g_value_dict[current] = rhs_value_dict[current]
                PERCOLATES += min_heap.heappop(open_set)[1]
                for node in current.neighbors:
                    if node is not end:
                        rhs_value_dict[node] = min(rhs_value_dict[node], g_value_dict[current] + NODE_TO_NODE_DISTANCE)
                        ACCESSES += 3

                    update_node_d_lite(node, start, g_value_dict, rhs_value_dict, open_set)

                if current is not end and not current.is_invisible_barrier() and not current.is_path():
                    current.make_closed()

            else:
                old_g_value = g_value_dict[current]
                g_value_dict[current] = float('inf')
                ACCESSES += 2

                neighbors_and_current = current.neighbors + [current]
                for node in neighbors_and_current:
                    ACCESSES += 1
                    if rhs_value_dict[node] == old_g_value + NODE_TO_NODE_DISTANCE:
                        if node is not end:
                            min_distance = float("inf")

                            for n in node.neighbors:
                                ACCESSES += 1
                                possible_new_min_distance = g_value_dict[n]
                                if possible_new_min_distance < min_distance:
                                    min_distance = possible_new_min_distance

                            ACCESSES += 1
                            rhs_value_dict[node] = min_distance + NODE_TO_NODE_DISTANCE

                    update_node_d_lite(node, start, g_value_dict, rhs_value_dict, open_set)

        EXPANSIONS += 1

        draw_func()

    ACCESSES += 1
    if g_value_dict[start] == float('inf'):
        return False
    else:
        return True


def perform_d_star_lite(draw_func, grid, start, end):
    global PERCOLATES
    global EXPANSIONS
    global ACCESSES
    global COUNT
    global K_VARIABLE
    PERCOLATES = 0
    EXPANSIONS = 0
    ACCESSES = 0
    COUNT = 0
    K_VARIABLE = 0

    # priority queue as a heap
    open_set_heap = []
    # create dictionary of g_values set to infinity
    g_value_dict = {node: float("inf") for row in grid for node in row}
    # create dictionary of rhs_values set to infinity
    rhs_value_dict = {node: float("inf") for row in grid for node in row}
    # set rhs_value of goal node to 0
    rhs_value_dict[end] = 0
    ACCESSES += 1

    open_set_heap.append(((heuristic(end.get_pos(), start.get_pos()), 0), COUNT, end))
    COUNT += 1

    origin = start
    true_origin = origin  # origin will be updated, so this is to remember first origin position

    # make all invisible barriers within range of agent's vision at start position visible
    for n in true_origin.neighbors:
        if n.is_invisible_barrier():
            n.make_visible_barrier()
            for n_neighbor in n.neighbors:
                n_neighbor.update_neighbors(grid)

    draw_func()

    # compute_shortest_path() returns True if a path, ultimately the shortest, to the goal is found
    if d_star_lite_compute_shortest_path(draw_func, g_value_dict, rhs_value_dict, open_set_heap, start, end):
        while start is not end:

            next_start_node = None
            min_distance = float("inf")
            for n in start.neighbors:
                ACCESSES += 1
                possible_new_min_distance = g_value_dict[n]
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
                if n.is_invisible_barrier():
                    n.make_visible_barrier()
                    nodes_changed.append(n)
                    # remove from heap if present
                    index = 0
                    for item in open_set_heap:
                        if n is item[2]:
                            PERCOLATES += min_heap.heapremove(open_set_heap, index)
                            break
                        index += 1

            if nodes_changed:
                K_VARIABLE += heuristic(origin.get_pos(), start.get_pos())
                origin = start

                # note that some code has been omitted here, as the code would not apply to an environment with
                # solely traversable and not traversable edges (edge is either some constant or infinity)
                for n in nodes_changed:
                    n.update_neighbors(grid)
                    for n_neighbor in n.neighbors:
                        n_neighbor.update_neighbors(grid)
                        ACCESSES += 2
                        if rhs_value_dict[n_neighbor] == g_value_dict[n] + NODE_TO_NODE_DISTANCE:
                            if n_neighbor is not end:
                                min_distance = float("inf")
                                for neighbor_of_n_neighbor in n_neighbor.neighbors:
                                    ACCESSES += 1
                                    possible_new_min_distance = g_value_dict[neighbor_of_n_neighbor]
                                    if possible_new_min_distance < min_distance:
                                        min_distance = possible_new_min_distance

                                ACCESSES += 1
                                rhs_value_dict[n_neighbor] = min_distance + NODE_TO_NODE_DISTANCE

                        update_node_d_lite(n_neighbor, start, g_value_dict, rhs_value_dict, open_set_heap)

                    for row in grid:
                        for node in row:
                            if node.is_closed() or node.is_open():
                                node.reset()

                if not d_star_lite_compute_shortest_path(draw_func, g_value_dict, rhs_value_dict, open_set_heap,
                                                         start, end):
                    break  # leave while loop if a path to start from goal does not exist

        if start is end:
            true_origin.make_original_start()
            print("Journey completed via D* Lite.")
        else:
            print("Journey unsuccessful - D* Lite failed to find path to goal.")
    else:
        print("Journey unsuccessful - D* Lite failed to find path to goal.")

    print("Heap percolates: " + str(PERCOLATES))
    print("Node expansions: " + str(EXPANSIONS))
    print("Node accesses: " + str(ACCESSES))


# _________________________________________________________________________________________________________________




# A*_______________________________________________________________________________________________________________


def a_star_compute_shortest_path(draw_func, start, end, g_value_dict, f_value_dict, open_set):
    global PERCOLATES
    global EXPANSIONS
    global ACCESSES
    global COUNT

    while open_set:  # heap is not empty
        item, percolates = min_heap.heappop(open_set)
        current = item[2]
        PERCOLATES += percolates

        if current is end:
            return True

        for neighbor in current.neighbors:
            temp_g_score = g_value_dict[current] + NODE_TO_NODE_DISTANCE
            ACCESSES += 2
            if temp_g_score < g_value_dict[neighbor]:
                ACCESSES += 2
                g_value_dict[neighbor] = temp_g_score
                f_value_dict[neighbor] = temp_g_score + heuristic(neighbor.get_pos(), end.get_pos())

                # check for presence of neighbor in heap
                neighbor_not_in_heap = True
                for item in open_set:
                    if neighbor is item[2]:
                        neighbor_not_in_heap = False

                if neighbor_not_in_heap:
                    ACCESSES += 1
                    PERCOLATES += min_heap.heappush(open_set, (f_value_dict[neighbor], COUNT, neighbor))
                    COUNT += 1
                    if not neighbor.is_invisible_barrier() and not neighbor.is_path() and not neighbor.is_start():
                        neighbor.make_open()

        draw_func()

        if current is not start and not current.is_invisible_barrier() and not current.is_path():
            current.make_closed()

        EXPANSIONS += 1

    return False


def perform_a_star(draw_func, grid, start, end, invisible_barriers, travel):
    global PERCOLATES
    global EXPANSIONS
    global ACCESSES
    global COUNT
    PERCOLATES = 0
    EXPANSIONS = 0
    ACCESSES = 0
    COUNT = 0

    g_value_dict = {node: float("inf") for row in grid for node in row}
    f_value_dict = {node: float("inf") for row in grid for node in row}

    # if travel true, the agent will attempt to move from start to goal position; else, the shortest path is computed
    if travel:
        a_star_with_travel(draw_func, start, end, grid, g_value_dict, f_value_dict)
    else:
        a_star_without_travel(draw_func, start, end, grid, g_value_dict, f_value_dict,
                              invisible_barriers)

    print("Heap percolates: " + str(PERCOLATES))
    print("Node expansions: " + str(EXPANSIONS))
    print("Node accesses: " + str(ACCESSES))


# compute shortest path with A*, and traverse to goal position
def a_star_with_travel(draw_func, start, end, grid, g_value_dict, f_value_dict):
    global PERCOLATES
    global ACCESSES
    global COUNT

    origin = start

    ACCESSES += 2
    g_value_dict[end] = 0
    f_value_dict[end] = heuristic(start.get_pos(), end.get_pos())

    open_set_heap = [(0, COUNT, end)]
    COUNT += 1

    # make all invisible barriers within range of agent's vision at start position visible
    for n in origin.neighbors:
        if n.is_invisible_barrier():
            n.make_visible_barrier()
            for n_neighbor in n.neighbors:
                n_neighbor.update_neighbors(grid)

    if a_star_compute_shortest_path(draw_func, end, start, g_value_dict, f_value_dict, open_set_heap):
        while start is not end:

            next_start_node = None
            min_distance = float("inf")
            for n in start.neighbors:
                ACCESSES += 1
                possible_new_min_distance = g_value_dict[n]
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
                if n.is_invisible_barrier():
                    n.make_visible_barrier()
                    nodes_changed.append(n)
                    # remove from heap if present
                    index = 0
                    for item in open_set_heap:
                        if n is item[2]:
                            PERCOLATES += min_heap.heapremove(open_set_heap, index)
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
                for row in grid:
                    for node in row:
                        if (node.is_closed() or node.is_open() or node.is_path() or node.is_start() or node.is_end()
                                or node.is_invisible_barrier()):
                            if node.is_open() or node.is_closed():
                                node.reset()
                            ACCESSES += 2
                            g_value_dict[node] = float("inf")
                            f_value_dict[node] = float("inf")
                COUNT = 0
                open_set_heap = [(0, COUNT, end)]
                COUNT += 1
                g_value_dict[end] = 0
                f_value_dict[end] = heuristic(start.get_pos(), end.get_pos())
                ACCESSES += 2

                if not a_star_compute_shortest_path(draw_func, end, start, g_value_dict, f_value_dict, open_set_heap):
                    break  # leave while loop if a path to start from goal does not exist

        if start is end:
            origin.make_original_start()
            print("Journey completed via A* (with travel).")
        else:
            print("Journey unsuccessful - A* (with travel) failed to find path to goal.")

    else:
        print("A* (with travel) failed to find path to goal.")


# compute shortest path with A*, and adapt path to changes if they occur (without traversing to goal)
def a_star_without_travel(draw_func, start, end, grid, g_value_dict, f_value_dict, invisible_barriers):
    global PERCOLATES
    global ACCESSES
    global COUNT

    ACCESSES += 2
    g_value_dict[start] = 0
    f_value_dict[start] = heuristic(start.get_pos(), end.get_pos())

    open_set_heap = [(0, COUNT, start)]
    COUNT += 1

    invisible_barriers_index = 0

    if a_star_compute_shortest_path(draw_func, start, end, g_value_dict, f_value_dict, open_set_heap):
        reconstruct_path(start, end, g_value_dict, draw_func)
        end.make_end()
        start.make_start()

        draw_func()

    while invisible_barriers_index < len(invisible_barriers):
        b = invisible_barriers[invisible_barriers_index]
        b.make_visible_barrier()

        invisible_barriers_index += 1

        index = 0
        for item in open_set_heap:
            if b is item[2]:
                PERCOLATES += min_heap.heapremove(open_set_heap, index)
                break
            index += 1
        b.update_neighbors(grid)

        for b_neighbor in b.neighbors:
            b_neighbor.update_neighbors(grid)

        # clear previous activity
        for row in grid:
            for node in row:
                if (node.is_path() or node.is_closed() or node.is_open() or node.is_start() or node.is_end()
                        or node.is_invisible_barrier()):
                    if node.is_closed() or node.is_open() or node.is_path():
                        node.reset()
                    ACCESSES += 2
                    g_value_dict[node] = float("inf")
                    f_value_dict[node] = float("inf")

        COUNT = 0
        open_set_heap = [(0, COUNT, start)]
        COUNT += 1
        g_value_dict[start] = 0
        f_value_dict[start] = heuristic(start.get_pos(), end.get_pos())

        ACCESSES += 2

        if a_star_compute_shortest_path(draw_func, start, end, g_value_dict, f_value_dict, open_set_heap):
            reconstruct_path(start, end, g_value_dict, draw_func)
            end.make_end()
            start.make_start()

            draw_func()

    else:
        print("A* (without travel) completed trials.")
