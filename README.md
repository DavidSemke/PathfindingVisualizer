# PathfindingVisualizer

Visualizes search algorithms, including A*, LPA*, and D* Lite.

To start the visualizer, run the visualizer.py file. Controls are included in the file, but they are also displayed here.

## Mouse Commands

LEFT MOUSE CLICK
- The first click places the start node.
- The second click places the goal node.
- Clicks after these nodes are placed will create barriers of the type determined by whatever type has been selected using the 'z' key toggle.
- Note that algorithms cannot start until a start and end node have been selected.

RIGHT MOUSE CLICK
- Clears a node of any unique status (becomes an empty node).

## Keyboard Commands

KEY 'z'
- Toggle between vis (visible) and invis (invisible) barriers.

KEY 'g'
- This will generate a graph with barriers. 
- The number of barriers corresponds to the set obstacle density of the program, which is by default 25%. 
- Toggle the default obstacle density by changing the BARRIER_RAND_CONST constant in grid.py. 
- In the program, you can toggle the type of barrier to be generated (visible or invisible) by pressing the key 'z'.

KEY 'm'
- This will generate a graph with barriers like when pressing key 'g', except in this case a mixture of invis and vis barriers will be generated. 
- When a node is selected to be a barrier during generation, there is a 50-50 chance of an invisible or visible barrier being generated.

KEY 'c'
- Clear the grid (all nodes become empty nodes).

KEY 's'
- Clear the grid of all nodes that are paths, open, or closed. 
- Node statuses are described in node.py (as well as below).

KEY 'b'
- Restore the graph to its state before the last search algorithm was initiated.

KEY '1'
- Initiate A* without travel.

KEY '2'
- Initiate A* with travel.

KEY '3'
- Initiate LPA*.

KEY '4'
- Initiate D* Lite.


## The following are descriptions of all node statuses. Corresponding colors are given further down.

CLOSED NODE
- In the case of A*, a closed node is a node that has been visited due to its f-value being the smallest.
- For LPA* and D* Lite, however, a closed node is a node that is locally consistent, where a node is locally consistent when its rhs-value is equal to its g-value.

OPEN NODE
- In all algorithms, the open node is a node that is in the open set, which is a priority queue.

PATH NODE
- A path node signifies that the node is part of the path to the goal node (may not be part of the shortest path, in the case of locating and moving to the goal).

INVISIBLE BARRIER
- A node which is not able to be traversed to, but the agent is not yet aware of this.

VISIBLE BARRIER
- A node which is not able to be traversed to, and the agent is aware of this.

ORIGINAL START NODE
- This node is simply the starting point for an algorithm which seeks to move the agent from the start node to the goal node.

EMPTY NODE
- A node which has no unique status.

START NODE
- The node at which the agent starts calculating the shortest path.
- This node changes when the algorithm demands the agent to traverse from the start node to the end node.

END/GOAL NODE
- The node that the agent is attempting to find a shortest path to, and perhaps traverse to.

## Node Status Colors

RED
- signifies a closed node
GREEN 
- signifies an open node
BLUE
- signifies an invisible barrier
YELLOW
- signifies the original start node if the start node is not static
WHITE
- signifies regular node (default color for a node)
BLACK
- signifies a visible barrier
PURPLE
- signifies a node which is part of the path taken to reach the goal
ORANGE
- signifies the start node
GREY
- for drawing the lines in the graph
TURQUOISE
- signifies the goal node