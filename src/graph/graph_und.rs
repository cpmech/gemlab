use crate::mesh::{Mesh, PointId};
use crate::StrError;
use russell_lab::AsArray2D;
use std::collections::{HashSet, VecDeque};

/// Defines an undirected graph structure
///
/// Potential uses include renumbering a mesh to reduce the associated bandwidth.
pub struct GraphUnd {
    /// Holds the adjacency (sparse) matrix (point connections)
    ///
    /// Note: If `calc_degree` is true, each row in this matrix is sorted in ascending order of degree,
    /// followed by point ID. Otherwise, each row is sorted in ascending order of point ID only.
    /// Degree here means the number of connections of a vertex.
    ///
    /// (nnode x variable nnode)
    adjacency: Vec<Vec<PointId>>,

    /// Defines an auxiliary queue for BFS (breadth-first-search) runs
    ///
    /// Note: this queue must be cleared before each use.
    ///
    /// (nnode)
    queue: VecDeque<PointId>,

    /// Holds the auxiliary list of bool indicating that a vertex has been explored
    ///
    /// Note: this array must be cleared before each use.
    ///
    /// (nnode)
    explored: Vec<bool>,

    /// Holds all point degrees (the number of connections of a vertex)
    ///
    /// Requires: `calc_degree == true`
    ///
    /// Optional(nnode)
    degree: Vec<usize>,

    /// Holds the id of the point with the minimum degree (the number of connections of a vertex)
    ///
    /// Requires: `calc_degree == true`
    ///
    /// Optional
    p_min_degree: usize,

    /// Holds the distance from the root to each point
    ///
    /// Optional(nnode)
    distance: Vec<usize>,

    /// Holds the parent node of each point in the BFS tree
    ///
    /// Optional(nnode)
    parent: Vec<PointId>,
}

impl GraphUnd {
    /// Allocates a new instance given an adjacency set
    fn from_adjacency_set(
        adjacency_set: &Vec<HashSet<PointId>>,
        calc_degree: bool,
        check_connectivity: bool,
    ) -> Result<Self, StrError> {
        // check the connectivity of the graph (all vertices must be explored)
        let nnode = adjacency_set.len();
        let mut queue = VecDeque::with_capacity(nnode);
        let mut explored = vec![false; nnode];
        if check_connectivity {
            explored[0] = true;
            queue.push_back(0);
            while queue.len() != 0 {
                if let Some(a) = queue.pop_front() {
                    for b in &adjacency_set[a] {
                        if !explored[*b] {
                            explored[*b] = true;
                            queue.push_back(*b);
                        }
                    }
                }
            }
            for exp in &explored {
                if !exp {
                    return Err("there are hanging vertices/edges in the mesh (graph is disconnected)");
                }
            }
        }

        // allocate adjacency matrix and degree array
        let mut adjacency = Vec::new();
        let mut degree = Vec::new();
        let mut p_min_degree = 0;
        if calc_degree {
            // compute the list of degrees (number of connections of a vertex)
            degree = vec![0; nnode];
            let mut min_degree = usize::MAX;
            for (p, row) in adjacency_set.iter().enumerate() {
                degree[p] = row.len();
                if degree[p] < min_degree {
                    p_min_degree = p;
                    min_degree = degree[p];
                }
            }
            // sort each row of the adjacency matrix by the degree, and then by id
            for row_set in adjacency_set.iter() {
                let mut row: Vec<_> = row_set.iter().copied().collect();
                row.sort_unstable_by_key(|p| (degree[*p], *p));
                adjacency.push(row);
            }
        } else {
            // sort each row of the adjacency matrix by point id
            for row_set in adjacency_set.iter() {
                let mut row: Vec<_> = row_set.iter().copied().collect();
                row.sort();
                adjacency.push(row);
            }
        }

        // results
        Ok(GraphUnd {
            adjacency,
            queue,
            explored,
            degree,
            p_min_degree,
            distance: Vec::new(),
            parent: Vec::new(),
        })
    }

    /// Allocates a new instance given a list of edges
    ///
    /// # Input
    ///
    /// * `calc_degree` -- calculates the degree (the number of connections of a vertex), as required by the
    ///    Cuthill-McKee algorithm. The degree affects the sorting of rows in the adjacency matrix.
    /// * `check_connectivity` -- checks if the graph is connected
    pub fn from_edges<'a, T>(edges: &'a T, calc_degree: bool, check_connectivity: bool) -> Result<Self, StrError>
    where
        T: AsArray2D<'a, usize>,
    {
        // find number of nodes
        let (nedge, ncorner) = edges.size();
        if ncorner < 2 {
            return Err("edges must have at least two nodes");
        }
        let mut nodes = HashSet::new();
        for e in 0..nedge {
            let (a, b) = (edges.at(e, 0), edges.at(e, 1));
            nodes.insert(a);
            nodes.insert(b);
        }
        let nnode = nodes.len();

        // find the adjacency (sparse) matrix of nodes' connections
        let mut adjacency_set = vec![HashSet::new(); nnode];
        for e in 0..nedge {
            let (a, b) = (edges.at(e, 0), edges.at(e, 1));
            adjacency_set[a].insert(b);
            adjacency_set[b].insert(a);
        }

        // allocate graph
        GraphUnd::from_adjacency_set(&adjacency_set, calc_degree, check_connectivity)
    }

    /// Allocates a new instance given a mesh
    ///
    /// # Input
    ///
    /// * `calc_degree` -- calculates the degree (the number of connections of a vertex), as required by the
    ///    Cuthill-McKee algorithm. The degree affects the sorting of rows in the adjacency matrix.
    /// * `check_connectivity` -- checks if the graph is connected
    pub fn from_mesh(mesh: &Mesh, calc_degree: bool, check_connectivity: bool) -> Result<Self, StrError> {
        // find the adjacency (sparse) matrix of nodes' connections
        let nnode = mesh.points.len();
        let mut adjacency_set = vec![HashSet::new(); nnode];
        for cell in &mesh.cells {
            for a in &cell.points {
                for b in &cell.points {
                    if *b != *a {
                        adjacency_set[*a].insert(*b);
                        adjacency_set[*b].insert(*a);
                    }
                }
            }
        }

        // allocate graph
        GraphUnd::from_adjacency_set(&adjacency_set, calc_degree, check_connectivity)
    }

    /// Computes the ordering array to renumber the vertices according to the (reverse) Cuthill-McKee algorithm
    ///
    /// **Note:** All nodes must be reachable from the root; i.e., the corresponding graph must be connected.
    ///
    /// # Input
    ///
    /// * `start_point` -- (root) the first point id, which will not be renumbered. Should have a low degree.
    ///   If None, a pseudo-peripheral point is determined and used as root.
    ///
    /// # Output
    ///
    /// Returns the ordering array such that `old = ordering[new]` where `old` is the original
    /// point id and `new` is the new point id. See the function [Graph::get_old_to_new_map]
    pub fn cuthill_mckee(&mut self, start_point: Option<PointId>) -> Result<Vec<PointId>, StrError> {
        // check if the degree of vertices is available
        let nnode = self.adjacency.len();
        if self.degree.len() != nnode {
            return Err("Cuthill-McKee algorithm requires the degree of vertices (calc_degree must be set to true)");
        }

        // root point
        let root = match start_point {
            Some(p) => p,
            None => self.pseudo_peripheral(None)?,
        };

        // clear auxiliary structures
        self.queue.clear();
        for i in 0..nnode {
            self.explored[i] = false;
        }

        // allocate auxiliary structures
        let mut ordering = vec![0; nnode];

        // first label
        let mut label = 0;
        self.explored[root] = true;
        ordering[label] = root;
        label += 1;
        self.queue.push_back(root);

        // execute a breadth-first search (BFS)
        while self.queue.len() != 0 {
            if let Some(a) = self.queue.pop_front() {
                for b in &self.adjacency[a] {
                    if !self.explored[*b] {
                        self.explored[*b] = true;
                        ordering[label] = *b;
                        label += 1;
                        self.queue.push_back(*b);
                    }
                }
            }
        }

        // reverse ordering
        ordering.reverse();
        Ok(ordering)
    }

    /// Finds the shortest path between two points using the breadth-first search (BFS) algorithm
    ///
    /// Returns a list of point ids representing the shortest path between the source and destination points
    pub fn shortest_path_bfs(&mut self, source: usize, destination: usize) -> Vec<usize> {
        // clear auxiliary structures
        self.queue.clear();
        let nnode = self.adjacency.len();
        if self.distance.len() != nnode {
            self.distance = vec![0; nnode];
        }
        if self.parent.len() != nnode {
            self.parent = vec![usize::MAX; nnode];
        }
        for i in 0..nnode {
            self.explored[i] = false;
            self.distance[i] = 0;
            self.parent[i] = usize::MAX;
        }

        // run BFS
        self.explored[source] = true;
        self.queue.push_back(source);
        while self.queue.len() != 0 {
            if let Some(a) = self.queue.pop_front() {
                for b in &self.adjacency[a] {
                    if !self.explored[*b] {
                        self.explored[*b] = true;
                        self.queue.push_back(*b);
                        self.distance[*b] = self.distance[a] + 1;
                        self.parent[*b] = a;
                    }
                }
            }
        }

        // run backwards to find the path
        let mut path = Vec::with_capacity(nnode);
        let mut current = destination;
        path.insert(0, current);
        while self.parent[current] != usize::MAX {
            path.insert(0, self.parent[current]);
            current = self.parent[current];
        }
        path
    }

    /// Runs a BFS to compute the distances (levels) from every vertex to the root vertex
    ///
    /// # Input
    ///
    /// * `root` -- The root point
    ///
    /// # Output
    ///
    /// Returns the `max_distance`
    pub fn calc_distance(&mut self, root: usize) -> usize {
        // clear auxiliary structures
        self.queue.clear();
        let nnode = self.adjacency.len();
        if self.distance.len() != nnode {
            self.distance = vec![0; nnode];
        }
        for i in 0..nnode {
            self.explored[i] = false;
            self.distance[i] = 0;
        }

        // run BFS
        self.explored[root] = true;
        self.queue.push_back(root);
        while self.queue.len() != 0 {
            if let Some(a) = self.queue.pop_front() {
                for b in &self.adjacency[a] {
                    if !self.explored[*b] {
                        self.explored[*b] = true;
                        self.queue.push_back(*b);
                        self.distance[*b] = self.distance[a] + 1;
                    }
                }
            }
        }

        // calculate max distance
        let mut max_distance = 0;
        for i in 0..nnode {
            if self.distance[i] > max_distance {
                max_distance = self.distance[i];
            }
        }
        max_distance
    }

    /// Finds a pseudo-peripheral point
    ///
    /// # Input
    ///
    /// * `start_point` -- (root) the first point id, which will not be renumbered. Should have a low degree.
    ///   If None, a point with a minimum degree will be used.
    pub fn pseudo_peripheral(&mut self, start_point: Option<PointId>) -> Result<usize, StrError> {
        // check if the degree of vertices is available
        let nnode = self.adjacency.len();
        if self.degree.len() != nnode {
            return Err("pseudo_peripheral requires the degree of vertices (calc_degree must be set to true)");
        }

        // root point
        let mut root = match start_point {
            Some(p) => p,
            None => self.p_min_degree,
        };

        // first distances
        let mut max_distance = self.calc_distance(root);

        // auxiliary
        const MAX_ITERATIONS: usize = 10;
        let mut success = false;

        // perform iterations
        for _ in 0..MAX_ITERATIONS {
            // loop over all points, consider only the points with max distance
            // equal to root's max distance, and select the point with the minimum degree
            let mut next_root = None;
            let mut min_deg = usize::MAX; // min degree among potential roots
            for i in 0..nnode {
                if self.distance[i] == max_distance {
                    let deg = self.degree[i];
                    if deg < min_deg {
                        min_deg = deg;
                        next_root = Some(i);
                    }
                }
            }
            // handle next root
            match next_root {
                None => {
                    // converged with no next root
                    success = true;
                    break;
                }
                Some(r) => {
                    root = r;
                    let next_max_distance = self.calc_distance(root);
                    if next_max_distance == max_distance {
                        // converged with next root having the same distance and ≤ degree
                        success = true;
                        break;
                    }
                    max_distance = next_max_distance;
                }
            }
        }
        if !success {
            return Err("INTERNAL ERROR: iterations did not converge");
        }
        Ok(root)
    }

    /// Calculates the (half) bandwidth (with diagonal) of the adjacency matrix
    ///
    /// ```text
    /// band = max{band_i, 0 ≤ n ≤ nnode-1}
    /// band_i = max{|i - j| + 1, any j > i}
    /// ```
    pub fn calc_bandwidth(&self) -> usize {
        let nnode = self.adjacency.len();
        let mut band = 0;
        for i in 0..nnode {
            let mut band_i = 1;
            for j in &self.adjacency[i] {
                if *j > i {
                    let delta = i.abs_diff(*j) + 1;
                    if delta > band_i {
                        band_i = delta
                    }
                }
            }
            if band_i > band {
                band = band_i;
            }
        }
        band
    }

    /// Prints the non-zero pattern of the laplacian matrix
    ///
    /// ```text
    /// L = D - A
    ///
    /// L: Laplacian matrix
    /// A: Adjacency matrix
    /// D: Diagonal matrix with the degrees
    /// ```
    pub fn print_non_zero_pattern(&self) {
        let nnode = self.adjacency.len();
        let mut non_zeros_pattern = vec![vec!["."; nnode]; nnode];
        for i in 0..nnode {
            non_zeros_pattern[i][i] = "D";
            for j in &self.adjacency[i] {
                non_zeros_pattern[i][*j] = "#";
            }
        }
        let width = nnode * 2 + 1;
        println!("\n┌{:1$}┐", " ", width);
        for i in 0..nnode {
            print!("│");
            for j in 0..nnode {
                print!(" {}", non_zeros_pattern[i][j])
            }
            print!(" │\n");
        }
        println!("└{:1$}┘", " ", width);
    }

    /// Converts the ordering array to the old_to_new map
    ///
    /// ```text
    /// old = ordering[new]
    /// new = old_to_new[old]
    /// Returns the ordering array such that `ordering[new_point_id] = old_point_id`
    /// ```
    pub fn get_old_to_new_map(ordering: &[PointId]) -> Vec<PointId> {
        let n = ordering.len();
        let mut old_to_new = vec![0; n];
        for new in 0..n {
            old_to_new[ordering[new]] = new;
        }
        old_to_new
    }

    /// Renumbers a mesh using Cuthill-McKey algorithm with pseudo-peripheral starting point
    ///
    /// # Input
    ///
    /// * `check_connectivity` -- checks if the associated graph is connected
    pub fn renumber_mesh(mesh: &mut Mesh, check_connectivity: bool) -> Result<(), StrError> {
        let calc_degree = true;
        let mut graph = GraphUnd::from_mesh(&mesh, calc_degree, check_connectivity)?;
        let ordering = graph.cuthill_mckee(None)?;
        let old_to_new = GraphUnd::get_old_to_new_map(&ordering);
        mesh.renumber_points(&old_to_new)
    }

    /// Returns the number of nodes in the graph
    pub fn get_nnode(&self) -> usize {
        self.adjacency.len()
    }

    /// Returns the number of edges in the graph
    pub fn get_nedge(&self) -> usize {
        let mut count = 0;
        for row in &self.adjacency {
            count += row.len();
        }
        // Divide by 2 since each edge is counted twice (undirected graph)
        count / 2
    }

    /// Returns the degree (number of connections) of a node
    ///
    /// # Arguments
    ///
    /// * `node` - The node index
    pub fn get_degree(&self, node: usize) -> Result<usize, StrError> {
        let nnode = self.adjacency.len();
        if self.degree.len() != nnode {
            return Err("degree information is not available (calc_degree must be set to true)");
        }
        if node >= nnode {
            return Err("node index out of bounds");
        }
        Ok(self.degree[node])
    }

    /// Checks if the graph contains a given edge
    ///
    /// # Arguments
    ///
    /// * `a` - First node
    /// * `b` - Second node
    pub fn has_edge(&self, a: usize, b: usize) -> bool {
        let nnode = self.adjacency.len();
        if a >= nnode || b >= nnode {
            return false;
        }
        self.adjacency[a].contains(&b)
    }

    /// Returns a list of all edges in the graph as pairs of node indices
    pub fn get_edges(&self) -> Vec<(usize, usize)> {
        let mut edges = Vec::new();
        for (a, neighbors) in self.adjacency.iter().enumerate() {
            for &b in neighbors {
                if a < b {
                    // Only add each edge once
                    edges.push((a, b));
                }
            }
        }
        edges
    }

    /// Returns the density of the graph
    ///
    /// The density is defined as: (2 * |E|) / (|V| * (|V| - 1))
    /// where |E| is the number of edges and |V| is the number of vertices
    pub fn density(&self) -> f64 {
        let nv = self.get_nnode() as f64;
        if nv <= 1.0 {
            return 0.0;
        }
        (2.0 * self.get_nedge() as f64) / (nv * (nv - 1.0))
    }

    /// Finds all connected components in the graph
    ///
    /// Returns a vector where each element is a vector of node indices belonging to that component
    pub fn connected_components(&mut self) -> Vec<Vec<usize>> {
        let nnode = self.get_nnode();
        let mut components = Vec::new();
        let mut unvisited: HashSet<usize> = (0..nnode).collect();

        while !unvisited.is_empty() {
            let start = *unvisited.iter().next().unwrap();
            let mut component = Vec::new();
            
            // Clear auxiliary structures
            self.queue.clear();
            for i in 0..nnode {
                self.explored[i] = false;
            }

            // Run BFS from the start node
            self.explored[start] = true;
            self.queue.push_back(start);
            component.push(start);
            unvisited.remove(&start);

            while let Some(node) = self.queue.pop_front() {
                for &neighbor in &self.adjacency[node] {
                    if !self.explored[neighbor] {
                        self.explored[neighbor] = true;
                        self.queue.push_back(neighbor);
                        component.push(neighbor);
                        unvisited.remove(&neighbor);
                    }
                }
            }

            components.push(component);
        }

        components
    }

    /// Checks if the graph is connected
    ///
    /// A graph is connected if there is a path between every pair of vertices
    pub fn is_connected(&mut self) -> bool {
        self.connected_components().len() == 1
    }

    /// Finds articulation points (cut vertices) in the graph
    ///
    /// An articulation point is a vertex whose removal increases the number of connected components
    pub fn find_articulation_points(&self) -> Vec<usize> {
        let nnode = self.get_nnode();
        if nnode <= 1 {
            return Vec::new();
        }

        let mut discovery = vec![0; nnode];
        let mut low = vec![0; nnode];
        let mut parent = vec![usize::MAX; nnode];
        let mut visited = vec![false; nnode];
        let mut is_articulation = vec![false; nnode];
        let mut time = 0;

        // DFS function
        fn dfs(u: usize, 
              discovery: &mut Vec<usize>,
              low: &mut Vec<usize>,
              parent: &mut Vec<usize>,
              visited: &mut Vec<bool>,
              is_articulation: &mut Vec<bool>,
              time: &mut usize,
              graph: &GraphUnd) {
            visited[u] = true;
            *time += 1;
            discovery[u] = *time;
            low[u] = *time;
            let mut children = 0;

            for &v in &graph.adjacency[u] {
                if !visited[v] {
                    children += 1;
                    parent[v] = u;
                    dfs(v, discovery, low, parent, visited, is_articulation, time, graph);

                    low[u] = low[u].min(low[v]);

                    if parent[u] == usize::MAX && children > 1 {
                        is_articulation[u] = true;
                    }
                    if parent[u] != usize::MAX && low[v] >= discovery[u] {
                        is_articulation[u] = true;
                    }
                } else if v != parent[u] {
                    low[u] = low[u].min(discovery[v]);
                }
            }
        }

        // Run DFS from each unvisited vertex
        for i in 0..nnode {
            if !visited[i] {
                dfs(i, &mut discovery, &mut low, &mut parent, &mut visited, 
                    &mut is_articulation, &mut time, self);
            }
        }

        // Collect articulation points
        is_articulation.iter()
            .enumerate()
            .filter(|(_, &is_art)| is_art)
            .map(|(i, _)| i)
            .collect()
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::GraphUnd;
    use crate::mesh::{Block, Cell, Figure, Mesh, Point, Samples};
    use crate::shapes::GeoKind;
    use russell_lab::NumMatrix;
    use std::collections::HashSet;

    const SAVE_FIGURE: bool = false;

    #[test]
    fn from_adjacency_set_works() {
        //          .1.
        //        .' | '.
        //      .'   |   '.
        //    0'     |     3
        //     '.    |    .'
        //       '.  |  .'
        //         '.2.'
        let adjacency_set = vec![
            HashSet::from([2, 1]),    // node 0 connects to 2 and 1
            HashSet::from([0, 3, 2]), // node 1 connects to 0, 3, and 2
            HashSet::from([0, 1, 3]), // node 2 connects to 0, 1, and 3
            HashSet::from([2, 1]),    // node 3 connects to 2 and 1
        ];
        let graph = GraphUnd::from_adjacency_set(&adjacency_set, false, false).unwrap();

        assert_eq!(graph.adjacency[0], &[1, 2]); // sorted by id
        assert_eq!(graph.adjacency[1], &[0, 2, 3]); // sorted by id
        assert_eq!(graph.adjacency[2], &[0, 1, 3]); // sorted by id
        assert_eq!(graph.adjacency[3], &[1, 2]); // sorted by id

        assert_eq!(graph.explored.len(), 4);
        assert_eq!(graph.degree.len(), 0);
        assert_eq!(graph.distance.len(), 0);
        assert_eq!(graph.parent.len(), 0);

        let graph = GraphUnd::from_adjacency_set(&adjacency_set, true, false).unwrap();

        assert_eq!(graph.adjacency[0], &[1, 2]); // sorted by id
        assert_eq!(graph.adjacency[1], &[0, 3, 2]); // sorted by degree, then id
        assert_eq!(graph.adjacency[2], &[0, 3, 1]); // sorted by degree, then id
        assert_eq!(graph.adjacency[3], &[1, 2]); // sorted by id

        assert_eq!(graph.explored.len(), 4);
        assert_eq!(graph.degree, &[2, 3, 3, 2]);
        assert_eq!(graph.distance.len(), 0);
        assert_eq!(graph.parent.len(), 0);
    }

    #[test]
    fn graph_from_edges_works_4nodes() {
        // 0 ––––––––––– 3
        // │      1      │
        // │             │
        // │ 0         3 │
        // │             │
        // │      2      |
        // 1 ––––––––––– 2

        // edge:       0       1       2       3
        let edges = [[0, 1], [0, 3], [1, 2], [2, 3]];
        let graph = GraphUnd::from_edges(&edges, true, true).unwrap();

        // node:                   0  1  2  3
        assert_eq!(graph.degree, &[2, 2, 2, 2]);
        assert_eq!(graph.p_min_degree, 0);

        // for (i, row) in graph.adjacency.iter().enumerate() { println!("{}: {:?}", i, row); }

        assert_eq!(graph.adjacency[0], &[1, 3]); // sorted by id
        assert_eq!(graph.adjacency[1], &[0, 2]); // sorted by id
        assert_eq!(graph.adjacency[2], &[1, 3]); // sorted by id
        assert_eq!(graph.adjacency[3], &[0, 2]); // sorted by id
    }

    #[test]
    fn graph_from_edges_works_6nodes() {
        // 4 ––––––––––––––– 5 .
        // │        0        │  `. 6
        // │                 │    `.
        // │                 │      `.
        // │ 1             4 │        3
        // │                 │     ,'
        // │                 │   ,'
        // │    2       3    │ ,'  5
        // 1 –––––– 0 –––––– 2

        // edge:       0       1       2       3       4       5       6
        let edges = [[4, 5], [1, 4], [0, 1], [0, 2], [5, 2], [2, 3], [5, 3]];
        let graph = GraphUnd::from_edges(&edges, true, true).unwrap();

        // node:                   0  1  2  3  4  5
        assert_eq!(graph.degree, &[2, 2, 3, 2, 2, 3]);
        assert_eq!(graph.p_min_degree, 0);

        // for (i, row) in graph.adjacency.iter().enumerate() { println!("{}: {:?}", i, row); }

        assert_eq!(graph.adjacency[0], &[1, 2]); // sorted by degree, then id
        assert_eq!(graph.adjacency[1], &[0, 4]); // sorted by id
        assert_eq!(graph.adjacency[2], &[0, 3, 5]); // sorted by degree, then id
        assert_eq!(graph.adjacency[3], &[2, 5]); // sorted by id
        assert_eq!(graph.adjacency[4], &[1, 5]); // sorted by degree, then id
        assert_eq!(graph.adjacency[5], &[3, 4, 2]); // sorted by degree, then id
    }

    #[test]
    fn graph_from_mesh_works_1() {
        // lin2_graph
        let mesh = Samples::graph_8_edges();
        let graph = GraphUnd::from_mesh(&mesh, true, false).unwrap();

        //                         0  1  2  3  4  5  6  7 (point)
        assert_eq!(graph.degree, &[1, 3, 2, 2, 3, 2, 1, 2]);
        assert_eq!(graph.p_min_degree, 0);

        // for (i, row) in graph.adjacency.iter().enumerate() { println!("{}: {:?}", i, row); }

        assert_eq!(graph.adjacency[1], &[2, 5, 7]); // sorted by id
        assert_eq!(graph.adjacency[2], &[1, 4]); // sorted by id
        assert_eq!(graph.adjacency[3], &[6, 4]); // sorted by degree
        assert_eq!(graph.adjacency[4], &[0, 2, 3]); // sorted by degree, then id
        assert_eq!(graph.adjacency[5], &[7, 1]); // sorted by degree
        assert_eq!(graph.adjacency[6], &[3]);
        assert_eq!(graph.adjacency[7], &[5, 1]); // sorted by degree
    }

    #[test]
    fn graph_from_mesh_works_2() {
        // 5------4------3
        // |      |      |
        // |      |      |
        // 0------1------2
        #[rustfmt::skip]
        let mesh = Mesh {
            ndim: 2,
            points: vec![
                Point { id: 0, marker: 0, coords: vec![0.0, 0.0] },
                Point { id: 1, marker: 0, coords: vec![1.0, 0.0] },
                Point { id: 2, marker: 0, coords: vec![2.0, 0.0] },
                Point { id: 3, marker: 0, coords: vec![2.0, 1.0] },
                Point { id: 4, marker: 0, coords: vec![1.0, 1.0] },
                Point { id: 5, marker: 0, coords: vec![0.0, 1.0] },
            ],
            cells: vec![
                Cell { id: 0, attribute: 1, kind: GeoKind::Qua4, points: vec![0, 1, 4, 5] },
                Cell { id: 1, attribute: 1, kind: GeoKind::Qua4, points: vec![1, 2, 3, 4] },
            ],
        };
        let graph = GraphUnd::from_mesh(&mesh, true, false).unwrap();

        //                         0  1  2  3  4  5 (point)
        assert_eq!(graph.degree, &[3, 5, 3, 3, 5, 3]);

        // for (i, row) in graph.adjacency.iter().enumerate() { println!("{}: {:?}", i, row); }

        assert_eq!(graph.adjacency[0], &[5, 1, 4]); // sorted by degree, then id
        assert_eq!(graph.adjacency[1], &[0, 2, 3, 5, 4]); // sorted by degree, then id
        assert_eq!(graph.adjacency[2], &[3, 1, 4]); // sorted by degree, then id
        assert_eq!(graph.adjacency[3], &[2, 1, 4]); // sorted by degree, then id
        assert_eq!(graph.adjacency[4], &[0, 2, 3, 5, 1]); // sorted by degree, then id
        assert_eq!(graph.adjacency[5], &[0, 1, 4]); // sorted by degree, then id

        let npoint = mesh.points.len();
        let mut incidence = NumMatrix::<usize>::new(npoint, npoint);
        for i in 0..npoint {
            for j in &graph.adjacency[i] {
                incidence.set(i, *j, 1);
            }
        }

        assert_eq!(
            format!("{}", incidence),
            // 0 1 2 3 4 5
            "┌             ┐\n\
             │ 0 1 0 0 1 1 │\n\
             │ 1 0 1 1 1 1 │\n\
             │ 0 1 0 1 1 0 │\n\
             │ 0 1 1 0 1 0 │\n\
             │ 1 1 1 1 0 1 │\n\
             │ 1 1 0 0 1 0 │\n\
             └             ┘"
        );
    }

    #[test]
    fn cuthill_mckee_requires_calc_degree() {
        // lin2_graph
        let mesh = Samples::graph_8_edges();
        let mut graph = GraphUnd::from_mesh(&mesh, false, false).unwrap();
        assert_eq!(
            graph.cuthill_mckee(Some(0)).err(),
            Some("Cuthill-McKee algorithm requires the degree of vertices (calc_degree must be set to true)")
        );
    }

    #[test]
    fn cuthill_mckee_works() {
        // lin2_graph
        let mesh = Samples::graph_8_edges();
        let mut graph = GraphUnd::from_mesh(&mesh, true, false).unwrap();
        let ordering = graph.cuthill_mckee(Some(0)).unwrap();
        // println!("ordering = {:?}", ordering);
        assert_eq!(ordering, &[7, 5, 6, 1, 3, 2, 4, 0]);
    }

    #[test]
    fn calc_distance_works_1() {
        // lin2_graph
        let mesh = Samples::graph_8_edges();
        let mut graph = GraphUnd::from_mesh(&mesh, false, false).unwrap();

        let max_distance = graph.calc_distance(0);
        assert_eq!(graph.distance, &[0, 3, 2, 2, 1, 4, 3, 4]);
        assert_eq!(max_distance, 4);

        let max_distance = graph.calc_distance(4);
        assert_eq!(graph.distance, &[1, 2, 1, 1, 0, 3, 2, 3]);
        assert_eq!(max_distance, 3);

        let max_distance = graph.calc_distance(5);
        assert_eq!(graph.distance, &[4, 1, 2, 4, 3, 0, 5, 1]);
        assert_eq!(max_distance, 5);

        let max_distance = graph.calc_distance(6);
        assert_eq!(graph.distance, &[3, 4, 3, 1, 2, 5, 0, 5]);
        assert_eq!(max_distance, 5);
    }

    #[test]
    fn calc_distance_works_2() {
        // 1-------0       7 ------6
        // |       |     .'|     .'|
        // |       |   .'  |   .'  |
        // |       | .'    | .'    |
        // 2       3-------4-------5
        let edges = [
            [0, 1], // 0
            [1, 2], // 1
            [3, 0], // 2
            [3, 4], // 3
            [4, 7], // 4
            [3, 7], // 5
            [7, 6], // 6
            [4, 5], // 7
            [6, 4], // 8
            [5, 6], // 9
        ];
        let mut graph = GraphUnd::from_edges(&edges, false, true).unwrap();

        let max_distance = graph.calc_distance(7);
        assert_eq!(graph.distance, &[2, 3, 4, 1, 1, 2, 1, 0]);
        assert_eq!(max_distance, 4);
    }

    #[test]
    fn shortest_path_bfs_works_1() {
        // 1-------0       7 ------6
        // |       |     .'|     .'|
        // |       |   .'  |   .'  |
        // |       | .'    | .'    |
        // 2       3-------4-------5
        let edges = [
            [0, 1], // 0
            [1, 2], // 1
            [3, 0], // 2
            [3, 4], // 3
            [4, 7], // 4
            [3, 7], // 5
            [7, 6], // 6
            [4, 5], // 7
            [6, 4], // 8
            [5, 6], // 9
        ];
        let mut graph = GraphUnd::from_edges(&edges, false, false).unwrap();

        let path = graph.shortest_path_bfs(0, 7);
        assert_eq!(path, &[0, 3, 7]);

        // Since calc_degree is false, the adjacency matrix is sorted by id only.
        // Thus, the results equal Mathematica's results.
        let path = graph.shortest_path_bfs(2, 6);
        assert_eq!(path, &[2, 1, 0, 3, 4, 6]);

        // Now, the adjacency matrix sorts by the degree first and then by the id, the node 7 is
        // selected instead of 4 because the degree of node 7 is 3 and the degree of node 4 is 4.
        let calc_degree = true; // will change the sorting of rows in the adjacency matrix
        let mut graph = GraphUnd::from_edges(&edges, calc_degree, false).unwrap();
        let path = graph.shortest_path_bfs(2, 6);
        assert_eq!(path, &[2, 1, 0, 3, 7, 6]);
    }

    #[test]
    fn pseudo_peripheral_requires_calc_degree() {
        // graph_8_edges
        let mesh = Samples::graph_8_edges();
        let mut graph = GraphUnd::from_mesh(&mesh, false, false).unwrap();
        assert_eq!(
            graph.pseudo_peripheral(None).err(),
            Some("pseudo_peripheral requires the degree of vertices (calc_degree must be set to true)")
        );
    }

    #[test]
    fn pseudo_peripheral_works() {
        // graph_8_edges
        let mesh = Samples::graph_8_edges();
        let mut graph = GraphUnd::from_mesh(&mesh, true, false).unwrap();
        assert_eq!(graph.pseudo_peripheral(None).unwrap(), 6);
        assert_eq!(graph.pseudo_peripheral(Some(4)).unwrap(), 6);
        assert_eq!(graph.pseudo_peripheral(Some(7)).unwrap(), 6);
        assert_eq!(graph.pseudo_peripheral(Some(6)).unwrap(), 5);

        // graph_12_edges
        let mesh = Samples::graph_12_edges();
        let mut graph = GraphUnd::from_mesh(&mesh, true, true).unwrap();
        assert_eq!(graph.pseudo_peripheral(Some(0)).unwrap(), 8);
        assert_eq!(graph.pseudo_peripheral(Some(4)).unwrap(), 2);
        assert_eq!(graph.pseudo_peripheral(None).unwrap(), 3);
    }

    #[test]
    fn test_density() {
        // Empty graph
        let edges: [[usize; 2]; 0] = [];
        let graph = GraphUnd::from_edges(&edges, false, false).unwrap();
        assert_eq!(graph.density(), 0.0);

        // Single node graph
        let edges = [[0, 0]];
        let graph = GraphUnd::from_edges(&edges, false, false).unwrap();
        assert_eq!(graph.density(), 0.0);

        // Complete graph K3
        let edges = [[0, 1], [1, 2], [0, 2]];
        let graph = GraphUnd::from_edges(&edges, false, false).unwrap();
        assert_eq!(graph.density(), 1.0);

        // Square graph (4 nodes, 4 edges)
        let edges = [[0, 1], [1, 2], [2, 3], [3, 0]];
        let graph = GraphUnd::from_edges(&edges, false, false).unwrap();
        assert_eq!(graph.density(), 4.0 / (4.0 * 3.0));
    }

    #[test]
    fn test_connected_components() {
        // Disconnected graph with three components
        let edges = [
            [0, 1], [1, 2],           // Component 1
            [3, 4],                   // Component 2
            [5, 6], [6, 7], [7, 5],  // Component 3
        ];
        let mut graph = GraphUnd::from_edges(&edges, false, false).unwrap();
        let components = graph.connected_components();
        
        assert_eq!(components.len(), 3);
        
        // Sort components and their contents for stable comparison
        let mut sorted_components: Vec<Vec<usize>> = components
            .into_iter()
            .map(|mut comp| {
                comp.sort_unstable();
                comp
            })
            .collect();
        sorted_components.sort_unstable();
        
        assert_eq!(sorted_components, vec![
            vec![0, 1, 2],
            vec![3, 4],
            vec![5, 6, 7],
        ]);
    }

    #[test]
    fn test_is_connected() {
        // Connected graph
        let edges = [[0, 1], [1, 2], [2, 3], [3, 0]];
        let mut graph = GraphUnd::from_edges(&edges, false, false).unwrap();
        assert!(graph.is_connected());

        // Disconnected graph
        let edges = [[0, 1], [2, 3]];
        let mut graph = GraphUnd::from_edges(&edges, false, false).unwrap();
        assert!(!graph.is_connected());

        // Single node is connected
        let edges = [[0, 0]];
        let mut graph = GraphUnd::from_edges(&edges, false, false).unwrap();
        assert!(graph.is_connected());
    }

    #[test]
    fn test_articulation_points() {
        // Graph with no articulation points (cycle)
        let edges = [[0, 1], [1, 2], [2, 3], [3, 0]];
        let graph = GraphUnd::from_edges(&edges, false, false).unwrap();
        assert!(graph.find_articulation_points().is_empty());

        // Graph with one articulation point
        let edges = [[0, 1], [1, 2], [2, 3], [1, 3], [1, 4]];
        let graph = GraphUnd::from_edges(&edges, false, false).unwrap();
        let mut art_points = graph.find_articulation_points();
        art_points.sort_unstable();
        assert_eq!(art_points, vec![1]);

        // Graph with multiple articulation points
        // 0 - 1 - 2 - 3 - 4
        //         |       |
        //         5 - 6 - 7
        let edges = [[0, 1], [1, 2], [2, 3], [3, 4], [2, 5], [5, 6], [6, 7], [7, 4]];
        let graph = GraphUnd::from_edges(&edges, false, false).unwrap();
        let mut art_points = graph.find_articulation_points();
        art_points.sort_unstable();
        assert_eq!(art_points, vec![2]);
    }

    #[test]
    fn gibbs_poole_stock_example() {
        // use graph example from:
        // Gibbs NW, Poole WG JR, and Stockmeyer PK (1976) An algorithm for reducing the bandwidth
        // and profile of a sparse matrix, SIAM Journal on Numerical Analysis, 13(2):236-250
        let mut block = Block::new(&[[0.0, 0.0], [5.0, 0.0], [5.0, 3.0], [0.0, 3.0]]).unwrap();
        block.set_ndiv(&[5, 3]).unwrap();
        let mut mesh = block.subdivide(GeoKind::Qua4).unwrap();
        let old_to_new = &[
            22, //  0
            20, //  1
            21, //  2
            2,  //  3
            10, //  4
            19, //  5
            3,  //  6
            18, //  7
            15, //  8
            7,  //  9
            16, // 10
            4,  // 11
            6,  // 12
            11, // 13
            0,  // 14
            17, // 15
            9,  // 16
            14, // 17
            12, // 18
            23, // 19
            5,  // 20
            13, // 21
            1,  // 22
            8,  // 23
        ];
        mesh.renumber_points(old_to_new).unwrap(); // this is to match the paper's numbers
        if SAVE_FIGURE {
            let mut fig = Figure::new();
            fig.show_point_ids(true);
            fig.draw(&mesh, "/tmp/gemlab/test_graph_gps_example.svg").unwrap();
        }
        let npoint = mesh.points.len();

        // original graph
        let mut graph = GraphUnd::from_mesh(&mesh, true, false).unwrap();
        let band = graph.calc_bandwidth();
        graph.print_non_zero_pattern();
        println!("band (original) = {}", band);
        assert_eq!(band, 22);

        // cuthill-mckee with fixed root = 8 (cm_8)
        let ordering = graph.cuthill_mckee(Some(8)).unwrap();

        // renumber mesh nodes
        let mut mesh_cm_8 = mesh.clone();
        let old_to_new = GraphUnd::get_old_to_new_map(&ordering);
        mesh_cm_8.renumber_points(&old_to_new).unwrap();

        // generate figure with levels/distance and mesh
        if SAVE_FIGURE {
            graph.calc_distance(8);
            for i in 0..npoint {
                mesh.points[i].marker = 1 + graph.distance[i] as i32; // use markers for the distance
            }
            let mut fig = Figure::new();
            fig.show_point_ids(true).show_point_marker(true);
            fig.draw(&mesh, "/tmp/gemlab/test_graph_gps_example_cm_8.svg").unwrap();
        }

        // print pattern with updated mesh (cm_8)
        let graph_cm_8 = GraphUnd::from_mesh(&mesh_cm_8, true, false).unwrap();
        let band = graph_cm_8.calc_bandwidth();
        graph_cm_8.print_non_zero_pattern();
        println!("band (cm_8) = {}", band);
        assert_eq!(band, 9);

        // CM algo with pseudo-peripheral root
        let mut graph = GraphUnd::from_mesh(&mesh, true, false).unwrap();

        // renumber mesh nodes (cuthill-mckee + pseudo-peripheral)
        let mut mesh_cm_pp = mesh.clone();
        GraphUnd::renumber_mesh(&mut mesh_cm_pp, false).unwrap();

        // generate figure with levels/distance and mesh
        if SAVE_FIGURE {
            let root = graph.pseudo_peripheral(None).unwrap();
            graph.calc_distance(root);
            for i in 0..npoint {
                mesh.points[i].marker = 1 + graph.distance[i] as i32; // use markers for the distance
            }
            let mut fig = Figure::new();
            fig.show_point_ids(true).show_point_marker(true);
            fig.draw(&mesh, "/tmp/gemlab/test_graph_gps_example_cm_pp.svg").unwrap();
        }

        // print pattern with updated mesh (cm_pp)
        let graph_cm_pp = GraphUnd::from_mesh(&mesh_cm_pp, true, false).unwrap();
        let band = graph_cm_pp.calc_bandwidth();
        graph_cm_pp.print_non_zero_pattern();
        println!("band (cm_pp) = {}", band);
        assert_eq!(band, 8);
    }
}
