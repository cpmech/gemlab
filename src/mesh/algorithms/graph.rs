use crate::mesh::{Mesh, PointId};
use crate::StrError;
use std::collections::{HashSet, VecDeque};

/// Defines a graph structure to assist in renumbering a mesh
pub struct Graph {
    /// Holds the adjacency (sparse) matrix (point connections)
    ///
    /// Note: each row in this matrix is sorted in ascending order of degree, followed by id
    pub adjacency: Vec<Vec<PointId>>,

    /// Holds all point degrees (the number of connections of a vertex)
    pub degree: Vec<usize>,

    /// Holds the id of the point with the minimum degree
    pub p_min_degree: usize,

    /// Holds the id of the point with the maximum degree
    pub p_max_degree: usize,

    /// Defines an auxiliary queue for BFS (breadth-first-search) runs
    ///
    /// Note: must be cleared before each use.
    queue: VecDeque<PointId>,

    /// Holds the auxiliary list of bool indicating that a vertex has been explored
    ///
    /// Note: must be cleared before each use.
    explored: Vec<bool>,

    /// Holds the distance from the root to each point
    ///
    /// Note: must be cleared before each use.
    distance: Vec<usize>,
}

impl Graph {
    /// Allocates a new instance
    ///
    /// **Note:** The graph must be connected.
    pub fn new(mesh: &Mesh) -> Result<Self, StrError> {
        // find the adjacency (sparse) matrix of nodes' connections
        let npoint = mesh.points.len();
        let mut adjacency_set = vec![HashSet::new(); npoint];
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

        // check the connectivity of the graph (all vertices must be explored)
        let mut queue = VecDeque::new();
        let mut explored = vec![false; npoint];
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

        // compute the list of degrees (number of connections of a vertex)
        let mut degree = vec![0; npoint];
        let mut p_min_degree = 0;
        let mut p_max_degree = 0;
        let mut min_degree = usize::MAX;
        let mut max_degree = 0;
        for (p, row) in adjacency_set.iter().enumerate() {
            degree[p] = row.len();
            if degree[p] < min_degree {
                p_min_degree = p;
                min_degree = degree[p];
            }
            if degree[p] > max_degree {
                p_max_degree = p;
                max_degree = degree[p];
            }
        }

        // sort each row of the adjacency matrix by the degree, and then by id
        let mut adjacency = Vec::new();
        for row_set in adjacency_set.iter() {
            let mut row: Vec<_> = row_set.iter().copied().collect();
            row.sort_unstable_by_key(|p| (degree[*p], *p));
            adjacency.push(row);
        }

        // results
        Ok(Graph {
            adjacency,
            degree,
            p_min_degree,
            p_max_degree,
            queue,
            explored,
            distance: vec![0; npoint],
        })
    }

    /// Computes the ordering array to renumber the vertices according to the (reverse) Cuthill-McKee algorithm
    ///
    /// **Note:** All nodes must be reachable from the root; i.e., the corresponding graph must be connected.
    ///
    /// # Input
    ///
    /// * `start_point` -- (root) the first point id, which will not be renumbered. Should have a low degree.
    ///   If None, a point with a minimum degree will be used.
    ///
    /// # Output
    ///
    /// Returns the ordering array such that `ordering[new_point_id] = old_point_id`
    pub fn cuthill_mckee(&mut self, start_point: Option<PointId>) -> Result<Vec<PointId>, StrError> {
        // root point
        let root = match start_point {
            Some(p) => p,
            None => self.p_min_degree,
        };

        // clear auxiliary structures
        self.queue.clear();
        let npoint = self.adjacency.len();
        for i in 0..npoint {
            self.explored[i] = false;
        }

        // allocate auxiliary structures
        let mut ordering = vec![0; npoint];

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
        let npoint = self.adjacency.len();
        for i in 0..npoint {
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
        for i in 0..npoint {
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
    pub fn pseudo_peripheral(&mut self, start_point: Option<PointId>) -> usize {
        // root point
        let mut root = match start_point {
            Some(p) => p,
            None => self.p_min_degree,
        };

        // first distances
        let mut max_distance = self.calc_distance(root);

        // auxiliary
        const MAX_ITERATIONS: usize = 10;
        let npoint = self.adjacency.len();
        let mut success = false;

        // perform iterations
        for _ in 0..MAX_ITERATIONS {
            // loop over all points, consider only the points with max distance
            // equal to root's max distance, and select the point with the minimum degree
            let mut next_root = None;
            let mut min_deg = usize::MAX; // min degree among potential roots
            for i in 0..npoint {
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
            panic!("INTERNAL ERROR: iterations did not converge");
        }
        root
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::Graph;
    use crate::mesh::{Cell, Mesh, Point, Samples};
    use crate::shapes::GeoKind;
    use russell_lab::NumMatrix;

    #[test]
    fn graph_new_works_1() {
        // lin2_graph
        let mesh = Samples::graph_8_edges();
        let graph = Graph::new(&mesh).unwrap();

        //                         0  1  2  3  4  5  6  7 (point)
        assert_eq!(graph.degree, &[1, 3, 2, 2, 3, 2, 1, 2]);
        assert_eq!(graph.p_min_degree, 0);
        assert_eq!(graph.p_max_degree, 1);

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
    fn graph_new_works_2() {
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
        let graph = Graph::new(&mesh).unwrap();

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
    fn cuthill_mckee_works() {
        // lin2_graph
        let mesh = Samples::graph_8_edges();
        let mut graph = Graph::new(&mesh).unwrap();
        let ordering = graph.cuthill_mckee(Some(0)).unwrap();
        // println!("ordering = {:?}", ordering);
        assert_eq!(ordering, &[7, 5, 6, 1, 3, 2, 4, 0]);
    }

    #[test]
    fn distance_works() {
        // lin2_graph
        let mesh = Samples::graph_8_edges();
        let mut graph = Graph::new(&mesh).unwrap();

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
    fn pseudo_peripheral_works() {
        // graph_8_edges
        let mesh = Samples::graph_8_edges();
        let mut graph = Graph::new(&mesh).unwrap();
        assert_eq!(graph.pseudo_peripheral(None), 6);
        assert_eq!(graph.pseudo_peripheral(Some(4)), 6);
        assert_eq!(graph.pseudo_peripheral(Some(7)), 6);
        assert_eq!(graph.pseudo_peripheral(Some(6)), 5);

        // graph_12_edges
        let mesh = Samples::graph_12_edges();
        let mut graph = Graph::new(&mesh).unwrap();
        assert_eq!(graph.pseudo_peripheral(Some(0)), 8);
        assert_eq!(graph.pseudo_peripheral(Some(4)), 2);
        assert_eq!(graph.pseudo_peripheral(None), 3);
    }
}
