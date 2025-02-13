use super::{Edges, PointId};
use crate::shapes::GeoKind;
use std::collections::HashMap;

impl<'a> Edges<'a> {
    /// Finds a sequence of edges by following connected points to create a path
    ///
    /// # Returns
    ///
    /// Returns a tuple `(path, points)` where:
    ///
    /// * `path` -- List of indices into `self.all` representing edges along the path
    /// * `points` -- List of all points along the path. The list includes the middle points of high-order edges.
    ///
    /// # Path Finding Algorithm
    ///
    /// 1. Start point selection:
    ///    * First tries to find points connected to only one edge (endpoints)
    ///    * If no endpoints exist (e.g., in a loop), uses point with lowest ID
    ///    * If multiple endpoints exist, uses endpoint with lowest ID
    ///
    /// 2. Path construction:
    ///    * Starts from selected point and follows connected edges
    ///    * At each step, takes first available edge not yet used
    ///    * Stops when no more edges can be followed or a loop is detected
    ///
    /// # Notes
    ///
    /// * Does not guarantee finding the longest path (NP-hard problem)
    /// * For simple connected paths (no branches), will find the only possible path
    /// * For branching paths, result depends on edge order in `self.all`
    /// * For loops, starts at lowest numbered point and follows first available edge
    /// * For disconnected components, follows path until no more connected edges
    ///
    /// # Special Cases
    ///
    /// * Empty edge list returns `(Vec::new(), Vec::new())`
    /// * Single edge returns path with just that edge
    /// * Loop starts from lowest numbered point
    /// * Branching paths follow first available edge at each step
    /// * Disconnected edges follow path until component ends
    ///
    /// # Examples
    ///
    /// ```
    /// use gemlab::mesh::{Edge, Edges};
    /// use gemlab::shapes::GeoKind;
    ///
    /// // Create some sample edges (a simple path 1-2-3)
    /// let e1 = Edge { kind: GeoKind::Lin2, points: vec![1, 2] };
    /// let e2 = Edge { kind: GeoKind::Lin2, points: vec![2, 3] };
    /// let edges = Edges { all: vec![&e1, &e2] };
    ///
    /// // Get path through edges
    /// let (path, points) = edges.any_path();
    /// assert_eq!(path, vec![0, 1]);      // Edge indices
    /// assert_eq!(points, vec![1, 2, 3]); // Point indices
    /// ```
    pub fn any_path(&self) -> (Vec<usize>, Vec<usize>) {
        // check if the list of edges is empty
        if self.all.is_empty() {
            return (Vec::new(), Vec::new());
        }

        // maps points to edges (indices in self.all)
        let mut nnode_per_edge_max = 2;
        let mut map: HashMap<PointId, Vec<usize>> = HashMap::new();
        for e in 0..self.all.len() {
            let edge = &self.all[e];
            let (a, b) = edge.key();
            map.entry(a).or_insert_with(Vec::new).push(e);
            map.entry(b).or_insert_with(Vec::new).push(e);
            if edge.points.len() > nnode_per_edge_max {
                nnode_per_edge_max = edge.points.len();
            }
        }

        // define an array of possible endpoints (those shared by a single edge)
        let mut endpoints: Vec<_> = map
            .iter()
            .filter_map(|(point, edges)| if edges.len() == 1 { Some(point) } else { None })
            .collect();

        // sort the list of possible endpoints by point id
        endpoints.sort();

        // start with the first point (e.g., a point in a loop will do)
        let endpoint = if endpoints.len() > 0 {
            // use the endpoint with the lowest point id among all extremities
            *endpoints[0]
        } else {
            // use the point with the lowest point id among all points
            let mut all_points: Vec<_> = map.keys().copied().collect();
            all_points.sort();
            all_points[0]
        };

        // define a helper function to get the next point, given the current point and edge
        let next_point = |point: PointId, edge: usize| {
            let (a, b) = (self.all[edge].points[0], self.all[edge].points[1]);
            if a == point {
                b
            } else {
                a
            }
        };

        // define a helper function to get the next array of edges, given the current point and edge
        let next_edges = |point: PointId, edge: usize| {
            let edges = map.get(&point).unwrap();
            let res: Vec<_> = edges.iter().filter(|e| **e != edge).collect();
            res
        };

        // select current edge, start, and next point
        let mut e = map.get(&endpoint).unwrap()[0];
        let mut a = endpoint;
        let mut b = next_point(endpoint, e);

        // allocate output arrays
        let mut path = Vec::with_capacity(self.all.len());
        let mut points = Vec::with_capacity(self.all.len() * nnode_per_edge_max);
        path.push(e);
        points.push(a);

        // helper function to update the points array with middle and last edge points of an edge
        let mut push_other_points = |e: usize, a: PointId| match self.all[e].kind {
            GeoKind::Lin2 => {
                if self.all[e].points[0] == a {
                    points.push(self.all[e].points[1]);
                } else {
                    points.push(self.all[e].points[0]);
                }
            }
            GeoKind::Lin3 => {
                if self.all[e].points[0] == a {
                    points.push(self.all[e].points[2]);
                    points.push(self.all[e].points[1]);
                } else {
                    points.push(self.all[e].points[2]);
                    points.push(self.all[e].points[0]);
                }
            }
            GeoKind::Lin4 => {
                if self.all[e].points[0] == a {
                    points.push(self.all[e].points[2]);
                    points.push(self.all[e].points[3]);
                    points.push(self.all[e].points[1]);
                } else {
                    points.push(self.all[e].points[3]);
                    points.push(self.all[e].points[2]);
                    points.push(self.all[e].points[0]);
                }
            }
            GeoKind::Lin5 => {
                if self.all[e].points[0] == a {
                    points.push(self.all[e].points[3]);
                    points.push(self.all[e].points[2]);
                    points.push(self.all[e].points[4]);
                    points.push(self.all[e].points[1]);
                } else {
                    points.push(self.all[e].points[4]);
                    points.push(self.all[e].points[2]);
                    points.push(self.all[e].points[3]);
                    points.push(self.all[e].points[0]);
                }
            }
            _ => panic!("the edge kind must be Lin"),
        };

        // push other points of the first edge, noting that the first point is already in the list
        push_other_points(e, a);

        // follow path
        for _ in 0..self.all.len() {
            let edges = next_edges(b, e);
            if edges.len() == 0 {
                // no more edges to follow
                break;
            }
            e = *edges[0];
            a = b;
            b = next_point(b, e);
            path.push(e);
            push_other_points(e, a);
            if b == endpoint {
                // loop detected
                break;
            }
        }

        (path, points)
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use crate::mesh::{Edge, Edges};
    use crate::shapes::GeoKind;

    #[rustfmt::skip]
    fn generate_sample() -> Vec<Edge> {
        //          (11)           (12)           (13)
        //    9--------------7-------------.4--------------6
        //    |              |           .' |              |
        //    |              |         .'   |              |
        // (4)|           (6)|   (7) .'     |              |
        //    |              |     .'       |              |
        //    |              |   .'         |              |
        //    |      (8)     | .'           |(9)           |(10)
        //    1--------------8'             |              |
        //    |              |              |              |
        //    |              |              |              |
        // (3)|           (5)|              |              |
        //    |              |              |              |
        //    |              |              |              |
        //    5--------------2--------------3-------------10
        //           (0)           (1)            (2)

        // Bottom horizontal edges
        let e0 = Edge { kind: GeoKind::Lin2, points: vec![5, 2] };
        let e1 = Edge { kind: GeoKind::Lin2, points: vec![2, 3] };
        let e2 = Edge { kind: GeoKind::Lin2, points: vec![3, 10] };

        // Left vertical edges
        let e3 = Edge { kind: GeoKind::Lin2, points: vec![5, 1] };
        let e4 = Edge { kind: GeoKind::Lin2, points: vec![1, 9] };

        // Central edges
        let e5 = Edge { kind: GeoKind::Lin2, points: vec![2, 8] };
        let e6 = Edge { kind: GeoKind::Lin2, points: vec![8, 7] };
        let e7 = Edge { kind: GeoKind::Lin2, points: vec![8, 4] };
        let e8 = Edge { kind: GeoKind::Lin2, points: vec![1, 8] };

        // Right vertical edges
        let e9  = Edge { kind: GeoKind::Lin2, points: vec![3, 4] };
        let e10 = Edge { kind: GeoKind::Lin2, points: vec![10, 6] };

        // Top horizontal edges
        let e11 = Edge { kind: GeoKind::Lin2, points: vec![9, 7] };
        let e12 = Edge { kind: GeoKind::Lin2, points: vec![7, 4] };
        let e13 = Edge { kind: GeoKind::Lin2, points: vec![4, 6] };

        vec![e0, e1, e2, e3, e4, e5, e6, e7, e8, e9, e10, e11, e12, e13]
    }

    #[test]
    fn test_any_path_1() {
        // Allocate edges
        let all = generate_sample();

        // Empty list of edges
        let empty = Edges { all: vec![] };
        assert_eq!(empty.any_path(), (Vec::new(), Vec::new()));

        // Branching
        // The selected endpoint is 1 on edge 8 because it is the lowest among [1,2,4,7]
        // The first edge is 8 with index 3 in `all`
        // The second edge is 5 with index 0 in `all` because it is the first in `all` which is not 8
        let branching = Edges {
            all: vec![&all[5], &all[7], &all[6], &all[8]],
        };
        assert_eq!(branching.any_path(), (vec![/*8*/ 3, /*5*/ 0], vec![1, 8, 2]));

        // Loop 1
        // The selected endpoint is 2 on edge 5 because it is the lowest among all points and there aren't extremities
        // The first edge is 5 with index 0 in `all` because it comes before edge 1 in `all`
        let loop1 = Edges {
            all: vec![&all[5], &all[9], &all[1], &all[7]], // 5 come first => clockwise loop
        };
        assert_eq!(
            loop1.any_path(),
            (vec![/*5*/ 0, /*7*/ 3, /*9*/ 1, /*1*/ 2], vec![2, 8, 4, 3, 2])
        );

        // Loop 2
        // The selected endpoint is 2 on edge 1 because it is the lowest among all points and there aren't extremities
        // The first edge is 1 with index 2 in `all` because it comes before edge 5 in `all`
        let loop2 = Edges {
            all: vec![&all[7], &all[9], &all[1], &all[5]], // 7 come first => counter-clockwise loop
        };
        assert_eq!(
            loop2.any_path(),
            (vec![/*1*/ 2, /*9*/ 1, /*7*/ 0, /*5*/ 3], vec![2, 3, 4, 8, 2])
        );

        // Disconnected
        // The selected endpoint is 5 on edge 3 because it is the lowest among [5,6,9,10]
        // The first edge is 3 with index 2 in `all`
        let disconnected = Edges {
            all: vec![&all[10], &all[4], &all[3]],
        };
        assert_eq!(disconnected.any_path(), (vec![/*3*/ 2, /*4*/ 1], vec![5, 1, 9]));

        // Bottom edges
        // The selected endpoint is 5 on edge 0 because it is the lowest among [5,10]
        // The first edge is 0 with index 1 in `all`
        let bottom = Edges {
            all: vec![&all[2], &all[0], &all[1]],
        };
        assert_eq!(bottom.any_path(), (vec![/*0*/ 1, /*1*/ 2, /*2*/ 0], vec![5, 2, 3, 10]));
    }

    #[rustfmt::skip]
    fn generate_sample_o2() -> Vec<Edge> {
        //           (11)           (12)           (13)
        //     9------105-----7-----106-----.4------107-----6
        //     |              |           .' |              |
        //     |              |         .'   |              |
        // (4)112         (6)113  (7) .'     |              |
        //     |              |     .104     |              |
        //     |              |   .'         |              |
        //     |      (8)     | .'           |(9)           |(10)
        //     1------103-----8'            110            111
        //     |              |              |              |
        //     |              |              |              |
        // (3)108         (5)109             |              |
        //     |              |              |              |
        //     |              |              |              |
        //     5------100-----2------101-----3-----102-----10
        //            (0)            (1)            (2)

        // Bottom horizontal edges
        let e0 = Edge { kind: GeoKind::Lin3, points: vec![5, 2, 100] };
        let e1 = Edge { kind: GeoKind::Lin3, points: vec![2, 3, 101] };
        let e2 = Edge { kind: GeoKind::Lin3, points: vec![3, 10, 102] };

        // Left vertical edges
        let e3 = Edge { kind: GeoKind::Lin3, points: vec![5, 1, 108] };
        let e4 = Edge { kind: GeoKind::Lin3, points: vec![1, 9, 112] };

        // Central edges
        let e5 = Edge { kind: GeoKind::Lin3, points: vec![2, 8, 109] };
        let e6 = Edge { kind: GeoKind::Lin3, points: vec![8, 7, 113] };
        let e7 = Edge { kind: GeoKind::Lin3, points: vec![8, 4, 104] };
        let e8 = Edge { kind: GeoKind::Lin3, points: vec![1, 8, 103] };

        // Right vertical edges
        let e9  = Edge { kind: GeoKind::Lin3, points: vec![3, 4, 110] };
        let e10 = Edge { kind: GeoKind::Lin3, points: vec![10, 6, 111] };

        // Top horizontal edges
        let e11 = Edge { kind: GeoKind::Lin3, points: vec![9, 7, 105] };
        let e12 = Edge { kind: GeoKind::Lin3, points: vec![7, 4, 106] };
        let e13 = Edge { kind: GeoKind::Lin3, points: vec![4, 6, 107] };

        vec![e0, e1, e2, e3, e4, e5, e6, e7, e8, e9, e10, e11, e12, e13]
    }

    #[test]
    fn test_any_path_2() {
        // Allocate edges
        let all = generate_sample_o2();

        // Branching
        let branching = Edges {
            all: vec![&all[5], &all[7], &all[6], &all[8]],
        };
        assert_eq!(branching.any_path(), (vec![/*8*/ 3, /*5*/ 0], vec![1, 103, 8, 109, 2]));

        // Loop 1
        let loop1 = Edges {
            all: vec![&all[5], &all[9], &all[1], &all[7]], // 5 come first => clockwise loop
        };
        assert_eq!(
            loop1.any_path(),
            (
                vec![/*5*/ 0, /*7*/ 3, /*9*/ 1, /*1*/ 2],
                vec![2, 109, 8, 104, 4, 110, 3, 101, 2]
            )
        );

        // Loop 2
        let loop2 = Edges {
            all: vec![&all[7], &all[9], &all[1], &all[5]], // 7 come first => counter-clockwise loop
        };
        assert_eq!(
            loop2.any_path(),
            (
                vec![/*1*/ 2, /*9*/ 1, /*7*/ 0, /*5*/ 3],
                vec![2, 101, 3, 110, 4, 104, 8, 109, 2]
            )
        );

        // Disconnected
        let disconnected = Edges {
            all: vec![&all[10], &all[4], &all[3]],
        };
        assert_eq!(
            disconnected.any_path(),
            (vec![/*3*/ 2, /*4*/ 1], vec![5, 108, 1, 112, 9])
        );

        // Bottom edges
        let bottom = Edges {
            all: vec![&all[2], &all[0], &all[1]],
        };
        assert_eq!(
            bottom.any_path(),
            (vec![/*0*/ 1, /*1*/ 2, /*2*/ 0], vec![5, 100, 2, 101, 3, 102, 10])
        );
    }
}
