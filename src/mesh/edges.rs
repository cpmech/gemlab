use super::Edges;
use crate::StrError;

impl<'a> Edges<'a> {
    /// Returns indices of edges in connected order, starting from index 0
    ///
    /// A connected path is a sequence of edges where:
    /// * Each edge shares exactly one vertex with the next edge
    /// * No edge connects to more than two other edges (no branches)
    /// * All edges must form a single continuous path
    ///
    /// # Returns
    ///
    /// * `Ok(Vec<usize>)` - Vector of indices into `self.all` representing the connected path
    /// * `Err(StrError)` - If the edges:
    ///   * Form a branching path (an edge connects to more than two others)
    ///   * Are disconnected (not all edges can be reached)
    ///   * Form a loop
    ///   * The list is empty
    ///
    /// # Examples
    ///
    /// ```
    /// use gemlab::mesh::Features;
    /// use gemlab::mesh::Samples;
    /// use gemlab::mesh::At;
    /// use gemlab::util::any_x;
    ///
    /// // 1.0  3-----------2-----------5
    /// //      |           |           |
    /// //      |           |           |
    /// //      |           |           |
    /// //      |           |           |
    /// // 0.0  0-----------1-----------4
    /// //     0.0         1.0         2.0
    /// let mesh = Samples::two_qua4();
    /// let feat = Features::new(&mesh, false);
    ///
    /// // Valid connected path along bottom edge (y=0)
    /// let edges = feat.search_edges(At::Y(0.0), any_x).unwrap();
    /// assert_eq!(edges.connected_path().unwrap(), vec![0, 1]);
    ///
    /// // Disconnected edges will error
    /// let edges = feat.search_many_edges(
    ///     &[At::Y(0.0), At::Y(1.0)],
    ///     any_x
    /// ).unwrap();
    /// assert!(edges.connected_path().is_err());
    ///
    /// // Branching edges will error
    /// let edges = feat.search_many_edges(
    ///     &[At::Y(0.0), At::X(0.0)],
    ///     any_x
    /// ).unwrap();
    /// assert!(edges.connected_path().is_err());
    /// ```
    pub fn connected_path(&self) -> Result<Vec<usize>, StrError> {
        // Algorithm:
        // 1. Start with edge at index 0
        // 2. Find the next unvisited edge that shares a vertex
        // 3. Continue until all edges are visited or an error is found
        // 4. Return error if:
        //    * Multiple edges connect to current edge (branch/loop)
        //    * Not all edges were reached (disconnected)

        // Check if empty
        if self.all.is_empty() {
            return Err("the edges list is empty");
        }

        // First edge
        let mut path = Vec::with_capacity(self.all.len());
        let mut visited = vec![false; self.all.len()];
        let mut current = 0;

        // Add first edge
        path.push(current);
        visited[current] = true;

        // Find connected edges
        while path.len() < self.all.len() {
            let (a, b) = self.all[current].key();
            let mut found = false;

            // Look for next connected edge
            for (i, edge) in self.all.iter().enumerate() {
                if !visited[i] {
                    let (c, d) = edge.key();
                    if a == c || a == d || b == c || b == d {
                        if found {
                            // More than one connection means there's a branch
                            return Err("found branching (or loop)");
                        }
                        current = i;
                        path.push(current);
                        visited[i] = true;
                        found = true;
                    }
                }
            }

            if !found {
                break;
            }
        }

        // All edges must be connected
        if path.len() != self.all.len() {
            return Err("found disconnected edges");
        }

        Ok(path)
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use crate::mesh::{At, Edge, Edges, Features, Samples};
    use crate::shapes::GeoKind;
    use crate::util::any_x;

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
    fn test_connected_path_empty() {
        let empty = Edges { all: vec![] };
        assert_eq!(empty.connected_path().err(), Some("the edges list is empty"));
    }

    #[test]
    fn test_connected_path_single_edge() {
        let all = generate_sample();
        let single = Edges { all: vec![&all[0]] };
        assert_eq!(single.connected_path().unwrap(), vec![0]);
    }

    #[test]
    fn test_connected_path_branching() {
        let all = generate_sample();

        // Test branching at node 8
        let branching = Edges {
            all: vec![&all[5], &all[6], &all[7], &all[8]],
        };
        assert_eq!(branching.connected_path().err(), Some("found branching (or loop)"));

        // Test branching at node 4
        let branching_at_4 = Edges {
            all: vec![&all[7], &all[9], &all[12], &all[13]],
        };
        assert_eq!(branching_at_4.connected_path().err(), Some("found branching (or loop)"));
    }

    #[test]
    fn test_connected_path_disconnected() {
        let all = generate_sample();

        // Test completely disconnected edges
        let disconnected = Edges {
            all: vec![&all[3], &all[4], &all[10]],
        };
        assert_eq!(disconnected.connected_path().err(), Some("found disconnected edges"));

        // Test partially connected edges
        let partially_connected = Edges {
            all: vec![&all[0], &all[1], &all[10], &all[11]],
        };
        assert_eq!(
            partially_connected.connected_path().err(),
            Some("found disconnected edges")
        );
    }

    #[test]
    fn test_connected_path_loops() {
        let all = generate_sample();

        // Test square loop
        let square_loop = Edges {
            all: vec![&all[5], &all[7], &all[9], &all[1]],
        };
        assert_eq!(square_loop.connected_path().err(), Some("found branching (or loop)"));

        // Test triangle loop - edges connecting nodes 8-4, 4-7, and 7-8
        let triangle_loop = Edges {
            all: vec![&all[7], &all[12], &all[6]],
        };
        assert_eq!(triangle_loop.connected_path().err(), Some("found branching (or loop)"));
    }

    #[test]
    fn test_connected_path_valid_paths() {
        let all = generate_sample();

        // Test bottom horizontal path
        let bottom = Edges {
            all: vec![&all[0], &all[1], &all[2]],
        };
        assert_eq!(bottom.connected_path().unwrap(), vec![0, 1, 2]);

        // Test left vertical path
        let left = Edges {
            all: vec![&all[3], &all[4]],
        };
        assert_eq!(left.connected_path().unwrap(), vec![0, 1]);

        // Test diagonal path
        let diagonal = Edges {
            all: vec![&all[5], &all[7]],
        };
        assert_eq!(diagonal.connected_path().unwrap(), vec![0, 1]);
    }

    #[test]
    fn test_connected_path_with_mesh() {
        let mesh = Samples::block_2d_four_qua8();
        let feat = Features::new(&mesh, true);

        // Test bottom horizontal path (y = 0.0)
        let edges = feat.search_edges(At::Y(0.0), any_x).unwrap();
        assert_eq!(edges.connected_path().unwrap(), vec![0, 1]);

        // Test left vertical path (x = 0.0)
        let edges = feat.search_edges(At::X(0.0), any_x).unwrap();
        assert_eq!(edges.connected_path().unwrap(), vec![0, 1]);

        // Test disconnected edges should error
        let edges = feat.search_many_edges(&[At::X(0.0), At::X(2.0)], any_x).unwrap();
        assert!(edges.connected_path().is_err());
    }
}
