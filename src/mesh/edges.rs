use super::Edges;

impl<'a> Edges<'a> {
    /// Returns indices of edges in connected order
    ///
    /// # Returns
    ///
    /// * Returns a vector of indices into `self.all` representing a path through connected edges
    /// * The path starts at index 0 and follows connected edges
    /// * Returns `None` if the edges are not all connected or if there are branches
    ///
    /// # Examples
    ///
    /// ```
    /// use gemlab::mesh::Features;
    /// use gemlab::mesh::Samples;
    /// use gemlab::mesh::At;
    /// use gemlab::util::any_x;
    ///
    /// let mesh = Samples::two_qua4();
    /// let features = Features::new(&mesh, false);
    ///
    /// // Get edges along bottom of mesh (y=0)
    /// let edges = features.search_edges(At::Y(0.0), any_x).unwrap();
    ///
    /// // Get connected path
    /// let path = edges.connected_path().unwrap();
    /// assert_eq!(path, vec![0, 1]); // Indices into edges.all
    /// ```
    pub fn connected_path(&self) -> Option<Vec<usize>> {
        if self.all.is_empty() {
            return Some(Vec::new());
        }

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
                            return None;
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
            return None;
        }

        Some(path)
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use crate::mesh::{At, Features, Samples};
    use crate::util::any_x;

    #[test]
    fn test_connected_edges_path() {
        // Create test mesh with some edges (QUA8 elements)
        // Note: Edges have 3 nodes each
        //
        //  14---16---13---20---18
        //   |         |         |
        //  17        15        19
        //   |         |         |
        //   3----6----2---12----9
        //   |         |         |
        //   7         5        11
        //   |         |         |
        //   0----4----1---10----8
        //
        let mesh = Samples::block_2d_four_qua8();
        let feat = Features::new(&mesh, true);

        // Test bottom horizontal path (y = 0.0)
        let edges = feat.search_edges(At::Y(0.0), any_x).unwrap();
        let path = edges.connected_path().unwrap();
        assert_eq!(path, vec![0, 1]); // Edge path: (0-4-1), (1-10-8)

        // Test top horizontal path (y = 2.0)
        let edges = feat.search_edges(At::Y(2.0), any_x).unwrap();
        let path = edges.connected_path().unwrap();
        assert_eq!(path, vec![0, 1]); // Edge path: (14-16-13), (13-20-18)

        // Test left vertical path (x = 0.0)
        let edges = feat.search_edges(At::X(0.0), any_x).unwrap();
        let path = edges.connected_path().unwrap();
        assert_eq!(path, vec![0, 1]); // Edge path: (0-7-3), (3-17-14)

        // Test disconnected edges should return None
        let edges = feat.search_many_edges(&[At::X(0.0), At::X(2.0)], any_x).unwrap();
        assert!(edges.connected_path().is_none());

        // Test branching edges should return None
        let edges = feat
            .search_many_edges(&[At::Y(0.0), At::Y(1.0), At::X(0.0)], any_x)
            .unwrap();
        assert!(edges.connected_path().is_none());
    }
}
