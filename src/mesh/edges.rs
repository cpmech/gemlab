#![allow(unused)]

use super::Edges;
use crate::StrError;

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
    pub fn connected_path(&self) -> Result<Vec<usize>, StrError> {
        if self.all.is_empty() {
            return Err("the edges list is empty");
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
                            return Err("found branching");
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
    use std::vec;

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
    fn test_connected_edges_path_1() {
        let all = generate_sample();

        let empty = Edges { all: vec![] };
        assert_eq!(empty.connected_path().err(), Some("the edges list is empty"));

        let branching = Edges {
            all: vec![&all[5], &all[6], &all[7], &all[8]],
        };
        assert_eq!(branching.connected_path().err(), Some("found branching"));

        let disconnected = Edges {
            all: vec![&all[3], &all[4], &all[10]],
        };
        assert_eq!(disconnected.connected_path().err(), Some("found disconnected edges"));

        let bottom = Edges {
            all: vec![&all[0], &all[1], &all[2]],
        };
        assert_eq!(bottom.connected_path().unwrap(), vec![0, 1, 2]);

        let loop1 = Edges {
            all: vec![&all[5], &all[7], &all[9], &all[1]],
        };
        assert_eq!(loop1.connected_path().unwrap(), vec![0, 1, 2, 3]);
    }
}
