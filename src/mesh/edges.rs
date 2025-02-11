#![allow(unused)]

use super::{Edge, Edges, PointId};
use crate::StrError;
use std::{collections::HashMap, vec};

impl<'a> Edges<'a> {
    /// Finds a sequence of points by traversing the edges
    ///
    /// Returns a list of points that are connected by the edges. The first point is found such that it has
    /// only one edge connected to it. If no such point is found (e.g., in a loop), the point with the
    /// **lowest number** among all points is selected as the first point. If there are multiple extremities
    /// (e.g. in the case of  branching), the point with the **lowest number** among all extremities is
    /// selected as the first point.
    ///
    /// Ideally, this function would find the longest path. However, the longest path problem is NP-hard!
    /// Therefore, we use a simple algorithm that finds **any path**. If the path is simply connected,
    /// i.e., without branching or disconnected edges, the path will be the longest path (obviously, because
    /// there is only one path).
    pub fn any_path(&self) -> Result<(Vec<PointId>, Vec<usize>), StrError> {
        // check if the list of edges is empty
        if self.all.is_empty() {
            return Err("the list of edges is empty");
        }

        // maps points to edges (indices in self.all)
        let mut map: HashMap<PointId, Vec<usize>> = HashMap::new();
        for e in 0..self.all.len() {
            let edge = &self.all[e];
            let (a, b) = edge.key();
            map.entry(a).or_insert_with(Vec::new).push(e);
            map.entry(b).or_insert_with(Vec::new).push(e);
        }

        println!();
        for (point, edges) in &map {
            println!(
                "point = {:?}, edges = {:?}",
                point,
                edges.iter().map(|e| self.all[*e].points.clone()).collect::<Vec<_>>()
            );
        }

        // define an array of possible endpoints (those shared by a single edge)
        let mut endpoints: Vec<_> = map
            .iter()
            .filter_map(|(point, edges)| if edges.len() == 1 { Some(point) } else { None })
            .collect();

        // sort the list of possible endpoints by point id
        endpoints.sort();

        println!("\nendpoints     = {:?}", endpoints);

        // start with the first point (e.g., a point in a loop will do)
        let mut endpoint = if endpoints.len() > 0 {
            // use the endpoint with the lowest point id among all extremities
            *endpoints[0]
        } else {
            // use the point with the lowest point id among all points
            let mut all_points: Vec<_> = map.keys().copied().collect();
            all_points.sort();
            all_points[0]
        };

        println!("endpoint      = {:?}", endpoint);

        let next_point = |point: PointId, edge: usize| {
            let (a, b) = (self.all[edge].points[0], self.all[edge].points[1]);
            if a == point {
                b
            } else {
                a
            }
        };

        let next_edges = |point: PointId, current_edge: usize| {
            let edges = map.get(&point).unwrap();
            let res: Vec<_> = edges.iter().filter(|e| **e != current_edge).collect();
            res
        };

        let mut current_edge = map.get(&endpoint).unwrap()[0];
        let mut current_point = next_point(endpoint, current_edge);
        println!("current_edge  = {:?}", self.all[current_edge].points);
        println!("current_point = {:?}", current_point);

        // follow path
        let mut path_points = Vec::with_capacity(self.all.len() + 1);
        let mut path_edges = Vec::with_capacity(self.all.len() + 1); // +1 for the loop case
        path_points.push(endpoint);
        path_points.push(current_point);
        path_edges.push(current_edge);
        loop {
            let edges = next_edges(current_point, current_edge);
            println!(
                "\nnext_edges    = {:?}",
                edges.iter().map(|e| self.all[**e].points.clone()).collect::<Vec<_>>()
            );
            if edges.len() == 0 {
                break;
            }
            current_edge = *edges[0];
            current_point = next_point(current_point, current_edge);
            if current_point == endpoint {
                path_edges.push(current_edge);
                break;
            }
            println!("current_edge  = {:?}", self.all[current_edge].points);
            println!("current_point = {:?}", current_point);
            path_points.push(current_point);
            path_edges.push(current_edge);
        }
        println!("\npath_points   = {:?}", path_points);
        println!("path_edges    = {:?}", path_edges);
        Ok((path_points, path_edges))
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
        // allocate edges
        let all = generate_sample();
        let empty = Edges { all: vec![] };
        let branching = Edges {
            all: vec![&all[5], &all[7], &all[6], &all[8]],
        };
        let disconnected = Edges {
            all: vec![&all[3], &all[4], &all[10]],
        };
        let loop1 = Edges {
            all: vec![&all[5], &all[9], &all[1], &all[7]], // 5 come first => clockwise loop
        };
        let loop2 = Edges {
            all: vec![&all[7], &all[9], &all[1], &all[5]], // 7 come first => counter-clockwise loop
        };
        let bottom = Edges {
            all: vec![&all[2], &all[1], &all[0]],
        };

        let (points, edges) = branching.any_path().unwrap();
        assert_eq!(points, &[1, 8, 2]);
        assert_eq!(edges, &[/*8*/ 3, /*5*/ 0]);

        println!("-----------------------------------");
        let (points, edges) = loop1.any_path().unwrap();
        assert_eq!(points, &[2, 8, 4, 3]);
        assert_eq!(edges, &[/*5*/ 0, /*7*/ 3, /*9*/ 1, /*1*/ 2]);

        println!("-----------------------------------");
        let (points, edges) = loop2.any_path().unwrap();
        assert_eq!(points, &[2, 3, 4, 8]);
        assert_eq!(edges, &[/*1*/ 2, /*9*/ 1, /*7*/ 0, /*5*/ 3]);
    }
}
