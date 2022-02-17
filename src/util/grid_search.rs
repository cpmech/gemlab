use crate::geometry::{point_circle_distance, point_cylinder_distance, point_line_distance, point_point_distance};
use crate::StrError;
use plotpy::{Curve, Plot, Shapes, Text};
use serde::{Deserialize, Serialize};
use std::collections::{HashMap, HashSet};
use std::fmt;

/// Holds the id and coordinates of an item
#[derive(Clone, Debug, Deserialize, Serialize)]
struct Item {
    id: usize,   // identification number
    x: Vec<f64>, // (ndim) coordinates
}

/// Holds items
#[derive(Clone, Debug, Deserialize, Serialize)]
struct Container {
    items: Vec<Item>,
}

/// Implements a grid for fast searching entries by coordinates
///
/// # Reference
///
/// * Durand, Farias, and Pedroso (2015) Computing intersections between
///   non-compatible curves and finite elements, Computational Mechanics;
///   DOI=10.1007/s00466-015-1181-y
#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct GridSearch {
    // constants
    ndim: usize,      // space dimension
    ndiv: Vec<usize>, // (ndim) number of divisions along each direction
    min: Vec<f64>,    // (ndim) min values
    max: Vec<f64>,    // (ndim) max values
    delta: Vec<f64>,  // (ndim) difference between max and min
    size: Vec<f64>,   // (ndim) side lengths of each container
    cf: Vec<usize>,   // (3) coefficients [1, ndiv[0], ndiv[0]*ndiv[1]] (Eq. 8)

    // square/cubic halo: bounding box corners around point, including the point
    tol: Vec<f64>,       // tolerances to compare coordinates and define the halo
    halo: Vec<Vec<f64>>, // (ncorner) 4 in 2D or 8 in 3D (each contains ndim coords)
    ncorner: usize,      // number of halo corners 4 in 2D or 8 in 3D

    // holds non-empty containers. maps container.index to container.data
    // a point may be located in more than one container (e.g., when at boundaries)
    containers: HashMap<usize, Container>,

    // indicates whether initialize has been called or not
    initialized: bool,

    // constants
    radius: f64,     // radius of the circumscribed circle of containers
    radius_tol: f64, // radius of the bounding box defined by the tolerances
}

impl GridSearch {
    // constants
    const DEFAULT_TOLERANCE: f64 = 1e-4; // default point-search tolerance for all directions

    /// Creates a new GridSearch
    pub fn new(ndim: usize) -> Result<Self, StrError> {
        // check
        if ndim < 2 || ndim > 3 {
            return Err("ndim must be 2 or 3");
        }

        // number of halo corners
        let ncorner = usize::pow(2, ndim as u32);

        // allocate grid
        let mut grid = GridSearch {
            ndim,
            ndiv: vec![0; ndim],
            min: vec![0.0; ndim],
            max: vec![0.0; ndim],
            delta: vec![0.0; ndim],
            size: vec![0.0; ndim],
            cf: vec![0; 3], // << must be 3
            tol: vec![0.0; ndim],
            halo: vec![vec![0.0; ndim]; ncorner],
            ncorner,
            containers: HashMap::new(),
            initialized: false,
            radius: 0.0,
            radius_tol: 0.0,
        };

        // set tolerances
        grid.set_tolerances(&vec![GridSearch::DEFAULT_TOLERANCE; ndim])?;

        // done
        Ok(grid)
    }

    /// Sets the tolerances to compare coordinates (len = ndim)
    ///
    /// # Note
    ///
    /// You must call `initialize` AFTER setting the tolerances.
    pub fn set_tolerances(&mut self, tol: &[f64]) -> Result<(), StrError> {
        if tol.len() != self.ndim {
            return Err("tol.len() must equal ndim");
        }
        self.tol.copy_from_slice(tol);
        self.radius_tol = 0.0;
        for i in 0..self.ndim {
            self.radius_tol += self.tol[i] * self.tol[i];
        }
        self.radius_tol = f64::sqrt(self.radius_tol);
        self.initialized = false;
        Ok(())
    }

    /// Initializes the GridSearch
    ///
    /// # Input
    ///
    /// * `ndiv` -- number of divisions along each dimension (len = 2 or 3 == ndim)
    /// * `min` -- min coordinates (len = ndim)
    /// * `max` -- max coordinates (len = ndim)
    pub fn initialize(&mut self, ndiv: &[usize], min: &[f64], max: &[f64]) -> Result<(), StrError> {
        // check input
        self.initialized = false;
        if ndiv.len() != self.ndim {
            return Err("ndiv.len() must equal ndim");
        }
        if min.len() != self.ndim {
            return Err("min.len() must equal ndim");
        }
        if max.len() != self.ndim {
            return Err("max.len() must equal ndim");
        }

        // check and compute sizes
        self.radius = 0.0;
        for i in 0..self.ndim {
            self.ndiv[i] = ndiv[i];
            if ndiv[i] < 1 {
                return Err("ndiv must be at least 1");
            }
            self.min[i] = min[i];
            self.max[i] = max[i];
            self.delta[i] = max[i] - min[i];
            if self.delta[i] <= 0.0 {
                return Err("max must be greater than min");
            }
            self.size[i] = self.delta[i] / (ndiv[i] as f64);
            self.tol[i] = self.tol[i];
            if self.size[i] <= 10.0 * self.tol[i] {
                return Err("ndiv is too high; ndiv must be such that the grid size is greater than 10 * tol");
            }
            self.radius += self.size[i] * self.size[i] / 4.0;
        }
        self.radius = f64::sqrt(self.radius);

        // coefficient
        self.cf[0] = 1;
        self.cf[1] = self.ndiv[0];
        self.cf[2] = self.ndiv[0] * self.ndiv[1];

        // done
        self.initialized = true;
        Ok(())
    }

    /// Inserts a new item to the right container in the grid
    ///
    /// # Input
    ///
    /// * `id` -- identification number for the item
    /// * `x` -- coordinates (ndim) of the item
    pub fn insert(&mut self, id: usize, x: &[f64]) -> Result<(), StrError> {
        // check
        if !self.initialized {
            return Err("initialize must be called first");
        }
        if x.len() != self.ndim {
            return Err("x.len() must equal ndim");
        }

        // add point to container
        let index = match self.container_index(x) {
            Some(i) => i,
            None => return Err("point is outside the grid"),
        };
        self.update_or_insert(index, id, x);

        // add point to containers touched by halo corners
        self.set_halo(x);
        let mut tmp = vec![0.0; self.ndim];
        for c in 0..self.ncorner {
            tmp.copy_from_slice(&self.halo[c][0..self.ndim]);
            if let Some(index_corner) = self.container_index(&tmp) {
                if index_corner != index {
                    self.update_or_insert(index_corner, id, x); // make sure to use original `x`
                }
            }
        }
        Ok(())
    }

    /// Find previously inserted item to the grid
    ///
    /// # Input
    ///
    /// * `x` -- coordinates (ndim) of the item
    ///
    /// # Output
    ///
    /// * `id` -- if found, returns the identification number of the item
    pub fn find(&self, x: &[f64]) -> Result<Option<usize>, StrError> {
        // check
        if !self.initialized {
            return Err("initialize must be called first");
        }
        if x.len() != self.ndim {
            return Err("x.len() must equal ndim");
        }

        // find index of container where x should be
        let index = match self.container_index(x) {
            Some(i) => i,
            None => return Err("point is outside the grid"),
        };

        // find container that should have a point close to `x`
        let container = match self.containers.get(&index) {
            Some(c) => c,
            None => return Ok(None), // no container has a point close to `x`
        };

        // find closest point to `x` in the container
        for item in &container.items {
            let distance = point_point_distance(&item.x, x)?;
            if distance <= self.radius_tol {
                return Ok(Some(item.id));
            }
        }
        Ok(None)
    }

    /// Finds points on a 2D or 3D line
    ///
    /// # Input
    ///
    /// * `a` -- (ndim) first point on the line
    /// * `b` -- (ndim) second point on the line (different than `a`)
    ///
    /// # Output
    ///
    /// Returns the ids of points.
    pub fn find_on_line(&self, a: &[f64], b: &[f64]) -> Result<HashSet<usize>, StrError> {
        // check
        if !self.initialized {
            return Err("initialize must be called first");
        }
        if a.len() != self.ndim {
            return Err("a.len() must equal ndim");
        }
        if b.len() != self.ndim {
            return Err("b.len() must equal ndim");
        }

        // find containers near the line
        let nearest_containers = self.containers_near_line(a, b)?;

        // find container points near the line
        let mut ids = HashSet::new();
        for index in nearest_containers {
            let container = self.containers.get(&index).unwrap();
            for item in &container.items {
                let distance = point_line_distance(a, b, &item.x)?;
                if distance <= self.radius_tol {
                    ids.insert(item.id.clone());
                }
            }
        }
        Ok(ids)
    }

    /// Finds points on the perimeter of a circle (2D only)
    ///
    /// # Input
    ///
    /// `center` -- 2D circle center
    /// `radius` -- circle radius
    ///
    /// # Output
    ///
    /// Returns the ids of points.
    ///
    /// # Note
    ///
    /// This works in 2D only.
    pub fn find_on_circle(&self, center: &[f64], radius: f64) -> Result<HashSet<usize>, StrError> {
        // check
        if !self.initialized {
            return Err("initialize must be called first");
        }
        if self.ndim != 2 {
            return Err("this works in 2D only");
        }
        if center.len() != self.ndim {
            return Err("center.len() must equal ndim");
        }

        // find containers near the circle
        let nearest_containers = self.containers_near_circle(center, radius)?;

        // find container points near the circle
        let mut ids = HashSet::new();
        for index in nearest_containers {
            let container = self.containers.get(&index).unwrap();
            for item in &container.items {
                let distance = point_circle_distance(center, radius, &item.x)?;
                if f64::abs(distance) <= self.radius_tol {
                    ids.insert(item.id.clone());
                }
            }
        }
        Ok(ids)
    }

    /// Finds points on the surface of a cylinder (3D only)
    ///
    /// # Input
    ///
    /// `a` -- 3D point on the cylinder axis
    /// `b` -- 3D point on the cylinder axis
    /// `radius` -- cylinder radius
    ///
    /// # Output
    ///
    /// Returns the ids of points.
    ///
    /// # Note
    ///
    /// This works in 3D only.
    pub fn find_on_cylinder(&self, a: &[f64], b: &[f64], radius: f64) -> Result<HashSet<usize>, StrError> {
        // check
        if !self.initialized {
            return Err("initialize must be called first");
        }
        if self.ndim != 3 {
            return Err("this works in 3D only");
        }
        if a.len() != self.ndim {
            return Err("a.len() must equal ndim");
        }
        if b.len() != self.ndim {
            return Err("b.len() must equal ndim");
        }

        // find containers near the cylinder
        let nearest_containers = self.containers_near_cylinder(a, b, radius)?;

        // find container points near the cylinder
        let mut ids = HashSet::new();
        for index in nearest_containers {
            let container = self.containers.get(&index).unwrap();
            for item in &container.items {
                let distance = point_cylinder_distance(a, b, radius, &item.x)?;
                if f64::abs(distance) <= self.radius_tol {
                    ids.insert(item.id.clone());
                }
            }
        }
        Ok(ids)
    }

    /// Find points on the x-y plane (3D only)
    ///
    /// # Input
    ///
    /// `z` -- the plane passes through `z`
    ///
    /// # Output
    ///
    /// Returns the ids of points.
    ///
    /// # Note
    ///
    /// This works in 3D only.
    pub fn find_on_plane_xy(&self, z: f64) -> Result<HashSet<usize>, StrError> {
        // check
        if !self.initialized {
            return Err("initialize must be called first");
        }
        if self.ndim != 3 {
            return Err("this works in 3D only");
        }

        // find containers near the plane
        let nearest_containers = self.containers_near_plane(2, z);

        // find container points near the plane
        let mut ids = HashSet::new();
        for index in nearest_containers {
            let container = self.containers.get(&index).unwrap();
            for item in &container.items {
                let distance = f64::abs(item.x[2] - z);
                if f64::abs(distance) <= self.radius_tol {
                    ids.insert(item.id.clone());
                }
            }
        }
        Ok(ids)
    }

    /// Find points on the y-z plane (3D only)
    ///
    /// # Input
    ///
    /// `x` -- the plane passes through `x`
    ///
    /// # Output
    ///
    /// Returns the ids of points.
    ///
    /// # Note
    ///
    /// This works in 3D only.
    pub fn find_on_plane_yz(&self, x: f64) -> Result<HashSet<usize>, StrError> {
        // check
        if !self.initialized {
            return Err("initialize must be called first");
        }
        if self.ndim != 3 {
            return Err("this works in 3D only");
        }

        // find containers near the plane
        let nearest_containers = self.containers_near_plane(0, x);

        // find container points near the plane
        let mut ids = HashSet::new();
        for index in nearest_containers {
            let container = self.containers.get(&index).unwrap();
            for item in &container.items {
                let distance = f64::abs(item.x[0] - x);
                if f64::abs(distance) <= self.radius_tol {
                    ids.insert(item.id.clone());
                }
            }
        }
        Ok(ids)
    }

    /// Find points on the x-z plane (3D only)
    ///
    /// # Input
    ///
    /// `y` -- the plane passes through `y`
    ///
    /// # Output
    ///
    /// Returns the ids of points.
    ///
    /// # Note
    ///
    /// This works in 3D only.
    pub fn find_on_plane_xz(&self, y: f64) -> Result<HashSet<usize>, StrError> {
        // check
        if !self.initialized {
            return Err("initialize must be called first");
        }
        if self.ndim != 3 {
            return Err("this works in 3D only");
        }

        // find containers near the plane
        let nearest_containers = self.containers_near_plane(1, y);

        // find container points near the plane
        let mut ids = HashSet::new();
        for index in nearest_containers {
            let container = self.containers.get(&index).unwrap();
            for item in &container.items {
                let distance = f64::abs(item.x[1] - y);
                if f64::abs(distance) <= self.radius_tol {
                    ids.insert(item.id.clone());
                }
            }
        }
        Ok(ids)
    }

    /// Inserts a new item if the coordinates are not found
    ///
    /// # Input
    ///
    /// * `id` -- identification number for the new item
    /// * `x` -- coordinates (ndim) of the item
    ///
    /// # Output
    ///
    /// * Returns the id of the item; either found or given as input.
    /// * Thus, if the output and input ids are different, it means that
    ///   the item exists already and the output id is the existing id.
    ///
    /// # Note
    ///
    /// This function is a combination of `find` and `insert`
    pub fn maybe_insert(&mut self, id: usize, x: &[f64]) -> Result<usize, StrError> {
        if let Some(id) = self.find(x)? {
            return Ok(id);
        }
        self.insert(id, x)?;
        Ok(id)
    }

    /// Returns a drawing of this object
    pub fn plot(&self) -> Result<Plot, StrError> {
        // check
        if !self.initialized {
            return Err("initialize must be called first");
        }

        // create plot
        let mut plot = Plot::new();

        // draw grid
        let mut xmin = vec![0.0; self.ndim];
        let mut xmax = vec![0.0; self.ndim];
        let mut ndiv = vec![0; self.ndim];
        for i in 0..self.ndim {
            xmin[i] = self.min[i];
            xmax[i] = self.max[i];
            ndiv[i] = self.ndiv[i];
        }
        let mut shapes = Shapes::new();
        shapes
            .set_alt_text_color("#5d5d5d")
            .draw_grid(&xmin, &xmax, &ndiv, false, true)?;
        plot.add(&shapes);

        // draw items
        let mut curve = Curve::new();
        let mut text = Text::new();
        curve
            .set_marker_style("o")
            .set_marker_color("#fab32faa")
            .set_marker_line_color("black")
            .set_marker_line_width(0.5);
        text.set_color("#cd0000");
        for container in self.containers.values() {
            for item in &container.items {
                let txt = format!("{}", item.id);
                if self.ndim == 2 {
                    curve.draw(&[item.x[0]], &[item.x[1]]);
                    text.draw(item.x[0], item.x[1], &txt);
                } else {
                    curve.draw_3d(&[item.x[0]], &[item.x[1]], &[item.x[2]]);
                    text.draw_3d(item.x[0], item.x[1], item.x[2], &txt);
                }
            }
        }
        plot.add(&curve);
        plot.add(&text);

        // done
        Ok(plot)
    }

    /// Calculates the container index where the point x should be located
    ///
    /// # Output
    ///
    /// * returns the index of the container or None if the point is out-of-range
    #[inline]
    fn container_index(&self, x: &[f64]) -> Option<usize> {
        let mut ratio = vec![0; self.ndim]; // (ndim) ratio = trunc(δx[i]/Δx[i]) (Eq. 8)
        let mut index = 0;
        for i in 0..self.ndim {
            if x[i] < self.min[i] || x[i] > self.max[i] {
                return None;
            }
            ratio[i] = ((x[i] - self.min[i]) / self.size[i]) as usize;
            if ratio[i] == self.ndiv[i] {
                // the point is exactly on the max edge, thus select inner container
                ratio[i] -= 1; // move to the inside
            }
            index += ratio[i] * self.cf[i];
        }
        Some(index)
    }

    /// Computes the i,j,k indices of the lower-left-bottom corner of the container
    #[inline]
    fn container_pivot_indices(&self, index: usize) -> (usize, usize, usize) {
        let i = index % self.ndiv[0];
        let j = (index % self.cf[2]) / self.ndiv[0];
        let k = index / self.cf[2];
        (i, j, k)
    }

    /// Computes the center coordinates of container
    #[inline]
    fn container_center(&self, cen: &mut [f64], i: usize, j: usize, k: usize) {
        cen[0] = self.min[0] + (i as f64) * self.size[0] + self.size[0] / 2.0;
        cen[1] = self.min[1] + (j as f64) * self.size[1] + self.size[1] / 2.0;
        if self.ndim == 3 {
            cen[2] = self.min[2] + (k as f64) * self.size[2] + self.size[2] / 2.0;
        }
    }

    /// Returns the indices of containers near a line
    #[inline]
    fn containers_near_line(&self, a: &[f64], b: &[f64]) -> Result<Vec<usize>, StrError> {
        let mut nearest_containers = Vec::new();
        let mut cen = vec![0.0; self.ndim];
        for index in self.containers.keys() {
            // compute container center
            let (i, j, k) = self.container_pivot_indices(*index);
            self.container_center(&mut cen, i, j, k);
            // check if the center of container is near the segment
            let distance = point_line_distance(a, b, &cen)?;
            if distance <= self.radius + self.radius_tol {
                nearest_containers.push(*index);
            }
        }
        Ok(nearest_containers)
    }

    /// Returns the indices of containers near a circle
    #[inline]
    fn containers_near_circle(&self, circle_center: &[f64], radius: f64) -> Result<Vec<usize>, StrError> {
        let mut nearest_containers = Vec::new();
        let mut cen = vec![0.0; self.ndim];
        for index in self.containers.keys() {
            // compute container center
            let (i, j, k) = self.container_pivot_indices(*index);
            self.container_center(&mut cen, i, j, k);
            // check if the center of container is near the circle
            let distance = point_circle_distance(circle_center, radius, &cen)?;
            if distance <= self.radius + self.radius_tol {
                nearest_containers.push(*index);
            }
        }
        Ok(nearest_containers)
    }

    /// Returns the indices of containers near a circle
    #[inline]
    fn containers_near_cylinder(&self, a: &[f64], b: &[f64], radius: f64) -> Result<Vec<usize>, StrError> {
        let mut nearest_containers = Vec::new();
        let mut cen = vec![0.0; self.ndim];
        for index in self.containers.keys() {
            // compute container center
            let (i, j, k) = self.container_pivot_indices(*index);
            self.container_center(&mut cen, i, j, k);
            // check if the center of container is near the cylinder
            let distance = point_cylinder_distance(a, b, radius, &cen)?;
            if distance <= self.radius + self.radius_tol {
                nearest_containers.push(*index);
            }
        }
        Ok(nearest_containers)
    }

    /// Returns the indices of containers near plane xy
    #[inline]
    fn containers_near_plane(&self, fixed_dim: usize, fixed_coord: f64) -> Vec<usize> {
        let mut nearest_containers = Vec::new();
        let mut cen = vec![0.0; self.ndim];
        for index in self.containers.keys() {
            // compute container center
            let (i, j, k) = self.container_pivot_indices(*index);
            self.container_center(&mut cen, i, j, k);
            // check if the center of container is near the plane
            let distance = f64::abs(cen[fixed_dim] - fixed_coord);
            if distance <= self.radius + self.radius_tol {
                nearest_containers.push(*index);
            }
        }
        nearest_containers
    }

    /// Update container or insert point in container
    #[inline]
    fn update_or_insert(&mut self, index: usize, id: usize, x: &[f64]) {
        let item = Item { id, x: Vec::from(x) };
        if self.containers.contains_key(&index) {
            let container = self.containers.get_mut(&index).unwrap();
            container.items.push(item);
        } else {
            self.containers.insert(index, Container { items: vec![item] });
        }
    }

    /// Sets square/cubic halo around point
    #[inline]
    fn set_halo(&mut self, x: &[f64]) {
        if self.ndim == 2 {
            self.halo[0][0] = x[0] - self.tol[0];
            self.halo[0][1] = x[1] - self.tol[1];

            self.halo[1][0] = x[0] + self.tol[0];
            self.halo[1][1] = x[1] - self.tol[1];

            self.halo[2][0] = x[0] + self.tol[0];
            self.halo[2][1] = x[1] + self.tol[1];

            self.halo[3][0] = x[0] - self.tol[0];
            self.halo[3][1] = x[1] + self.tol[1];
        } else {
            self.halo[0][0] = x[0] - self.tol[0];
            self.halo[0][1] = x[1] - self.tol[1];
            self.halo[0][2] = x[2] - self.tol[2];

            self.halo[1][0] = x[0] + self.tol[0];
            self.halo[1][1] = x[1] - self.tol[1];
            self.halo[1][2] = x[2] - self.tol[2];

            self.halo[2][0] = x[0] + self.tol[0];
            self.halo[2][1] = x[1] + self.tol[1];
            self.halo[2][2] = x[2] - self.tol[2];

            self.halo[3][0] = x[0] - self.tol[0];
            self.halo[3][1] = x[1] + self.tol[1];
            self.halo[3][2] = x[2] - self.tol[2];

            self.halo[4][0] = x[0] - self.tol[0];
            self.halo[4][1] = x[1] - self.tol[1];
            self.halo[4][2] = x[2] + self.tol[2];

            self.halo[5][0] = x[0] + self.tol[0];
            self.halo[5][1] = x[1] - self.tol[1];
            self.halo[5][2] = x[2] + self.tol[2];

            self.halo[6][0] = x[0] + self.tol[0];
            self.halo[6][1] = x[1] + self.tol[1];
            self.halo[6][2] = x[2] + self.tol[2];

            self.halo[7][0] = x[0] - self.tol[0];
            self.halo[7][1] = x[1] + self.tol[1];
            self.halo[7][2] = x[2] + self.tol[2];
        }
    }
}

impl fmt::Display for GridSearch {
    /// Shows info about the items in the grid containers
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        // items
        let mut unique_items: HashMap<usize, bool> = HashMap::new();
        let mut indices: Vec<_> = self.containers.keys().collect();
        indices.sort();
        for index in indices {
            let container = self.containers.get(index).unwrap();
            let mut ids: Vec<_> = container.items.iter().map(|item| item.id).collect();
            ids.sort();
            write!(f, "{}: {:?}\n", index, ids).unwrap();
            for id in ids {
                unique_items.insert(id, true);
            }
        }
        // summary
        let mut ids: Vec<_> = unique_items.keys().collect();
        ids.sort();
        write!(f, "ids = {:?}\n", ids).unwrap();
        write!(f, "nitem = {}\n", unique_items.len()).unwrap();
        write!(f, "ncontainer = {}\n", self.containers.len()).unwrap();
        write!(f, "ndiv = {:?}\n", self.ndiv).unwrap();
        Ok(())
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::*;
    use crate::util::SQRT_2;
    use plotpy::{Curve, Shapes, Surface};
    use russell_chk::{assert_approx_eq, assert_vec_approx_eq};

    struct TestData<'a> {
        id: usize,
        x: &'a [f64],
        container: usize,        // container index calc with x only
        containers: &'a [usize], // container index calc with halo and x
    }

    fn get_test_grid_2d() -> GridSearch {
        let mut grid = GridSearch::new(2).unwrap();
        grid.initialize(&[5, 5], &[-0.2, -0.2], &[0.8, 1.8]).unwrap();
        grid
    }

    fn get_test_grid_3d() -> GridSearch {
        let mut grid = GridSearch::new(3).unwrap();
        grid.initialize(&[3, 3, 3], &[-1.0, -1.0, -1.0], &[1.0, 1.0, 1.0])
            .unwrap();
        grid
    }

    fn get_test_data_2d<'a>() -> Vec<TestData<'a>> {
        vec![
            TestData {
                id: 100,
                x: &[0.0, 0.2],
                container: 6,
                containers: &[0, 1, 5, 6],
            },
            TestData {
                id: 200,
                x: &[0.2, 0.6000000000000001],
                container: 12,
                containers: &[6, 7, 11, 12],
            },
            TestData {
                id: 300,
                x: &[0.4, 1.0],
                // in this case, the ratio is:
                // 3.0000000000000004, 2.9999999999999996
                // thus, the point falls in #13 instead of #18
                container: 13,
                containers: &[12, 13, 17, 18],
            },
            TestData {
                id: 400,
                x: &[0.6, 1.4],
                // in this case, the ratio is:
                // 4, 3.9999999999999996
                // thus, the point falls in #19 instead of #24
                container: 19,
                containers: &[18, 19, 23, 24],
            },
            TestData {
                id: 500,
                x: &[0.8, 1.8],
                container: 24,
                containers: &[24],
            },
            TestData {
                id: 600,
                x: &[0.6, 0.0],
                container: 4,
                containers: &[3, 4], // will be on 3 twice
            },
            // circle @ (-0.2,1.8) radius 0.3
            TestData {
                id: 101,
                x: &[-0.2, 1.5],
                container: 20,
                containers: &[20],
            },
            TestData {
                id: 102,
                x: &[-0.2 + 0.3 * SQRT_2 / 2.0, 1.8 - 0.3 * SQRT_2 / 2.0],
                container: 21,
                containers: &[21],
            },
            TestData {
                id: 103,
                x: &[0.1, 1.8],
                container: 21,
                containers: &[21],
            },
        ]
    }

    fn get_test_data_3d<'a>() -> Vec<TestData<'a>> {
        const L: f64 = 2.0;
        vec![
            TestData {
                id: 100,
                x: &[-1.0, -1.0, -1.0],
                container: 0,
                containers: &[0],
            },
            TestData {
                id: 200,
                x: &[0.0, 0.0, 0.0],
                container: 13,
                containers: &[13],
            },
            TestData {
                id: 300,
                x: &[-1.0 + L * 2.0 / 3.0, -1.0 + L * 2.0 / 3.0, -1.0 + L * 2.0 / 3.0],
                container: 26,
                containers: &[13, 14, 16, 17, 22, 23, 25, 26],
            },
            TestData {
                id: 400,
                x: &[-1.0 + L * 2.0 / 3.0, -0.7, -0.7],
                container: 2,
                containers: &[1, 2], // will be on 1 four times
            },
            TestData {
                id: 500,
                x: &[0.75, -0.75, -0.75],
                container: 2,
                containers: &[2],
            },
            // cylinder with 0.4 radius parallel to y passing through (1,-1,-1) and (1,1,-1)
            TestData {
                id: 101,
                x: &[1.0 - 0.4, -1.0, -1.0],
                container: 2,
                containers: &[2],
            },
            TestData {
                id: 102,
                x: &[1.0 - 0.4 * SQRT_2 / 2.0, -1.0, -1.0 + 0.4 * SQRT_2 / 2.0],
                container: 2,
                containers: &[2],
            },
            TestData {
                id: 103,
                x: &[1.0, -1.0, -1.0 + 0.4],
                container: 2,
                containers: &[2],
            },
            TestData {
                id: 104,
                x: &[1.0 - 0.4, 0.75, -1.0],
                container: 8,
                containers: &[8],
            },
            TestData {
                id: 105,
                x: &[1.0 - 0.4 * SQRT_2 / 2.0, 0.75, -1.0 + 0.4 * SQRT_2 / 2.0],
                container: 8,
                containers: &[8],
            },
            TestData {
                id: 106,
                x: &[1.0, 0.75, -1.0 + 0.4],
                container: 8,
                containers: &[8],
            },
        ]
    }

    #[test]
    fn new_fails_on_wrong_input() {
        let grid = GridSearch::new(1);
        assert_eq!(grid.err(), Some("ndim must be 2 or 3"));
    }

    #[test]
    fn set_tolerances_fails_on_wrong_input() -> Result<(), StrError> {
        let mut grid = GridSearch::new(2)?;
        let res = grid.set_tolerances(&[0.1]);
        assert_eq!(res.err(), Some("tol.len() must equal ndim"));
        Ok(())
    }

    #[test]
    fn set_tolerances_works() -> Result<(), StrError> {
        let mut grid = GridSearch::new(2)?;
        grid.set_tolerances(&[0.2, 0.2])?;
        assert_eq!(grid.tol, &[0.2, 0.2]);
        assert_eq!(grid.radius_tol, f64::sqrt(0.2 * 0.2 + 0.2 * 0.2));
        Ok(())
    }

    #[test]
    fn initialize_fails_on_wrong_input() -> Result<(), StrError> {
        let mut grid = GridSearch::new(2)?;
        let res = grid.initialize(&[1], &[0.0, 0.0], &[1.0, 1.0]);
        assert_eq!(res.err(), Some("ndiv.len() must equal ndim"));

        let mut grid = GridSearch::new(2)?;
        let res = grid.initialize(&[1, 1], &[0.0], &[1.0, 1.0]);
        assert_eq!(res.err(), Some("min.len() must equal ndim"));

        let mut grid = GridSearch::new(2)?;
        let res = grid.initialize(&[1, 1], &[0.0, 0.0], &[1.0]);
        assert_eq!(res.err(), Some("max.len() must equal ndim"));

        let mut grid = GridSearch::new(2)?;
        let res = grid.initialize(&[0, 1], &[0.0, 0.0], &[1.0, 1.0]);
        assert_eq!(res.err(), Some("ndiv must be at least 1"));

        let mut grid = GridSearch::new(2)?;
        let res = grid.initialize(&[1, 1], &[0.0, 0.0], &[0.0, 1.0]);
        assert_eq!(res.err(), Some("max must be greater than min"));

        let mut grid = GridSearch::new(2)?;
        grid.set_tolerances(&[0.1, 0.1])?;
        let res = grid.initialize(&[10, 10], &[0.0, 0.0], &[1.0, 1.0]);
        assert_eq!(
            res.err(),
            Some("ndiv is too high; ndiv must be such that the grid size is greater than 10 * tol")
        );
        Ok(())
    }

    #[test]
    fn new_works() {
        let g2d = get_test_grid_2d();
        assert_eq!(g2d.ndim, 2);
        assert_eq!(g2d.ndiv, [5, 5]);
        assert_eq!(g2d.min, [-0.2, -0.2]);
        assert_eq!(g2d.max, [0.8, 1.8]);
        assert_eq!(g2d.delta, [1.0, 2.0]);
        assert_eq!(g2d.size, [0.2, 0.4]);
        assert_eq!(g2d.cf, [1, 5, 25]);
        assert_eq!(g2d.tol, [1e-4, 1e-4]);
        assert_eq!(g2d.halo.len(), 4);
        assert_eq!(g2d.ncorner, 4);
        assert_eq!(g2d.containers.len(), 0);
        assert_eq!(g2d.initialized, true);
        assert_approx_eq!(g2d.radius, f64::sqrt(0.1 * 0.1 + 0.2 * 0.2), 1e-15);
        assert_approx_eq!(g2d.radius_tol, f64::sqrt(2.0 * 1e-4 * 1e-4), 1e-15);

        let g3d = get_test_grid_3d();
        assert_eq!(g3d.ndim, 3);
        assert_eq!(g3d.ndiv, [3, 3, 3]);
        assert_eq!(g3d.min, [-1.0, -1.0, -1.0]);
        assert_eq!(g3d.max, [1.0, 1.0, 1.0]);
        assert_eq!(g3d.delta, [2.0, 2.0, 2.0]);
        assert_eq!(g3d.size, [2.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0]);
        assert_eq!(g3d.cf, [1, 3, 9]);
        assert_eq!(g3d.tol, [1e-4, 1e-4, 1e-4]);
        assert_eq!(g3d.halo.len(), 8);
        assert_eq!(g3d.ncorner, 8);
        assert_eq!(g3d.containers.len(), 0);
        assert_eq!(g2d.initialized, true);
        assert_approx_eq!(g3d.radius, f64::sqrt(3.0 * (1.0 / 3.0) * (1.0 / 3.0)), 1e-15);
        assert_approx_eq!(g3d.radius_tol, f64::sqrt(3.0 * 1e-4 * 1e-4), 1e-15);
    }

    #[test]
    fn display_trait_works() {
        let g2d = get_test_grid_2d();
        assert_eq!(
            format!("{}", g2d),
            "ids = []\n\
             nitem = 0\n\
             ncontainer = 0\n\
             ndiv = [5, 5]\n"
        );

        let g3d = get_test_grid_3d();
        assert_eq!(
            format!("{}", g3d),
            "ids = []\n\
             nitem = 0\n\
             ncontainer = 0\n\
             ndiv = [3, 3, 3]\n"
        );
    }

    #[test]
    fn set_halo_works() {
        let mut g2d = get_test_grid_2d();
        g2d.set_halo(&[0.5, 0.5]);
        assert_eq!(g2d.halo[0], [0.4999, 0.4999]);
        assert_eq!(g2d.halo[1], [0.5001, 0.4999]);
        assert_eq!(g2d.halo[2], [0.5001, 0.5001]);
        assert_eq!(g2d.halo[3], [0.4999, 0.5001]);

        let mut g3d = get_test_grid_3d();
        g3d.set_halo(&[0.5, 0.5, 0.5]);
        assert_eq!(g3d.halo[0], [0.4999, 0.4999, 0.4999]);
        assert_eq!(g3d.halo[1], [0.5001, 0.4999, 0.4999]);
        assert_eq!(g3d.halo[2], [0.5001, 0.5001, 0.4999]);
        assert_eq!(g3d.halo[3], [0.4999, 0.5001, 0.4999]);
        assert_eq!(g3d.halo[4], [0.4999, 0.4999, 0.5001]);
        assert_eq!(g3d.halo[5], [0.5001, 0.4999, 0.5001]);
        assert_eq!(g3d.halo[6], [0.5001, 0.5001, 0.5001]);
        assert_eq!(g3d.halo[7], [0.4999, 0.5001, 0.5001]);
    }

    #[test]
    fn container_index_works() -> Result<(), StrError> {
        let g2d = get_test_grid_2d();
        for data in get_test_data_2d() {
            let index = g2d.container_index(data.x).unwrap();
            assert_eq!(index, data.container);
        }
        let index = g2d.container_index(&[0.80001, 0.0]);
        assert_eq!(index, None); // outside

        let g3d = get_test_grid_3d();
        for data in get_test_data_3d() {
            let index = g3d.container_index(data.x).unwrap();
            assert_eq!(index, data.container);
        }
        let index = g3d.container_index(&[1.00001, 0.0, 0.0]);
        assert_eq!(index, None); // outside
        Ok(())
    }

    #[test]
    fn container_pivot_indices_works() -> Result<(), StrError> {
        let g2d = get_test_grid_2d();
        assert_eq!(g2d.container_pivot_indices(0), (0, 0, 0));
        assert_eq!(g2d.container_pivot_indices(4), (4, 0, 0));
        assert_eq!(g2d.container_pivot_indices(20), (0, 4, 0));
        assert_eq!(g2d.container_pivot_indices(24), (4, 4, 0));
        assert_eq!(g2d.container_pivot_indices(8), (3, 1, 0));
        assert_eq!(g2d.container_pivot_indices(17), (2, 3, 0));
        assert_eq!(g2d.container_pivot_indices(23), (3, 4, 0));

        let g3d = get_test_grid_3d();
        assert_eq!(g3d.container_pivot_indices(0), (0, 0, 0));
        assert_eq!(g3d.container_pivot_indices(2), (2, 0, 0));
        assert_eq!(g3d.container_pivot_indices(6), (0, 2, 0));
        assert_eq!(g3d.container_pivot_indices(8), (2, 2, 0));
        assert_eq!(g3d.container_pivot_indices(18), (0, 0, 2));
        assert_eq!(g3d.container_pivot_indices(20), (2, 0, 2));
        assert_eq!(g3d.container_pivot_indices(24), (0, 2, 2));
        assert_eq!(g3d.container_pivot_indices(26), (2, 2, 2));
        assert_eq!(g3d.container_pivot_indices(13), (1, 1, 1));
        Ok(())
    }

    #[test]
    fn container_center_works() -> Result<(), StrError> {
        let g2d = get_test_grid_2d();
        let mut cen2d = vec![0.0; 2];
        g2d.container_center(&mut cen2d, 0, 0, 0);
        assert_vec_approx_eq!(cen2d, &[-0.1, 0.0], 1e-15);
        g2d.container_center(&mut cen2d, 4, 0, 0);
        assert_vec_approx_eq!(cen2d, &[0.7, 0.0], 1e-15);
        g2d.container_center(&mut cen2d, 0, 4, 0);
        assert_vec_approx_eq!(cen2d, &[-0.1, 1.6], 1e-15);
        g2d.container_center(&mut cen2d, 4, 4, 0);
        assert_vec_approx_eq!(cen2d, &[0.7, 1.6], 1e-15);
        g2d.container_center(&mut cen2d, 3, 1, 0);
        assert_vec_approx_eq!(cen2d, &[0.5, 0.4], 1e-15);
        g2d.container_center(&mut cen2d, 2, 3, 0);
        assert_vec_approx_eq!(cen2d, &[0.3, 1.2], 1e-15);
        g2d.container_center(&mut cen2d, 3, 4, 0);
        assert_vec_approx_eq!(cen2d, &[0.5, 1.6], 1e-15);

        let g3d = get_test_grid_3d();
        let mut cen3d = vec![0.0; 3];
        let l = 2.0 / 3.0;
        let hl = l / 2.0;
        let l2 = l * 2.0;
        g3d.container_center(&mut cen3d, 0, 0, 0);
        assert_vec_approx_eq!(cen3d, &[-1.0 + hl, -1.0 + hl, -1.0 + hl], 1e-14);
        g3d.container_center(&mut cen3d, 2, 0, 0);
        assert_vec_approx_eq!(cen3d, &[-1.0 + l2 + hl, -1.0 + hl, -1.0 + hl], 1e-14);
        g3d.container_center(&mut cen3d, 0, 2, 0);
        assert_vec_approx_eq!(cen3d, &[-1.0 + hl, -1.0 + l2 + hl, -1.0 + hl], 1e-14);
        g3d.container_center(&mut cen3d, 2, 2, 2);
        assert_vec_approx_eq!(cen3d, &[-1.0 + l2 + hl, -1.0 + l2 + hl, -1.0 + l2 + hl], 1e-14);
        Ok(())
    }

    #[test]
    fn insert_fails_on_wrong_input() {
        let mut g2d = get_test_grid_2d();
        let res = g2d.insert(0, &[0.0, 0.0, 0.0]);
        assert_eq!(res, Err("x.len() must equal ndim"));
        let res = g2d.insert(1000, &[0.80001, 0.0]);
        assert_eq!(res, Err("point is outside the grid"));

        let mut g3d = get_test_grid_3d();
        let res = g3d.insert(0, &[0.0, 0.0]);
        assert_eq!(res, Err("x.len() must equal ndim"));
        let res = g3d.insert(1000, &[1.00001, 0.0, 0.0]);
        assert_eq!(res, Err("point is outside the grid"));
    }

    #[test]
    fn insert_2d_works() -> Result<(), StrError> {
        let mut grid = get_test_grid_2d();
        for data in get_test_data_2d() {
            grid.insert(data.id, data.x)?;
            for index in data.containers {
                let container = grid.containers.get(index).unwrap();
                container.items.iter().find(|item| item.id == data.id).unwrap();
            }
        }
        assert_eq!(
            format!("{}", grid),
            "0: [100]\n\
             1: [100]\n\
             3: [600, 600]\n\
             4: [600]\n\
             5: [100]\n\
             6: [100, 200]\n\
             7: [200]\n\
             11: [200]\n\
             12: [200, 300]\n\
             13: [300]\n\
             17: [300]\n\
             18: [300, 400]\n\
             19: [400]\n\
             20: [101]\n\
             21: [102, 103]\n\
             23: [400]\n\
             24: [400, 500]\n\
             ids = [100, 101, 102, 103, 200, 300, 400, 500, 600]\n\
             nitem = 9\n\
             ncontainer = 17\n\
             ndiv = [5, 5]\n"
        );
        let mut indices: Vec<_> = grid.containers.into_keys().collect();
        indices.sort();
        assert_eq!(indices, &[0, 1, 3, 4, 5, 6, 7, 11, 12, 13, 17, 18, 19, 20, 21, 23, 24]);
        Ok(())
    }

    #[test]
    fn insert_3d_works() -> Result<(), StrError> {
        let mut grid = get_test_grid_3d();
        for data in get_test_data_3d() {
            grid.insert(data.id, data.x)?;
            for index in data.containers {
                let container = grid.containers.get(index).unwrap();
                container.items.iter().find(|item| item.id == data.id).unwrap();
            }
        }
        assert_eq!(
            format!("{}", grid),
            "0: [100]\n\
             1: [400, 400, 400, 400]\n\
             2: [101, 102, 103, 400, 500]\n\
             8: [104, 105, 106]\n\
             13: [200, 300]\n\
             14: [300]\n\
             16: [300]\n\
             17: [300]\n\
             22: [300]\n\
             23: [300]\n\
             25: [300]\n\
             26: [300]\n\
             ids = [100, 101, 102, 103, 104, 105, 106, 200, 300, 400, 500]\n\
             nitem = 11\n\
             ncontainer = 12\n\
             ndiv = [3, 3, 3]\n"
        );
        let mut indices: Vec<_> = grid.containers.into_keys().collect();
        indices.sort();
        assert_eq!(indices, &[0, 1, 2, 8, 13, 14, 16, 17, 22, 23, 25, 26]);
        Ok(())
    }

    #[test]
    fn find_fails_on_wrong_input() {
        let mut g2d = get_test_grid_2d();
        g2d.initialized = false;
        let res = g2d.find(&[0.0, 0.0, 0.0]);
        assert_eq!(res, Err("initialize must be called first"));
        g2d.initialized = true;
        let res = g2d.find(&[0.0, 0.0, 0.0]);
        assert_eq!(res, Err("x.len() must equal ndim"));
        let res = g2d.find(&[0.80001, 0.0]);
        assert_eq!(res, Err("point is outside the grid"));

        let g3d = get_test_grid_3d();
        let res = g3d.find(&[0.0, 0.0]);
        assert_eq!(res, Err("x.len() must equal ndim"));
        let res = g3d.find(&[1.00001, 0.0, 0.0]);
        assert_eq!(res, Err("point is outside the grid"));
    }

    #[test]
    fn find_works() -> Result<(), StrError> {
        let mut g2d = get_test_grid_2d();
        for data in get_test_data_2d() {
            g2d.insert(data.id, data.x)?;
            let id = g2d.find(data.x)?;
            assert_eq!(id, Some(data.id));
        }
        let id = g2d.find(&[0.5, 0.5])?;
        assert_eq!(id, None);

        let mut g3d = get_test_grid_3d();
        for data in get_test_data_3d() {
            g3d.insert(data.id, data.x)?;
            let id = g3d.find(data.x)?;
            assert_eq!(id, Some(data.id));
        }
        let id = g3d.find(&[0.5, 0.5, 0.5])?;
        assert_eq!(id, None);
        Ok(())
    }

    #[test]
    fn containers_near_line_2d_works() -> Result<(), StrError> {
        let mut g2d = get_test_grid_2d();
        for data in get_test_data_2d() {
            g2d.insert(data.id, data.x)?;
        }
        // vertical right-most
        let mut indices = g2d.containers_near_line(&[0.6, 0.0], &[0.6, 1.8])?;
        indices.sort();
        assert_eq!(indices, &[3, 4, 13, 18, 19, 23, 24]);

        // vertical middle
        let mut indices = g2d.containers_near_line(&[0.1 + g2d.radius, 0.0], &[0.1 + g2d.radius, 1.8])?;
        indices.sort();
        assert_eq!(indices, &[1, 3, 6, 7, 11, 12, 13, 17, 18, 21, 23]);

        // horizontal top-most
        let mut indices = g2d.containers_near_line(&[-0.2, 1.8], &[0.8, 1.8])?;
        indices.sort();
        assert_eq!(indices, &[20, 21, 23, 24]);

        // sloped
        let mut indices = g2d.containers_near_line(&[0.2, -0.2], &[0.8, 0.1])?;
        indices.sort();
        assert_eq!(indices, &[1, 3, 4]);
        Ok(())
    }

    #[test]
    fn containers_near_line_3d_works() -> Result<(), StrError> {
        let mut g3d = get_test_grid_3d();
        for data in get_test_data_3d() {
            g3d.insert(data.id, data.x)?;
        }
        let mut indices = g3d.containers_near_line(&[-1.0, -1.0, -1.0], &[1.0, -1.0, -1.0])?;
        indices.sort();
        assert_eq!(indices, &[0, 1, 2]);

        let mut indices = g3d.containers_near_line(&[-1.0, -1.0, -1.0], &[1.0, 1.0, 1.0])?;
        indices.sort();
        assert_eq!(indices, &[0, 1, 13, 14, 16, 17, 22, 23, 25, 26]);
        Ok(())
    }

    #[test]
    fn find_on_line_fails_on_wrong_input() {
        let mut g2d = get_test_grid_2d();
        g2d.initialized = false;
        let res = g2d.find_on_line(&[0.0, 0.0], &[1.0, 1.0]);
        assert_eq!(res, Err("initialize must be called first"));
        g2d.initialized = true;
        let res = g2d.find_on_line(&[0.0], &[1.0, 1.0]);
        assert_eq!(res, Err("a.len() must equal ndim"));
        let res = g2d.find_on_line(&[0.0, 0.0], &[1.0]);
        assert_eq!(res, Err("b.len() must equal ndim"));
    }

    #[test]
    fn find_on_line_2d_works() -> Result<(), StrError> {
        let mut g2d = get_test_grid_2d();
        for data in get_test_data_2d() {
            g2d.insert(data.id, data.x)?;
        }
        // vertical right-most
        let map = g2d.find_on_line(&[0.6, 0.0], &[0.6, 1.8])?;
        let mut ids: Vec<_> = map.iter().collect();
        ids.sort();
        assert_eq!(ids, [&400, &600]);

        // vertical middle
        let map = g2d.find_on_line(&[0.1 + g2d.radius, 0.0], &[0.1 + g2d.radius, 1.8])?;
        let ids: Vec<_> = map.iter().collect();
        assert_eq!(ids.len(), 0);

        // horizontal top-most
        let map = g2d.find_on_line(&[-0.2, 1.8], &[0.8, 1.8])?;
        let mut ids: Vec<_> = map.iter().collect();
        ids.sort();
        assert_eq!(ids, [&103, &500]);

        // sloped
        let map = g2d.find_on_line(&[0.2, -0.2], &[0.8, 0.1])?;
        let mut ids: Vec<_> = map.iter().collect();
        ids.sort();
        assert_eq!(ids, [&600]);
        Ok(())
    }

    #[test]
    fn find_on_line_3d_works() -> Result<(), StrError> {
        let mut g3d = get_test_grid_3d();
        for data in get_test_data_3d() {
            g3d.insert(data.id, data.x)?;
        }
        let map = g3d.find_on_line(&[-1.0, -1.0, -1.0], &[1.0, -1.0, -1.0])?;
        let mut ids: Vec<_> = map.iter().collect();
        ids.sort();
        assert_eq!(ids, [&100, &101]);

        let map = g3d.find_on_line(&[-1.0, -1.0, -1.0], &[1.0, 1.0, 1.0])?;
        let mut ids: Vec<_> = map.iter().collect();
        ids.sort();
        assert_eq!(ids, [&100, &200, &300]);
        Ok(())
    }

    #[test]
    fn containers_near_circle_works() -> Result<(), StrError> {
        let mut g2d = get_test_grid_2d();
        for data in get_test_data_2d() {
            g2d.insert(data.id, data.x)?;
        }
        let mut indices = g2d.containers_near_circle(&[-0.2, 1.8], 0.3)?;
        indices.sort();
        assert_eq!(indices, &[20, 21]);
        Ok(())
    }

    #[test]
    fn find_on_circle_fails_on_wrong_input() {
        let mut g2d = get_test_grid_2d();
        g2d.initialized = false;
        let res = g2d.find_on_circle(&[0.0, 0.0], 1.0);
        assert_eq!(res, Err("initialize must be called first"));
        g2d.initialized = true;
        let res = g2d.find_on_circle(&[-0.2], 0.3);
        assert_eq!(res, Err("center.len() must equal ndim"));

        let g3d = get_test_grid_3d();
        let res = g3d.find_on_circle(&[0.0, 0.0, 0.0], 1.0);
        assert_eq!(res, Err("this works in 2D only"));
    }

    #[test]
    fn find_on_circle_works() -> Result<(), StrError> {
        let mut g2d = get_test_grid_2d();
        for data in get_test_data_2d() {
            g2d.insert(data.id, data.x)?;
        }
        let map = g2d.find_on_circle(&[-0.2, 1.8], 0.3)?;
        let mut ids: Vec<_> = map.iter().collect();
        ids.sort();
        assert_eq!(ids, [&101, &102, &103]);
        Ok(())
    }

    #[test]
    fn containers_near_cylinder_works() -> Result<(), StrError> {
        let mut g3d = get_test_grid_3d();
        for data in get_test_data_3d() {
            g3d.insert(data.id, data.x)?;
        }
        let mut indices = g3d.containers_near_cylinder(&[1.0, -1.0, -1.0], &[1.0, 1.0, -1.0], 0.4)?;
        indices.sort();
        assert_eq!(indices, &[2, 8]);
        Ok(())
    }

    #[test]
    fn find_on_cylinder_fails_on_wrong_input() {
        let mut g3d = get_test_grid_3d();
        g3d.initialized = false;
        let res = g3d.find_on_cylinder(&[0.0, 0.0, 0.0], &[1.0, 0.0, 0.0], 1.0);
        assert_eq!(res, Err("initialize must be called first"));
        g3d.initialized = true;
        let res = g3d.find_on_cylinder(&[0.0, 0.0], &[1.0, 0.0, 0.0], 1.0);
        assert_eq!(res, Err("a.len() must equal ndim"));
        let res = g3d.find_on_cylinder(&[0.0, 0.0, 0.0], &[1.0, 0.0], 1.0);
        assert_eq!(res, Err("b.len() must equal ndim"));

        let g2d = get_test_grid_2d();
        let res = g2d.find_on_cylinder(&[0.0, 0.0, 0.0], &[1.0, 0.0, 0.0], 1.0);
        assert_eq!(res, Err("this works in 3D only"));
    }

    #[test]
    fn find_on_cylinder_works() -> Result<(), StrError> {
        let mut g3d = get_test_grid_3d();
        for data in get_test_data_3d() {
            g3d.insert(data.id, data.x)?;
        }
        let map = g3d.find_on_cylinder(&[1.0, -1.0, -1.0], &[1.0, 1.0, -1.0], 0.4)?;
        let mut ids: Vec<_> = map.iter().collect();
        ids.sort();
        assert_eq!(ids, [&101, &102, &103, &104, &105, &106]);
        Ok(())
    }

    #[test]
    fn containers_near_plane_works() -> Result<(), StrError> {
        let mut g3d = get_test_grid_3d();
        for data in get_test_data_3d() {
            g3d.insert(data.id, data.x)?;
        }
        let mut indices = g3d.containers_near_plane(0, -1.0);
        indices.sort();
        assert_eq!(indices, &[0]);

        let mut indices = g3d.containers_near_plane(1, -1.0);
        indices.sort();
        assert_eq!(indices, &[0, 1, 2]);

        let mut indices = g3d.containers_near_plane(2, -1.0);
        indices.sort();
        assert_eq!(indices, &[0, 1, 2, 8]);

        let mut indices = g3d.containers_near_plane(2, 0.0);
        indices.sort();
        assert_eq!(indices, &[13, 14, 16, 17]);
        Ok(())
    }

    #[test]
    fn find_on_plane_fails_on_wrong_input() {
        let g2d = get_test_grid_2d();
        let mut g3d = get_test_grid_3d();
        g3d.initialized = false;
        let res = g3d.find_on_plane_xy(-1.0);
        assert_eq!(res, Err("initialize must be called first"));
        g3d.initialized = true;
        let res = g2d.find_on_plane_xy(-1.0);
        assert_eq!(res, Err("this works in 3D only"));

        g3d.initialized = false;
        let res = g3d.find_on_plane_yz(-1.0);
        assert_eq!(res, Err("initialize must be called first"));
        g3d.initialized = true;
        let res = g2d.find_on_plane_yz(-1.0);
        assert_eq!(res, Err("this works in 3D only"));

        g3d.initialized = false;
        let res = g3d.find_on_plane_xz(-1.0);
        assert_eq!(res, Err("initialize must be called first"));
        g3d.initialized = true;
        let res = g2d.find_on_plane_xz(-1.0);
        assert_eq!(res, Err("this works in 3D only"));
    }

    #[test]
    fn find_on_plane_works() -> Result<(), StrError> {
        let mut g3d = get_test_grid_3d();
        for data in get_test_data_3d() {
            g3d.insert(data.id, data.x)?;
        }
        let map = g3d.find_on_plane_xy(-1.0)?;
        let mut ids: Vec<_> = map.iter().collect();
        ids.sort();
        assert_eq!(ids, [&100, &101, &104]);

        let map = g3d.find_on_plane_yz(-1.0)?;
        let mut ids: Vec<_> = map.iter().collect();
        ids.sort();
        assert_eq!(ids, [&100]);

        let map = g3d.find_on_plane_xz(-1.0)?;
        let mut ids: Vec<_> = map.iter().collect();
        ids.sort();
        assert_eq!(ids, [&100, &101, &102, &103]);
        Ok(())
    }

    #[test]
    fn maybe_insert_works() -> Result<(), StrError> {
        let mut grid = GridSearch::new(2)?;
        grid.initialize(&[3, 3], &[0.0, 0.0], &[1.0, 1.0])?;
        let id = grid.maybe_insert(111, &[0.5, 0.5])?;
        assert_eq!(id, 111);
        let id = grid.maybe_insert(222, &[0.75, 0.75])?;
        assert_eq!(id, 222);
        let id = grid.maybe_insert(333, &[0.5, 0.5])?;
        assert_eq!(id, 111); // existent
        Ok(())
    }

    #[test]
    fn plot_works() -> Result<(), StrError> {
        let mut g2d = get_test_grid_2d();
        for data in get_test_data_2d() {
            g2d.insert(data.id, data.x)?;
        }
        let mut plot = g2d.plot()?;
        let mut shapes = Shapes::new();
        shapes.set_face_color("None").set_edge_color("magenta");
        shapes.draw_circle(0.5, 0.8, g2d.radius);
        shapes.draw_circle(0.1, 0.4, g2d.radius);
        shapes.draw_circle(0.5, 0.0, g2d.radius);
        shapes.draw_circle(0.1, 1.6, g2d.radius);
        shapes.set_edge_color("#fab32f").set_line_width(1.5);
        shapes.draw_polyline(&vec![vec![0.6, -0.2], vec![0.6, 1.8]], false);
        shapes.draw_polyline(&vec![vec![-0.2, 1.8], vec![0.8, 1.8]], false);
        shapes.draw_polyline(&vec![vec![0.2, -0.2], vec![0.8, 0.1]], false);
        shapes.draw_polyline(&vec![vec![0.1 + g2d.radius, -0.2], vec![0.1 + g2d.radius, 1.8]], false);
        shapes.set_edge_color("green").set_line_width(0.5);
        shapes.draw_circle(-0.2, 1.8, 0.3);
        plot.add(&shapes);
        plot.set_equal_axes(true)
            .set_range(-0.4, 1.0, -0.4, 2.0)
            .set_ticks_x(0.1, 0.0, "")
            .set_num_ticks_y(12)
            .grid_and_labels("x", "y")
            .set_figure_size_points(400.0, 800.0);
        plot.save("/tmp/gemlab/search_grid_plot_2d_works.svg")?;

        let mut g3d = get_test_grid_3d();
        for data in get_test_data_3d() {
            g3d.insert(data.id, data.x)?;
        }
        let mut plot = g3d.plot()?;
        let mut curve = Curve::new();
        curve.set_line_color("#fab32f").set_line_width(1.5);
        curve.draw_3d(&[-1.0, 1.0], &[-1.0, -1.0], &[-1.0, -1.0]);
        curve.draw_3d(&[-1.0, 1.0], &[-1.0, 1.0], &[-1.0, 1.0]);
        let mut surface = Surface::new();
        surface.set_with_surface(false).set_with_wireframe(true);
        surface.draw_cylinder(&[1.0, -1.0, -1.0], &[1.0, 1.0, -1.0], 0.4, 1, 20)?;
        plot.add(&curve);
        plot.add(&surface);
        plot.set_figure_size_points(600.0, 600.0);
        plot.save("/tmp/gemlab/search_grid_plot_3d_works.svg")?;
        // plot.save_and_show("/tmp/gemlab/search_grid_plot_3d_works.svg")?;
        Ok(())
    }

    #[test]
    fn initialize_flag_works() -> Result<(), StrError> {
        let mut grid = GridSearch::new(2)?;
        let res = grid.insert(0, &[0.0, 0.0]);
        assert_eq!(res.err(), Some("initialize must be called first"));

        grid.initialize(&[1, 1], &[0.0, 0.0], &[1.0, 1.0])?;
        let res = grid.insert(0, &[0.0, 0.0]);
        assert_eq!(res.err(), None);
        assert_eq!(grid.initialized, true);

        grid.set_tolerances(&[0.1, 0.1])?;
        assert_eq!(grid.initialized, false);
        let res = grid.insert(0, &[0.0, 0.0]);
        assert_eq!(res.err(), Some("initialize must be called first"));

        let res = grid.find(&[0.0, 0.0]);
        assert_eq!(res.err(), Some("initialize must be called first"));

        let res = grid.plot();
        assert_eq!(res.err(), Some("initialize must be called first"));

        Ok(())
    }

    #[test]
    fn clone_and_serialize_work() -> Result<(), StrError> {
        let mut grid = GridSearch::new(2)?;
        grid.initialize(&[3, 3], &[0.0, 0.0], &[1.0, 1.0])?;
        grid.insert(0, &[0.5, 0.5])?;
        // clone
        let mut cloned = grid.clone();
        cloned.insert(1, &[0.8, 0.8])?;
        assert_eq!(
            format!("{}", grid),
            "4: [0]\n\
             ids = [0]\n\
             nitem = 1\n\
             ncontainer = 1\n\
             ndiv = [3, 3]\n"
        );
        assert_eq!(
            format!("{}", cloned),
            "4: [0]\n\
             8: [1]\n\
             ids = [0, 1]\n\
             nitem = 2\n\
             ncontainer = 2\n\
             ndiv = [3, 3]\n"
        );
        // serialize
        let mut serialized = Vec::new();
        let mut serializer = rmp_serde::Serializer::new(&mut serialized);
        grid.serialize(&mut serializer).map_err(|_| "grid serialize failed")?;
        assert!(serialized.len() > 0);
        // deserialize
        let mut deserializer = rmp_serde::Deserializer::new(&serialized[..]);
        let grid_read: GridSearch =
            Deserialize::deserialize(&mut deserializer).map_err(|_| "cannot deserialize grid data")?;
        assert_eq!(
            format!("{}", grid_read),
            "4: [0]\n\
             ids = [0]\n\
             nitem = 1\n\
             ncontainer = 1\n\
             ndiv = [3, 3]\n"
        );
        Ok(())
    }
}
