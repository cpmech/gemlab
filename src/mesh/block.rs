use super::{Cell, Mesh, Point};
use crate::mesh::PointId;
use crate::shapes::{GeoClass, GeoKind, Scratchpad, HEX_EDGE_TO_FACE};
use crate::util::{AsArray2D, GridSearch};
use crate::StrError;
use plotpy::{Canvas, Plot};
use russell_lab::math::PI;
use russell_lab::Vector;
use std::collections::HashMap;

/// Arguments to transform the Block generated mesh into a ring
///
/// ```text
///   |            /
///   |           / αmax
///   ***=---__  /
///   |         / _
///   |        / ? *._          ,
///   |       / ????? *.     ,-'
///   ***=-_ / ???????? *.,-' αmin
///   |     / - ?????? ,-'*
///   |    /    *.? ,-'    *
///   |   /      ,*'        *
///   |  /    ,-'  *         *
///   | /  ,-'      *         *
///   |/.-'         #         #
///   o ----------- # ------- # --> r
///               rmin       rmax
/// ```
#[derive(Clone, Debug)]
pub struct ArgsRing {
    pub amin: f64,
    pub amax: f64,
    pub rmin: f64,
    pub rmax: f64,
    pub zmin: f64,
    pub zmax: f64,
}

/// Defines constraints for a side of a block (2D)
#[derive(Clone, Debug)]
pub enum Constraint2D {
    /// Circumference specified by (xc,yc,radius)
    Circle(f64, f64, f64),
}

/// Defines constraints for a side of a block (3D)
#[derive(Clone, Debug)]
pub enum Constraint3D {
    /// Surface of a cylinder parallel to x specified by (yc,zc,radius)
    CylinderX(f64, f64, f64),

    /// Surface of a cylinder parallel to y specified by (xc,zc,radius)
    CylinderY(f64, f64, f64),

    /// Surface of a cylinder parallel to z specified by (xc,yc,radius)
    CylinderZ(f64, f64, f64),
}

/// Defines a zero distance to consider point coincident with centre or circle/cylinder
const TOL_DISTANCE: f64 = 1e-8;

/// Defines a polygon on polyhedron that can be split into smaller shapes
///
/// # Block shape
///
/// ## 2D: Local indices of points and edges
///
/// ```text
///             Points (4 or 8)                        Edges (4)
///   y
///   |        3      6      2                             2
///   +--x      @-----@-----@                        +-----------+
///             |           |                        |           |
///             |           |                        |           |
///           7 @           @ 5                     3|           |1
///             |           |                        |           |
///             |           |                        |           |
///             @-----@-----@                        +-----------+
///            0      4      1                             0
/// ```
///
/// ## 3D: Local indices of points, edges, and faces
///
/// ```text
///                   Points (8 or 20)                   Faces (6)
///    z
///    |           4        15        7
///   ,+--y         @-------@--------@                 +----------------+
/// x'            ,'|              ,'|               ,'|              ,'|
///          12 @'  |         14 ,'  |             ,'  |  ___       ,'  |
///           ,'    |16        ,@    |19         ,'    |,'5,'  [0],'    |
///     5   ,'      @      6 ,'      @         ,'      |~~~     ,'      |
///       @'=======@=======@'        |       +'===============+'  ,'|   |
///       |      13 |      |         |       |   ,'|   |      |   |3|   |
///       |         |      |  11     |       |   |2|   |      |   |,'   |
///    17 |       0 @- - - | @- - - -@       |   |,'   +- - - | +- - - -+
///       @       ,'       @       ,' 3      |       ,'       |       ,'
///       |   8 @'      18 |     ,'          |     ,' [1]  ___|     ,'
///       |   ,'           |   ,@ 10         |   ,'      ,'4,'|   ,'
///       | ,'             | ,'              | ,'        ~~~  | ,'
///       @-------@--------@'                +----------------+'
///     1         9         2
///
///      
///                  Edges (12)
///      
///                         7
///                 +----------------+
///               ,'|              ,'|
///           4 ,'  |8          6,'  |
///           ,'    |          ,'    |   
///         ,'      |   5    ,'      |11
///       +'===============+'        |
///       |         |      |         |
///       |         |      |  3      |
///       |         .- - - | -  - - -+
///      9|       ,'       |       ,'
///       |    0,'         |10   ,'
///       |   ,'           |   ,' 2
///       | ,'             | ,'
///       +----------------+'
///               1
/// ```
///
/// # Examples
///
/// ```
/// use gemlab::mesh::Block;
/// use gemlab::shapes::GeoKind;
/// use gemlab::StrError;
///
/// fn main() -> Result<(), StrError> {
///     let mut block = Block::new(&[
///         [0.0, 0.0],
///         [2.0, 0.0],
///         [2.0, 2.0],
///         [0.0, 2.0],
///     ])?;
///     let mesh = block.subdivide(GeoKind::Qua4)?;
///     // 7---------6---------8
///     // |         |         |
///     // |   [2]   |   [3]   |
///     // |         |         |
///     // 3---------2---------5
///     // |         |         |
///     // |   [0]   |   [1]   |
///     // |         |         |
///     // 0---------1---------4
///     assert_eq!(mesh.points.len(), 9);
///     assert_eq!(mesh.cells.len(), 4);
///     Ok(())
/// }
/// ```
pub struct Block {
    /// Attribute of all elements in this block
    attribute: usize,

    /// Space dimension
    ndim: usize,

    /// number of divisions along each dim (ndim)
    ndiv: Vec<usize>,

    /// Delta ksi along each dim (ndim, {ndiv0,ndiv1,ndiv2})
    delta_ksi: Vec<Vec<f64>>,

    /// Constraints on edges (nedge) (2D only)
    edge_constraints: HashMap<usize, Constraint2D>,

    /// Constraints on faces (nface) (3D only)
    face_constraints: HashMap<usize, Constraint3D>,

    /// Has constraints?
    has_constraints: bool,

    /// Scratchpad for the block
    pad: Scratchpad,

    /// Arguments to transform the generated mesh into a ring
    args_ring: ArgsRing,

    /// Transform the generated mesh into a ring
    transform_into_ring: bool,
}

impl Block {
    /// Length of shape along each direction in reference coords space
    const NAT_LENGTH: f64 = 2.0;

    /// Default number of divisions
    const NDIV: usize = 2;

    /// Default attribute
    const ATTRIBUTE: usize = 1;

    /// Allocate a new instance
    ///
    /// # Input (2D)
    ///
    /// * `coords` -- (nrow = 4 or 8, ncol = ndim = 2) matrix with all coordinates
    ///
    /// # Input (3D)
    ///
    /// * `coords` -- (nrow = 8 or 20, ncol = ndim = 3) matrix with all coordinates
    pub fn new<'a, T>(coords: &'a T) -> Result<Self, StrError>
    where
        T: AsArray2D<'a, f64>,
    {
        // check
        let (nrow, ndim) = coords.size();
        if ndim < 2 || ndim > 3 {
            return Err("ncol = ndim must be 2 or 3");
        }
        if ndim == 2 {
            if nrow != 4 && nrow != 8 {
                return Err("in 2D, nrow must be 4 or 8");
            }
        } else {
            if nrow != 8 && nrow != 20 {
                return Err("in 3D, nrow must be 8 or 20");
            }
        }

        // internally, we save the block as a Qua8 or a Hex20
        let kind = if ndim == 2 { GeoKind::Qua8 } else { GeoKind::Hex20 };
        let nnode = kind.nnode();

        // coordinates of the block
        let mut pad = Scratchpad::new(ndim, kind)?;
        if ndim == 2 {
            if nrow == 8 {
                // all vertices given, ok
                for m in 0..nnode {
                    for j in 0..ndim {
                        pad.set_xx(m, j, coords.at(m, j));
                    }
                }
            } else {
                // copy "corner" vertices
                for m in 0..nrow {
                    for j in 0..ndim {
                        pad.set_xx(m, j, coords.at(m, j));
                    }
                }
                // generate mid vertices
                for j in 0..ndim {
                    pad.set_xx(4, j, (coords.at(0, j) + coords.at(1, j)) / 2.0);
                    pad.set_xx(5, j, (coords.at(1, j) + coords.at(2, j)) / 2.0);
                    pad.set_xx(6, j, (coords.at(2, j) + coords.at(3, j)) / 2.0);
                    pad.set_xx(7, j, (coords.at(3, j) + coords.at(0, j)) / 2.0);
                }
            }
        } else {
            if nrow == 20 {
                // all vertices given, ok
                for m in 0..nnode {
                    for j in 0..ndim {
                        pad.set_xx(m, j, coords.at(m, j));
                    }
                }
            } else {
                // copy "corner" vertices
                for m in 0..nrow {
                    for j in 0..ndim {
                        pad.set_xx(m, j, coords.at(m, j));
                    }
                }
                // generate mid vertices
                for j in 0..ndim {
                    pad.set_xx(8, j, (coords.at(0, j) + coords.at(1, j)) / 2.0);
                    pad.set_xx(9, j, (coords.at(1, j) + coords.at(2, j)) / 2.0);
                    pad.set_xx(10, j, (coords.at(2, j) + coords.at(3, j)) / 2.0);
                    pad.set_xx(11, j, (coords.at(3, j) + coords.at(0, j)) / 2.0);

                    pad.set_xx(12, j, (coords.at(4, j) + coords.at(5, j)) / 2.0);
                    pad.set_xx(13, j, (coords.at(5, j) + coords.at(6, j)) / 2.0);
                    pad.set_xx(14, j, (coords.at(6, j) + coords.at(7, j)) / 2.0);
                    pad.set_xx(15, j, (coords.at(7, j) + coords.at(4, j)) / 2.0);

                    pad.set_xx(16, j, (coords.at(0, j) + coords.at(4, j)) / 2.0);
                    pad.set_xx(17, j, (coords.at(1, j) + coords.at(5, j)) / 2.0);
                    pad.set_xx(18, j, (coords.at(2, j) + coords.at(6, j)) / 2.0);
                    pad.set_xx(19, j, (coords.at(3, j) + coords.at(7, j)) / 2.0);
                }
            }
        };

        // block
        Ok(Block {
            attribute: Block::ATTRIBUTE,
            ndim,
            ndiv: vec![Block::NDIV; ndim],
            delta_ksi: vec![vec![1.0; Block::NDIV]; ndim],
            edge_constraints: HashMap::new(),
            face_constraints: HashMap::new(),
            has_constraints: false,
            pad,
            args_ring: ArgsRing {
                amin: 0.0,
                amax: PI / 2.0,
                rmin: 1.0,
                rmax: 2.0,
                zmin: 0.0,
                zmax: 1.0,
            },
            transform_into_ring: false,
        })
    }

    /// Allocates a new square block
    ///
    /// # Input
    ///
    /// * `l` -- the edge length
    #[rustfmt::skip]
    pub fn new_square(l: f64) -> Self {
        Block::new(&[
            [0.0, 0.0],
            [  l, 0.0],
            [  l,   l],
            [0.0,   l],
        ]).unwrap()
    }

    /// Allocates a new cubic block
    ///
    /// # Input
    ///
    /// * `l` -- the edge length
    #[rustfmt::skip]
    pub fn new_cube(l: f64) -> Self {
        Block::new(&[
            [0.0, 0.0, 0.0],
            [  l, 0.0, 0.0],
            [  l,   l, 0.0],
            [0.0,   l, 0.0],
            [0.0, 0.0,   l],
            [  l, 0.0,   l],
            [  l,   l,   l],
            [0.0,   l,   l],
        ]).unwrap()
    }

    /// Draws this block
    pub fn draw(&self, plot: &mut Plot, with_ids: bool, set_range: bool) -> Result<(), StrError> {
        if self.ndim == 2 && self.has_constraints {
            for ct in self.edge_constraints.values() {
                match ct {
                    Constraint2D::Circle(xc, yc, r) => {
                        let mut circle = Canvas::new();
                        circle
                            .set_face_color("None")
                            .set_edge_color("#bfbfbf")
                            .set_line_width(7.0)
                            .draw_circle(xc, yc, r);
                        plot.add(&circle);
                    }
                }
            }
        }
        self.pad.draw_shape(plot, "#046002", with_ids, set_range)
    }

    /// Sets group
    pub fn set_attribute(&mut self, attribute: usize) -> &mut Self {
        self.attribute = attribute;
        self
    }

    /// Sets the number of equal divisions
    ///
    /// For each direction:
    ///
    /// ```text
    /// Δξᵐ = wᵐ ⋅ L / Σ_m wᵐ
    /// ```
    ///
    /// where `L=2` is the edge-length (in reference coordinates) and `wᵐ` are
    /// the weights for each division `m`.
    pub fn set_ndiv(&mut self, ndiv: &[usize]) -> Result<&mut Self, StrError> {
        if ndiv.len() != self.ndim {
            return Err("ndiv.len() must be equal to ndim");
        }
        for i in 0..self.ndim {
            if ndiv[i] < 1 {
                return Err("ndiv must be ≥ 1");
            }
            self.ndiv[i] = ndiv[i];
            let w = 1.0;
            let sum_w = ndiv[i] as f64;
            self.delta_ksi[i] = vec![w * Block::NAT_LENGTH / sum_w; ndiv[i]];
        }
        Ok(self)
    }

    /// Sets (or removes) a constraint to an edge of this block
    ///
    /// # Input
    ///
    /// * `e` -- index of edge; must be < 4 (nedge)
    /// * `constraint` -- the constraint
    ///
    /// # Warning
    ///
    /// **Only one constraint is applied to a point.**
    ///
    /// Thus, if a point is under multiple constraints, one constraint, selected in an unknown order,
    /// will be applied to this point and the other constraints will be ignored for this particular point.
    ///
    /// Therefore, make sure that multiple constraints are compatible one with another. In other words,
    /// the effect of two constraints on a point (e.g. at a corner) must be equal.
    pub fn set_edge_constraint(&mut self, e: usize, value: Option<Constraint2D>) -> Result<&mut Self, StrError> {
        if self.ndim != 2 {
            return Err("set_edge_constraint requires ndim = 2");
        }
        if e > 3 {
            return Err("edge index must be < 4 (nedge)");
        }
        match value {
            Some(v) => {
                self.edge_constraints.insert(e, v);
            }
            None => {
                self.edge_constraints.remove(&e);
            }
        }
        self.has_constraints = self.edge_constraints.len() > 0;
        Ok(self)
    }

    /// Sets (or removes) a constraint to a face of this block
    ///
    /// # Input
    ///
    /// * `f` -- index of face; must be < 6 (nface)
    /// * `constraint` -- the constraint
    ///
    /// # Warning
    ///
    /// **Only one constraint is applied to a point.**
    ///
    /// Thus, if a point is under multiple constraints, one constraint, selected in an unknown order,
    /// will be applied to this point and the other constraints will be ignored for this particular point.
    ///
    /// Therefore, make sure that multiple constraints are compatible one with another. In other words,
    /// the effect of two constraints on a point (e.g. at a corner) must be equal.
    pub fn set_face_constraint(&mut self, f: usize, value: Option<Constraint3D>) -> Result<&mut Self, StrError> {
        if self.ndim != 3 {
            return Err("set_face_constraint requires ndim = 3");
        }
        if f > 5 {
            return Err("face index must be < 6 (nface)");
        }
        match value {
            Some(v) => {
                self.face_constraints.insert(f, v);
            }
            None => {
                self.face_constraints.remove(&f);
            }
        }
        self.has_constraints = self.face_constraints.len() > 0;
        Ok(self)
    }

    /// Sets flag to transform the generated mesh into a ring
    ///
    /// The resulting domain will be the area indicated with "?".
    /// In 3D, an extrusion is applied along the out-of-plane direction.
    ///
    /// ```text
    ///   |            /
    ///   |           / αmax
    ///   ***=---__  /
    ///   |         / _
    ///   |        / ? *._          ,
    ///   |       / ????? *.     ,-'
    ///   ***=-_ / ???????? *.,-' αmin
    ///   |     / - ?????? ,-'*
    ///   |    /    *.? ,-'    *
    ///   |   /      ,*'        *
    ///   |  /    ,-'  *         *
    ///   | /  ,-'      *         *
    ///   |/.-'         #         #
    ///   o ----------- # ------- # --> r
    ///               rmin       rmax
    /// ```
    ///
    /// Intermediary mapping:
    ///
    /// r(ξ₀,ξ₁,ξ₂) = rmin + (ξ₀ - ξ₀min) · Δr / Δξ₀
    /// α(ξ₀,ξ₁,ξ₂) = αmin + (ξ₁ - ξ₁min) · Δα / Δξ₁
    /// z(ξ₀,ξ₁,ξ₂) = ξ₂
    ///
    /// Cylindrical coordinates:
    ///
    /// x₀ := r · cos(α)
    /// x₁ := r · sin(α)
    /// x₂ := z
    ///
    /// # Input
    ///
    /// * `args` -- the parameters shown in the figure above with the **angles in radians**
    ///   If args is None, then transform-into-ring will be disabled
    pub fn set_transform_into_ring(&mut self, args: Option<ArgsRing>) -> Result<&mut Self, StrError> {
        match args {
            Some(a) => {
                if a.amax <= a.amin {
                    return Err("amax must be greater than amin");
                }
                if a.rmax <= a.rmin {
                    return Err("rmax must be greater than rmin");
                }
                if a.zmax <= a.zmin {
                    return Err("zmax must be greater than zmin");
                }
                if a.amin < 0.0 || a.amin > PI / 2.0 {
                    return Err("amin must be given in radians and must be ≥ 0 and ≤ PI/2");
                }
                if a.amax < 0.0 || a.amax > PI / 2.0 {
                    return Err("amax must be given in radians and must be ≥ 0 and ≤ PI/2");
                }
                self.args_ring = a;
                self.transform_into_ring = true;
            }
            None => self.transform_into_ring = false,
        }
        Ok(self)
    }

    /// Subdivide block into vertices and cells (mesh)
    ///
    /// # Input
    ///
    /// * `target` -- If 2D, a quadrilateral (must have [GeoClass::Qua])
    ///               If 3D, a hexahedron (must have [GeoClass::Hex])
    pub fn subdivide(&mut self, target: GeoKind) -> Result<Mesh, StrError> {
        // check
        let ndim = self.ndim;
        if ndim == 2 {
            if target.class() != GeoClass::Qua {
                return Err("in 2D, the GeoClass of target must be Qua");
            }
        } else {
            if target.class() != GeoClass::Hex {
                return Err("in 3D, the GeoClass of target must be Hex");
            }
        }

        // grid to search reference coordinates
        let mut grid_ksi = GridSearch::new(&vec![-1.0; ndim], &vec![1.0; ndim], None, None, None)?;

        // resulting mesh
        let mut mesh = Mesh {
            ndim,
            points: Vec::new(),
            cells: Vec::new(),
        };

        // constants
        let target_nnode = target.nnode();

        // constants used only if there are constraints and
        // the constraint causes some middle edge nodes to move
        let mut pad_lin2 = Scratchpad::new(ndim, GeoKind::Lin2)?;
        let target_nedge = target.nedge();
        let target_edge_kind = target.edge_kind().unwrap();
        let target_edge_nnode = target.edge_nnode();
        let target_n_interior_nodes = target.n_interior_nodes();
        let mut serendipity = if self.has_constraints {
            match target {
                GeoKind::Qua9 => Some(Scratchpad::new(ndim, GeoKind::Qua8)?),
                GeoKind::Qua16 => Some(Scratchpad::new(ndim, GeoKind::Qua12)?),
                GeoKind::Qua17 => Some(Scratchpad::new(ndim, GeoKind::Qua8)?), // only option available
                _ => None,
            }
        } else {
            None
        };

        // maps the id of a constrained point to the side (edge/face) where
        // the constrained has been applied (i.e., the point coordinates were modified)
        let mut constrained_point_to_side: HashMap<PointId, usize> = HashMap::new(); // point_id => side

        // The idea here is to map the TARGET CELL (shown on the right)
        // to the BLOCK's reference (natural) space. (right-to-left mapping)
        // We keep track of the points already added in the natural space of
        // the block to avoid duplicates (using the GridSearch tool).
        // Then we use the isoparametric formula to convert the reference
        // coordinates in the block space to the real space (that's why
        // fn_interp and pad are defined for the self.kind of the block).
        //          BLOCK                      TARGET CELL
        // (-1,+1)-----------(+1,+1)    (-1,+1)-----------(+1,+1)
        //    |                 |   _.~'   |                 |
        //    |                 _.~'       |                 |
        //    |      o------o ~'|   scale  |                 |
        //    |      |      |   | <<  &    |                 |
        //    |      |      |   | translate|                 |
        //    |      o------o .___         |                 |
        //    |                 | '~~~._   |                 |
        // (-1,-1)-----------(+1,-1)    (-1,-1)-----------(+1,-1)

        // center of new cell in reference space (natural coordinates)
        let mut center = vec![0.0; ndim];

        // size ratios between new cell and block in natural coordinates
        let mut scale = vec![0.0; ndim];

        // natural coordinates of new points
        let mut ksi = vec![0.0; ndim];

        // real coordinates of new points
        let mut x = Vector::new(ndim);

        // number of divisions along each direction
        let (nx, ny, nz) = (self.ndiv[0], self.ndiv[1], if ndim == 2 { 1 } else { self.ndiv[2] });

        // for each z-division
        if ndim == 3 {
            center[2] = -1.0 + self.delta_ksi[2][0] / 2.0;
        }
        for k in 0..nz {
            if ndim == 3 {
                scale[2] = self.delta_ksi[2][k] / Block::NAT_LENGTH;
            }

            // for each y-division
            center[1] = -1.0 + self.delta_ksi[1][0] / 2.0;
            for j in 0..ny {
                scale[1] = self.delta_ksi[1][j] / Block::NAT_LENGTH;

                // for each x-division
                center[0] = -1.0 + self.delta_ksi[0][0] / 2.0;
                for i in 0..nx {
                    scale[0] = self.delta_ksi[0][i] / Block::NAT_LENGTH;

                    // new cell id
                    let cell_id = mesh.cells.len();

                    // for each node of the new (target) cell
                    // (there may be more nodes than the block; e.g., internal nodes)
                    let mut points = vec![0; target_nnode];
                    for m in 0..target_nnode {
                        // reference natural coordinates of the new cell nodes
                        let ksi_ref = target.reference_coords(m);

                        // scale and translate the reference coordinates
                        // to the space of the block
                        for a in 0..ndim {
                            ksi[a] = center[a] + scale[a] * ksi_ref[a];
                        }

                        // look up existent point in the reference space of the block
                        // or create new point
                        let point_id = match grid_ksi.find(&ksi)? {
                            Some(point_id) => point_id,
                            None => {
                                // insert point id in grid-search
                                let point_id = mesh.points.len();
                                grid_ksi.insert(point_id, &ksi)?;

                                // reference coordinates are now scaled and translated
                                if self.transform_into_ring {
                                    self.map_coords_cylindrical(&mut x, &ksi);
                                } else {
                                    // compute real coordinates of point using the block's
                                    self.pad.calc_coords(&mut x, &ksi)?;
                                    if self.has_constraints {
                                        if let Some(side) = self.handle_constraints(&mut x, &ksi)? {
                                            constrained_point_to_side.insert(point_id, side);
                                        }
                                    }
                                }

                                // add new point to mesh
                                mesh.points.push(Point {
                                    id: point_id,
                                    marker: 0,
                                    coords: x.as_data().clone(),
                                });
                                point_id
                            }
                        };
                        points[m] = point_id;
                    }

                    // fix middle edge nodes after corner movement due to constraints
                    // (only process edges with more than 2 nodes; i.e., those with middle/interior nodes)
                    if self.has_constraints && target_edge_nnode > 2 {
                        for (point_id, side) in &constrained_point_to_side {
                            for e in 0..target_nedge {
                                // only process an edge that is orthogonal to the constrained side (edge/face)
                                if ndim == 2 {
                                    if e == *side {
                                        continue; // skip edge that is constrained (nothing should be changed on that edge)
                                    }
                                } else {
                                    if HEX_EDGE_TO_FACE[e][0] == *side || HEX_EDGE_TO_FACE[e][1] == *side {
                                        continue; // skip edge that belongs to the constrained side (face) (nothing should be changed on that edge)
                                    }
                                }
                                // check if the edge contains the constrained point
                                let mut edge_contains_point = false;
                                for idx in 0..target_edge_nnode {
                                    if points[target.edge_node_id(e, idx)] == *point_id {
                                        edge_contains_point = true;
                                        break;
                                    }
                                }
                                if !edge_contains_point {
                                    continue; // skip edge that doesn't contain the constrained point
                                }
                                // set Lin2 with the first two (corner) points of target's edge
                                mesh.set_pad(
                                    &mut pad_lin2,
                                    &[points[target.edge_node_id(e, 0)], points[target.edge_node_id(e, 1)]],
                                );
                                // only change the coordinates of the middle edge nodes: idx >= 2
                                for idx in 2..target_edge_nnode {
                                    let p = points[target.edge_node_id(e, idx)];
                                    // use a Lin2 with the natural coords of target to calculate
                                    // the real position of the interior/middle node
                                    let ksi_edge = target_edge_kind.reference_coords(idx);
                                    pad_lin2.calc_coords(&mut x, ksi_edge)?;
                                    for dim in 0..ndim {
                                        mesh.points[p].coords[dim] = x[dim];
                                    }
                                }
                            }
                        }
                        // now get the serendipity version of target and use it to calculate the
                        // interior coordinates of a cell with interior points (e.g. Qua16)
                        // (this is required to maintain the interior nodes "aligned" with the
                        // just moved middle edge nodes)
                        if target_n_interior_nodes > 0 {
                            if let Some(ref mut pad_ser) = serendipity {
                                let nn = pad_ser.interp.dim();
                                let pts = points[0..nn].to_vec();
                                mesh.set_pad(pad_ser, &pts);
                                for idx in 0..target_n_interior_nodes {
                                    let m = target.interior_node(idx);
                                    let ksi_interior = target.reference_coords(m);
                                    pad_ser.calc_coords(&mut x, ksi_interior)?;
                                    for dim in 0..ndim {
                                        mesh.points[points[m]].coords[dim] = x[dim];
                                    }
                                }
                            }
                        }
                    }

                    // new cell
                    let cell = Cell {
                        id: cell_id,
                        attribute: self.attribute,
                        kind: target,
                        points,
                    };
                    mesh.cells.push(cell);

                    // next x-center
                    center[0] += self.delta_ksi[0][i];
                }

                // next y-center
                center[1] += self.delta_ksi[1][j];
            }

            // next z-center
            if ndim == 3 {
                center[2] += self.delta_ksi[2][k];
            }
        }
        Ok(mesh)
    }

    /// Maps reference coordinates to real coordinates on a cylinder
    #[inline]
    fn map_coords_cylindrical(&self, x: &mut Vector, ksi: &[f64]) {
        assert_eq!(x.dim(), ksi.len());
        const KSI_MIN: f64 = -1.0;
        const KSI_DEL: f64 = 2.0;
        let (amin, amax) = (self.args_ring.amin, self.args_ring.amax);
        let (rmin, rmax) = (self.args_ring.rmin, self.args_ring.rmax);
        let r = rmin + (ksi[0] - KSI_MIN) * (rmax - rmin) / KSI_DEL;
        let a = amin + (ksi[1] - KSI_MIN) * (amax - amin) / KSI_DEL;
        x[0] = r * f64::cos(a);
        x[1] = r * f64::sin(a);
        if x.dim() == 3 {
            let (zmin, zmax) = (self.args_ring.zmin, self.args_ring.zmax);
            x[2] = zmin + (ksi[2] - KSI_MIN) * (zmax - zmin) / KSI_DEL;
        }
    }

    /// Handles constraints
    ///
    /// Returns the local index of the edge or face corresponding to the side
    /// where the constraint has been applied.
    /// Returns None if the point is not on a constrained side.
    fn handle_constraints(&self, x: &mut Vector, ksi: &[f64]) -> Result<Option<usize>, StrError> {
        const TOL: f64 = 1e-13;
        if self.ndim == 2 {
            if f64::abs(-1.0 - ksi[1]) < TOL {
                // bottom => e = 0
                if let Some(ct) = self.edge_constraints.get(&0) {
                    if self.apply_constraint_2d(x, ct)? {
                        return Ok(Some(0));
                    }
                }
            }
            if f64::abs(1.0 - ksi[0]) < TOL {
                // right => e = 1
                if let Some(ct) = self.edge_constraints.get(&1) {
                    if self.apply_constraint_2d(x, ct)? {
                        return Ok(Some(1));
                    }
                }
            }
            if f64::abs(1.0 - ksi[1]) < TOL {
                // top => e = 2
                if let Some(ct) = self.edge_constraints.get(&2) {
                    if self.apply_constraint_2d(x, ct)? {
                        return Ok(Some(2));
                    }
                }
            }
            if f64::abs(-1.0 - ksi[0]) < TOL {
                // left => e = 3
                if let Some(ct) = self.edge_constraints.get(&3) {
                    if self.apply_constraint_2d(x, ct)? {
                        return Ok(Some(3));
                    }
                }
            }
        } else {
            if f64::abs(-1.0 - ksi[0]) < TOL {
                // negative x => f = 0
                if let Some(ct) = self.face_constraints.get(&0) {
                    if self.apply_constraint_3d(x, ct)? {
                        return Ok(Some(0));
                    }
                }
            }
            if f64::abs(1.0 - ksi[0]) < TOL {
                // positive x => f = 1
                if let Some(ct) = self.face_constraints.get(&1) {
                    if self.apply_constraint_3d(x, ct)? {
                        return Ok(Some(1));
                    }
                }
            }
            if f64::abs(-1.0 - ksi[1]) < TOL {
                // negative y => f = 2
                if let Some(ct) = self.face_constraints.get(&2) {
                    if self.apply_constraint_3d(x, ct)? {
                        return Ok(Some(2));
                    }
                }
            }
            if f64::abs(1.0 - ksi[1]) < TOL {
                // positive y => f = 3
                if let Some(ct) = self.face_constraints.get(&3) {
                    if self.apply_constraint_3d(x, ct)? {
                        return Ok(Some(3));
                    }
                }
            }
            if f64::abs(-1.0 - ksi[2]) < TOL {
                // negative z => f = 4
                if let Some(ct) = self.face_constraints.get(&4) {
                    if self.apply_constraint_3d(x, ct)? {
                        return Ok(Some(4));
                    }
                }
            }
            if f64::abs(1.0 - ksi[2]) < TOL {
                // positive z => f = 5
                if let Some(ct) = self.face_constraints.get(&5) {
                    if self.apply_constraint_3d(x, ct)? {
                        return Ok(Some(5));
                    }
                }
            }
        }
        Ok(None)
    }

    /// Applies 2D constraint
    ///
    /// Returns true if the constraint has been applied; otherwise, returns false
    fn apply_constraint_2d(&self, x: &mut Vector, ct: &Constraint2D) -> Result<bool, StrError> {
        match ct {
            Constraint2D::Circle(xc, yc, r) => {
                let dx = x[0] - xc;
                let dy = x[1] - yc;
                let d = f64::sqrt(dx * dx + dy * dy);
                if f64::abs(d) <= TOL_DISTANCE {
                    return Err("cannot apply constraint because a point is at the center of the circle");
                }
                let gap = r - d;
                if f64::abs(gap) > 0.0 {
                    let move_x = gap * dx / d;
                    let move_y = gap * dy / d;
                    x[0] += move_x;
                    x[1] += move_y;
                    return Ok(true);
                }
            }
        }
        Ok(false)
    }

    /// Applies 3D constraint
    ///
    /// Returns true if the constraint has been applied; otherwise, returns false
    fn apply_constraint_3d(&self, x: &mut Vector, ct: &Constraint3D) -> Result<bool, StrError> {
        match ct {
            Constraint3D::CylinderX(yc, zc, r) => {
                let dy = x[1] - yc;
                let dz = x[2] - zc;
                let d = f64::sqrt(dy * dy + dz * dz);
                if f64::abs(d) <= TOL_DISTANCE {
                    return Err("cannot apply constraint because a point is at the center of the cylinder-x");
                }
                let gap = r - d;
                if f64::abs(gap) > 0.0 {
                    let move_y = gap * dy / d;
                    let move_z = gap * dz / d;
                    x[1] += move_y;
                    x[2] += move_z;
                    return Ok(true);
                }
            }
            Constraint3D::CylinderY(xc, zc, r) => {
                let dx = x[0] - xc;
                let dz = x[2] - zc;
                let d = f64::sqrt(dx * dx + dz * dz);
                if f64::abs(d) <= TOL_DISTANCE {
                    return Err("cannot apply constraint because a point is at the center of the cylinder-y");
                }
                let gap = r - d;
                if f64::abs(gap) > 0.0 {
                    let move_x = gap * dx / d;
                    let move_z = gap * dz / d;
                    x[0] += move_x;
                    x[2] += move_z;
                    return Ok(true);
                }
            }
            Constraint3D::CylinderZ(xc, yc, r) => {
                let dx = x[0] - xc;
                let dy = x[1] - yc;
                let d = f64::sqrt(dx * dx + dy * dy);
                if f64::abs(d) <= TOL_DISTANCE {
                    return Err("cannot apply constraint because a point is at the center of the cylinder-z");
                }
                let gap = r - d;
                if f64::abs(gap) > 0.0 {
                    let move_x = gap * dx / d;
                    let move_y = gap * dy / d;
                    x[0] += move_x;
                    x[1] += move_y;
                    return Ok(true);
                }
            }
        }
        Ok(false)
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::{ArgsRing, Block, Constraint2D, Constraint3D};
    use crate::geometry::point_point_distance;
    use crate::mesh::{Figure, Mesh, Samples};
    use crate::shapes::GeoKind;
    use plotpy::{Canvas, Plot, Surface};
    use russell_chk::{approx_eq, vec_approx_eq};
    use russell_lab::math::{PI, SQRT_2};

    const SAVE_FIGURE: bool = false;

    fn draw_ring(mesh: &Mesh, args: &ArgsRing, filename: &str) {
        let mut circle_in = Canvas::new();
        let mut circle_mid = Canvas::new();
        let mut circle_out = Canvas::new();
        circle_in
            .set_face_color("None")
            .set_edge_color("#bfbfbfbb")
            .set_line_width(7.0)
            .draw_circle(0.0, 0.0, args.rmin);
        circle_mid
            .set_face_color("None")
            .set_edge_color("#bfbfbfbb")
            .set_line_width(7.0)
            .draw_circle(0.0, 0.0, (args.rmax + args.rmin) / 2.0);
        circle_out
            .set_face_color("None")
            .set_edge_color("#bfbfbfbb")
            .set_line_width(7.0)
            .draw_circle(0.0, 0.0, args.rmax);
        mesh.draw(None, filename, |plot, before| {
            if !before {
                plot.add(&circle_in);
                plot.add(&circle_mid);
                plot.add(&circle_out);
            }
        })
        .unwrap();
    }

    fn draw<F>(mesh: Mesh, block: &Block, set_range: bool, filename: &str, mut pre: F)
    where
        F: FnMut(&mut Plot),
    {
        let mut fig = Figure::new();
        fig.cell_ids = true;
        fig.point_ids = true;
        fig.point_dots = true;
        fig.figure_size = Some((600.0, 600.0));
        fig.canvas_point_ids
            .set_bbox(false)
            .set_align_horizontal("left")
            .set_align_vertical("bottom");
        mesh.draw(Some(fig), filename, |plot, before| {
            if !before {
                pre(plot);
                block.draw(plot, false, set_range).unwrap();
            }
        })
        .unwrap();
    }

    #[test]
    fn derive_works() {
        let args = ArgsRing {
            amin: 1.0,
            amax: 2.0,
            rmin: 3.0,
            rmax: 4.0,
            zmin: 5.0,
            zmax: 6.0,
        };
        let clone = args.clone();
        let correct = "ArgsRing { amin: 1.0, amax: 2.0, rmin: 3.0, rmax: 4.0, zmin: 5.0, zmax: 6.0 }";
        assert_eq!(format!("{:?}", clone), correct);

        let constraint = Constraint2D::Circle(2.0, 3.0, 1.0);
        let clone = constraint.clone();
        let correct = "Circle(2.0, 3.0, 1.0)";
        assert_eq!(format!("{:?}", constraint), correct);
        assert_eq!(format!("{:?}", clone), correct);

        let constraint = Constraint3D::CylinderZ(2.0, 3.0, 1.0);
        let clone = constraint.clone();
        let correct = "CylinderZ(2.0, 3.0, 1.0)";
        assert_eq!(format!("{:?}", constraint), correct);
        assert_eq!(format!("{:?}", clone), correct);
    }

    #[test]
    fn new_fails_on_wrong_input() {
        assert_eq!(Block::new(&[[]]).err(), Some("ncol = ndim must be 2 or 3"));
        assert_eq!(Block::new(&[[0.0, 0.0]]).err(), Some("in 2D, nrow must be 4 or 8"));
        assert_eq!(
            Block::new(&[[0.0, 0.0, 0.0]]).err(),
            Some("in 3D, nrow must be 8 or 20")
        );
    }

    #[test]
    fn new_works() {
        let b2d = Block::new(&[[0.0, 0.0], [2.0, 0.0], [2.0, 2.0], [0.0, 2.0]]).unwrap();
        assert_eq!(b2d.attribute, 1);
        assert_eq!(b2d.ndim, 2);
        assert_eq!(b2d.ndiv, &[2, 2]);
        assert_eq!(format!("{:?}", b2d.delta_ksi), "[[1.0, 1.0], [1.0, 1.0]]");
        assert_eq!(
            format!("{}", b2d.pad.xxt),
            "┌                 ┐\n\
             │ 0 2 2 0 1 2 1 0 │\n\
             │ 0 0 2 2 0 1 2 1 │\n\
             └                 ┘"
        );

        let b2d = Block::new(&[
            [0.0, 0.0],
            [2.0, 0.0],
            [2.0, 2.0],
            [0.0, 2.0],
            [1.0, 0.0],
            [2.0, 1.0],
            [1.0, 2.0],
            [0.0, 1.0],
        ])
        .unwrap();
        assert_eq!(
            format!("{}", b2d.pad.xxt),
            "┌                 ┐\n\
             │ 0 2 2 0 1 2 1 0 │\n\
             │ 0 0 2 2 0 1 2 1 │\n\
             └                 ┘"
        );

        let b3d = Block::new(&[
            [0.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [2.0, 2.0, 0.0],
            [0.0, 2.0, 0.0],
            [0.0, 0.0, 2.0],
            [2.0, 0.0, 2.0],
            [2.0, 2.0, 2.0],
            [0.0, 2.0, 2.0],
        ])
        .unwrap();
        assert_eq!(b3d.attribute, 1);
        assert_eq!(b3d.ndim, 3);
        assert_eq!(
            format!("{}", b3d.pad.xxt),
            "┌                                         ┐\n\
             │ 0 2 2 0 0 2 2 0 1 2 1 0 1 2 1 0 0 2 2 0 │\n\
             │ 0 0 2 2 0 0 2 2 0 1 2 1 0 1 2 1 0 0 2 2 │\n\
             │ 0 0 0 0 2 2 2 2 0 0 0 0 2 2 2 2 1 1 1 1 │\n\
             └                                         ┘"
        );
        assert_eq!(b3d.ndiv, &[2, 2, 2]);
        assert_eq!(format!("{:?}", b3d.delta_ksi), "[[1.0, 1.0], [1.0, 1.0], [1.0, 1.0]]");

        let b3d = Block::new(&[
            [0.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [2.0, 2.0, 0.0],
            [0.0, 2.0, 0.0],
            [0.0, 0.0, 2.0],
            [2.0, 0.0, 2.0],
            [2.0, 2.0, 2.0],
            [0.0, 2.0, 2.0],
            [1.0, 0.0, 0.0],
            [2.0, 1.0, 0.0],
            [1.0, 2.0, 0.0],
            [0.0, 1.0, 0.0],
            [1.0, 0.0, 2.0],
            [2.0, 1.0, 2.0],
            [1.0, 2.0, 2.0],
            [0.0, 1.0, 2.0],
            [0.0, 0.0, 1.0],
            [2.0, 0.0, 1.0],
            [2.0, 2.0, 1.0],
            [0.0, 2.0, 1.0],
        ])
        .unwrap();
        assert_eq!(
            format!("{}", b3d.pad.xxt),
            "┌                                         ┐\n\
             │ 0 2 2 0 0 2 2 0 1 2 1 0 1 2 1 0 0 2 2 0 │\n\
             │ 0 0 2 2 0 0 2 2 0 1 2 1 0 1 2 1 0 0 2 2 │\n\
             │ 0 0 0 0 2 2 2 2 0 0 0 0 2 2 2 2 1 1 1 1 │\n\
             └                                         ┘"
        );
    }

    #[test]
    fn set_attribute_works() {
        let mut block = Block::new_square(1.0);
        block.set_attribute(2);
        assert_eq!(block.attribute, 2);
    }

    #[test]
    fn set_ndiv_works() {
        let mut block = Block::new_square(1.0);
        assert_eq!(block.set_ndiv(&[1]).err(), Some("ndiv.len() must be equal to ndim"));
        assert_eq!(block.set_ndiv(&[0, 1]).err(), Some("ndiv must be ≥ 1"));
        block.set_ndiv(&[2, 4]).unwrap();
        assert_eq!(block.ndiv, &[2, 4]);
        assert_eq!(format!("{:?}", block.delta_ksi), "[[1.0, 1.0], [0.5, 0.5, 0.5, 0.5]]");
    }

    #[test]
    fn set_edge_constraint_works() {
        let mut block = Block::new_square(1.0);
        let ct = Constraint2D::Circle(-1.0, -1.0, 2.0);
        assert_eq!(format!("{:?}", block.edge_constraints), "{}");
        assert_eq!(block.has_constraints, false);
        block.set_edge_constraint(0, Some(ct)).unwrap();
        assert_eq!(format!("{:?}", block.edge_constraints), "{0: Circle(-1.0, -1.0, 2.0)}");
        assert_eq!(block.has_constraints, true);
        block.set_edge_constraint(0, None).unwrap();
        assert_eq!(format!("{:?}", block.edge_constraints), "{}");
        assert_eq!(block.has_constraints, false);
        assert_eq!(
            block.set_edge_constraint(4, None).err(),
            Some("edge index must be < 4 (nedge)")
        );
        let mut block_3d = Block::new_cube(1.0);
        assert_eq!(
            block_3d.set_edge_constraint(0, None).err(),
            Some("set_edge_constraint requires ndim = 2")
        );
    }

    #[test]
    fn set_face_constraint_works() {
        let mut block = Block::new_cube(1.0);
        let ct = Constraint3D::CylinderZ(-1.0, -1.0, 2.0);
        assert_eq!(format!("{:?}", block.face_constraints), "{}");
        assert_eq!(block.has_constraints, false);
        block.set_face_constraint(0, Some(ct)).unwrap();
        assert_eq!(
            format!("{:?}", block.face_constraints),
            "{0: CylinderZ(-1.0, -1.0, 2.0)}"
        );
        assert_eq!(block.has_constraints, true);
        block.set_face_constraint(0, None).unwrap();
        assert_eq!(format!("{:?}", block.face_constraints), "{}");
        assert_eq!(block.has_constraints, false);
        assert_eq!(
            block.set_face_constraint(6, None).err(),
            Some("face index must be < 6 (nface)")
        );
        let mut block_2d = Block::new_square(1.0);
        assert_eq!(
            block_2d.set_face_constraint(0, None).err(),
            Some("set_face_constraint requires ndim = 3")
        );
    }

    #[test]
    fn set_transform_into_ring_works() {
        let mut block = Block::new_square(1.0);
        assert_eq!(block.transform_into_ring, false);
        let mut args = ArgsRing {
            amin: PI,
            amax: PI / 2.0,
            rmin: 3.0,
            rmax: 4.0,
            zmin: 5.0,
            zmax: 6.0,
        };
        assert_eq!(
            block.set_transform_into_ring(Some(args.clone())).err(),
            Some("amax must be greater than amin")
        );
        assert_eq!(block.transform_into_ring, false);
        args.amin = PI / 3.0;
        args.rmin = 5.0;
        assert_eq!(
            block.set_transform_into_ring(Some(args.clone())).err(),
            Some("rmax must be greater than rmin")
        );
        assert_eq!(block.transform_into_ring, false);
        args.rmin = 3.0;
        args.zmin = 7.0;
        assert_eq!(
            block.set_transform_into_ring(Some(args.clone())).err(),
            Some("zmax must be greater than zmin")
        );
        assert_eq!(block.transform_into_ring, false);
        args.zmin = 5.0;
        args.amin = -PI / 2.0;
        assert_eq!(
            block.set_transform_into_ring(Some(args.clone())).err(),
            Some("amin must be given in radians and must be ≥ 0 and ≤ PI/2")
        );
        assert_eq!(block.transform_into_ring, false);
        args.amin = PI / 3.0;
        args.amax = 2.0 * PI;
        assert_eq!(
            block.set_transform_into_ring(Some(args.clone())).err(),
            Some("amax must be given in radians and must be ≥ 0 and ≤ PI/2")
        );
        assert_eq!(block.transform_into_ring, false);
        args.amax = PI / 2.0;
        block.set_transform_into_ring(Some(args)).unwrap();
        assert_eq!(block.transform_into_ring, true);
        block.set_transform_into_ring(None).unwrap();
        assert_eq!(block.transform_into_ring, false);
    }

    #[test]
    fn subdivide_fails_on_wrong_input() {
        let mut b2d = Block::new_square(1.0);
        assert_eq!(
            b2d.subdivide(GeoKind::Tri3).err(),
            Some("in 2D, the GeoClass of target must be Qua")
        );
        let mut b3d = Block::new_cube(1.0);
        assert_eq!(
            b3d.subdivide(GeoKind::Tet4).err(),
            Some("in 3D, the GeoClass of target must be Hex")
        );
    }

    #[test]
    fn subdivide_2d_qua4_works() {
        // 7---------------6---------------8
        // |               |               |
        // |               |               |
        // |      [2]      |      [3]      |
        // |               |               |
        // |               |               |
        // 3---------------2---------------5
        // |               |               |
        // |               |               |
        // |      [0]      |      [1]      |
        // |               |               |
        // |               |               |
        // 0---------------1---------------4
        #[rustfmt::skip]
        let mut block = Block::new(&[
            [0.0, 0.0],
            [2.0, 0.0],
            [2.0, 2.0],
            [0.0, 2.0],
        ]).unwrap();
        let mesh = block.subdivide(GeoKind::Qua4).unwrap();
        let correct = Samples::block_2d_four_qua4();
        assert_eq!(format!("{:?}", mesh), format!("{:?}", correct));
        mesh.check_all().unwrap();
    }

    #[test]
    fn subdivide_works_on_the_same_block_twice() {
        #[rustfmt::skip]
        let mut block = Block::new(&[
            [0.0, 0.0],
            [2.0, 0.0],
            [2.0, 2.0],
            [0.0, 2.0],
        ]).unwrap();
        let mesh = block.subdivide(GeoKind::Qua4).unwrap();
        let correct = Samples::block_2d_four_qua4();
        assert_eq!(format!("{:?}", mesh), format!("{:?}", correct));
        let mesh = block.subdivide(GeoKind::Qua4).unwrap();
        assert_eq!(format!("{:?}", mesh), format!("{:?}", correct));
    }

    #[test]
    fn subdivide_2d_qua8_works() {
        // 14------16------13------20------18
        //  |               |               |
        //  |               |               |
        // 17      [2]     15      [3]     19
        //  |               |               |
        //  |               |               |
        //  3-------6-------2------12-------9
        //  |               |               |
        //  |               |               |
        //  7      [0]      5      [1]     11
        //  |               |               |
        //  |               |               |
        //  0-------4-------1------10-------8
        #[rustfmt::skip]
        let mut block = Block::new(&[
            [0.0, 0.0],
            [2.0, 0.0],
            [2.0, 2.0],
            [0.0, 2.0],
        ]).unwrap();
        let mesh = block.subdivide(GeoKind::Qua8).unwrap();
        let correct = Samples::block_2d_four_qua8();
        assert_eq!(format!("{:?}", mesh), format!("{:?}", correct));
        mesh.check_all().unwrap();
    }

    #[test]
    fn subdivide_2d_qua9_works() {
        // 16------18------15------23------21
        //  |               |               |
        //  |               |               |
        // 19      20      17      24      22
        //  |               |               |
        //  |               |               |
        //  3-------6-------2------13------10
        //  |               |               |
        //  |               |               |
        //  7       8       5      14      12
        //  |               |               |
        //  |               |               |
        //  0-------4-------1------11-------9
        #[rustfmt::skip]
        let mut block = Block::new(&[
            [0.0, 0.0],
            [2.0, 0.0],
            [2.0, 2.0],
            [0.0, 2.0],
        ]).unwrap();
        let mesh = block.subdivide(GeoKind::Qua9).unwrap();
        let correct = Samples::block_2d_four_qua9();
        assert_eq!(format!("{:?}", mesh), format!("{:?}", correct));
        mesh.check_all().unwrap();
    }

    #[test]
    fn subdivide_2d_qua12_works() {
        // 21---26---23----20---32---30----28
        //  |               |               |
        // 24              25              31
        //  |               |               |
        // 27              22              29
        //  |               |               |
        //  3---10-----6----2---19---16----13
        //  |               |               |
        //  7               9              18
        //  |               |               |
        // 11               5              15
        //  |               |               |
        //  0----4-----8----1---14---17----12
        #[rustfmt::skip]
        let mut block = Block::new(&[
            [0.0, 0.0],
            [3.0, 0.0],
            [3.0, 3.0],
            [0.0, 3.0],
        ]).unwrap();
        let mesh = block.subdivide(GeoKind::Qua12).unwrap();
        let correct = Samples::block_2d_four_qua12();
        assert_eq!(mesh.points.len(), correct.points.len());
        assert_eq!(mesh.cells.len(), correct.cells.len());
        for point in &correct.points {
            vec_approx_eq(&point.coords, &correct.points[point.id].coords, 1e-15);
        }
        for cell in &correct.cells {
            assert_eq!(cell.points, correct.cells[cell.id].points);
        }
        mesh.check_all().unwrap();
    }

    #[test]
    fn subdivide_2d_qua16_works() {
        // 29---34----31---28---44---42----40
        //  |               |               |
        // 32   39    38   33   48   47    43
        //  |               |               |
        // 35   36    37   30   45   46    41
        //  |               |               |
        //  3---10-----6----2---23---20----17
        //  |               |               |
        //  7   15    14    9   27   26    22
        //  |               |               |
        // 11   12    13    5   24   25    19
        //  |               |               |
        //  0----4-----8----1---18---21----16
        #[rustfmt::skip]
        let mut block = Block::new(&[
            [0.0, 0.0],
            [3.0, 0.0],
            [3.0, 3.0],
            [0.0, 3.0],
        ]).unwrap();
        let mesh = block.subdivide(GeoKind::Qua16).unwrap();
        let correct = Samples::block_2d_four_qua16();
        assert_eq!(mesh.points.len(), correct.points.len());
        assert_eq!(mesh.cells.len(), correct.cells.len());
        for point in &correct.points {
            vec_approx_eq(&point.coords, &correct.points[point.id].coords, 1e-15);
        }
        for cell in &correct.cells {
            assert_eq!(cell.points, correct.cells[cell.id].points);
        }
        mesh.check_all().unwrap();
    }

    #[test]
    fn subdivide_2d_qua17_works() {
        // 30---38---32---37---29---48---43---47---41
        //  |                   |                   |
        // 39                  36                  46
        //  |                   |                   |
        // 33        34        31        44        42
        //  |                   |                   |
        // 40                  35                  45
        //  |                   |                   |
        //  3---14----6---13----2---28---21---27---18
        //  |                   |                   |
        // 15                  12                  26
        //  |                   |                   |
        //  7         8         5        22        20
        //  |                   |                   |
        // 16                  11                  25
        //  |                   |                   |
        //  0----9----4---10----1---23---19---24---17
        #[rustfmt::skip]
        let mut block = Block::new(&[
            [0.0, 0.0],
            [4.0, 0.0],
            [4.0, 4.0],
            [0.0, 4.0],
        ]).unwrap();
        let mesh = block.subdivide(GeoKind::Qua17).unwrap();
        let correct = Samples::block_2d_four_qua17();
        assert_eq!(format!("{:?}", mesh), format!("{:?}", correct));
        mesh.check_all().unwrap();
    }

    #[test]
    fn subdivide_3d_works() {
        //              18------------------21------------------25
        //              /.                  /.                  /|
        //             / .                 / .                 / |
        //            /  .                /  .                /  |
        //           /   .               /   .               /   |
        //          /    .              /    .              /    |
        //        19------------------20------------------24     |
        //        /.     .            /.     .            /|     |
        //       / .     .           / .     .           / |     |
        //      /  .     .          /  .     .          /  |     |
        //     /   .     .         /   .     .         /   |     |
        //    /    .     .        /    .     .        /    |     |
        //  22==================23==================26     |     |
        //   |     .     .       |     .     .       |     |     |
        //   |     .     .       |     .     .       |     |     |
        //   |     .     4 - - - | - - . - - 7 - - - | - - | - -15
        //   |     .    /.       |     .    /.       |     |    /|
        //   |     .   / .       |     .   / .       |     |   / |
        //   |     .  /  .       |     .  /  .       |     |  /  |
        //   |     . /   .       |     . /   .       |     | /   |
        //   |     ./    .       |     ./    .       |     |/    |
        //   |     5 - - - - - - | - - 6 - - - - - - | - -14     |
        //   |    /.     .       |    /.     .       |    /|     |
        //   |   / .     .       |   / .     .       |   / |     |
        //   |  /  .     .       |  /  .     .       |  /  |     |
        //   | /   .     .       | /   .     .       | /   |     |
        //   |/    .     .       |/    .     .       |/    |     |
        //  10==================11==================17     |     |
        //   |     .     .       |     .     .       |     |     |
        //   |     .     .       |     .     .       |     |     |
        //   |     .     0 - - - | - - . - - 3 - - - | - - | - -13
        //   |     .    /        |     .    /        |     |    /
        //   |     .   /         |     .   /         |     |   /
        //   |     .  /          |     .  /          |     |  /
        //   |     . /           |     . /           |     | /
        //   |     ./            |     ./            |     |/
        //   |     1 - - - - - - | - - 2 - - - - - - | - -12
        //   |    /              |    /              |    /
        //   |   /               |   /               |   /
        //   |  /                |  /                |  /
        //   | /                 | /                 | /
        //   |/                  |/                  |/
        //   8===================9==================16
        #[rustfmt::skip]
        let mut block = Block::new(&[
            [0.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [2.0, 2.0, 0.0],
            [0.0, 2.0, 0.0],
            [0.0, 0.0, 4.0],
            [2.0, 0.0, 4.0],
            [2.0, 2.0, 4.0],
            [0.0, 2.0, 4.0],
        ]).unwrap();
        let mesh = block.subdivide(GeoKind::Hex8).unwrap();
        let correct = Samples::block_3d_eight_hex8();
        assert_eq!(format!("{:?}", mesh), format!("{:?}", correct));
        mesh.check_all().unwrap();
    }

    #[test]
    fn subdivide_3d_o2_works() {
        //              51--------58--------54--------74--------71
        //              /.                  /.                  /|
        //             / .                 / .                 / |
        //           55  .               57  .               73  |
        //           /   .               /   .               /   |
        //          /    .              /    .              /    |
        //        52--------56--------53--------72--------70     |
        //        /.     .            /.     .            /|     |
        //       / .    59           / .    62           / |    76
        //     65  .     .         67  .     .         79  |     |
        //     /   .     .         /   .     .         /   |     |
        //    /    .     .        /    .     .        /    |     |
        //  63========66========64========78========77     |     |
        //   |     .     .       |     .     .       |     |     |
        //   |    60     .       |    61     .       |    75     |
        //   |     .     4 - - - |15 - . - - 7 - - - |41 - | - -35
        //   |     .    /.       |     .    /.       |     |    /|
        //   |     .   / .       |     .   / .       |     |   / |
        //   |     . 12  .       |     . 14  .       |     | 40  |
        //   |     . /   .       |     . /   .       |     | /   |
        //  68     ./    .      69     ./    .      80     |/    |
        //   |     5 - - - -13 - | - - 6 - - - -39 - | - -34     |
        //   |    /.     .       |    /.     .       |    /|     |
        //   |   / .    16       |   / .    19       |   / |    43
        //   | 27  .     .       | 29  .     .       | 49  |     |
        //   | /   .     .       | /   .     .       | /   |     |
        //   |/    .     .       |/    .     .       |/    |     |
        //  22========28========23========48========45     |     |
        //   |     .     .       |     .     .       |     |     |
        //   |    17     .       |    18     .       |    42     |
        //   |     .     0 - - - |11 - . - - 3 - - - |38 - | - -33
        //   |     .    /        |     .    /        |     |    /
        //   |     .   /         |     .   /         |     |   /
        //   |     .  8          |     . 10          |     | 37
        //   |     . /           |     . /           |     | /
        //  30     ./           31     ./           50     |/
        //   |     1 - - - - 9 - | - - 2 - - - -36 - | - -32
        //   |    /              |    /              |    /
        //   |   /               |   /               |   /
        //   | 24                | 26                | 47
        //   | /                 | /                 | /
        //   |/                  |/                  |/
        //  20========25========21========46========44
        #[rustfmt::skip]
        let mut block = Block::new(&[
            [0.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [2.0, 2.0, 0.0],
            [0.0, 2.0, 0.0],
            [0.0, 0.0, 4.0],
            [2.0, 0.0, 4.0],
            [2.0, 2.0, 4.0],
            [0.0, 2.0, 4.0],
        ]).unwrap();
        let mesh = block.subdivide(GeoKind::Hex20).unwrap();
        let correct = Samples::block_3d_eight_hex20();
        assert_eq!(format!("{:?}", mesh), format!("{:?}", correct));
        mesh.check_all().unwrap();
    }

    #[test]
    fn transform_into_ring_works_2d() {
        let mut block = Block::new_square(1.0);
        block.set_ndiv(&[2, 2]).unwrap();
        block
            .set_transform_into_ring(Some(ArgsRing {
                amin: 0.0,
                amax: PI / 2.0,
                rmin: 3.0,
                rmax: 8.0,
                zmin: 0.0,
                zmax: 1.0,
            }))
            .unwrap();
        let mesh = block.subdivide(GeoKind::Qua4).unwrap();
        mesh.check_all().unwrap();
        assert_eq!(mesh.points.len(), 9);
        assert_eq!(mesh.cells.len(), 4);
        assert_eq!(mesh.cells[0].points, &[0, 1, 2, 3]);
        assert_eq!(mesh.cells[1].points, &[1, 4, 5, 2]);
        assert_eq!(mesh.cells[2].points, &[3, 2, 6, 7]);
        assert_eq!(mesh.cells[3].points, &[2, 5, 8, 6]);
        for point in &mesh.points {
            let mut radius = 0.0;
            for i in 0..2 {
                assert!(point.coords[i] >= 0.0);
                assert!(point.coords[i] <= block.args_ring.rmax);
                radius += point.coords[i] * point.coords[i];
            }
            radius = f64::sqrt(radius);
            if [0, 3, 7].contains(&point.id) {
                approx_eq(radius, block.args_ring.rmin, 1e-15);
            }
            if [1, 2, 6].contains(&point.id) {
                approx_eq(radius, (block.args_ring.rmin + block.args_ring.rmax) / 2.0, 1e-17);
            }
            if [4, 5, 8].contains(&point.id) {
                approx_eq(radius, block.args_ring.rmax, 1e-17);
            }
        }
        if SAVE_FIGURE {
            draw_ring(
                &mesh,
                &block.args_ring,
                "/tmp/gemlab/test_block_transform_into_ring_2d.svg",
            );
        }
    }

    #[test]
    fn transform_into_ring_works_2d_qua16() {
        let mut block = Block::new_square(1.0);
        block.set_ndiv(&[2, 2]).unwrap();
        block
            .set_transform_into_ring(Some(ArgsRing {
                amin: 0.0,
                amax: PI / 2.0,
                rmin: 3.0,
                rmax: 8.0,
                zmin: 0.0,
                zmax: 1.0,
            }))
            .unwrap();
        let mesh = block.subdivide(GeoKind::Qua16).unwrap();
        mesh.check_all().unwrap();
        for point in &mesh.points {
            let mut radius = 0.0;
            for i in 0..2 {
                assert!(point.coords[i] >= 0.0);
                assert!(point.coords[i] <= block.args_ring.rmax);
                radius += point.coords[i] * point.coords[i];
            }
            radius = f64::sqrt(radius);
            if [0, 11, 7, 3, 35, 32, 29].contains(&point.id) {
                approx_eq(radius, block.args_ring.rmin, 1e-15);
            }
            if [1, 5, 9, 2, 30, 33, 28].contains(&point.id) {
                approx_eq(radius, (block.args_ring.rmin + block.args_ring.rmax) / 2.0, 1e-17);
            }
            if [16, 19, 22, 17, 41, 43, 40].contains(&point.id) {
                approx_eq(radius, block.args_ring.rmax, 1e-17);
            }
        }
        if SAVE_FIGURE {
            draw_ring(
                &mesh,
                &block.args_ring,
                "/tmp/gemlab/test_block_transform_into_ring_2d_qua16.svg",
            );
        }
    }

    #[test]
    fn transform_into_ring_works_3d() {
        let mut block = Block::new_cube(1.0);
        block.set_ndiv(&[2, 2, 2]).unwrap();
        block
            .set_transform_into_ring(Some(ArgsRing {
                amin: 0.0,
                amax: PI / 2.0,
                rmin: 3.0,
                rmax: 8.0,
                zmin: 0.0,
                zmax: 2.0,
            }))
            .unwrap();
        let mesh = block.subdivide(GeoKind::Hex8).unwrap();
        mesh.check_all().unwrap();
        assert_eq!(mesh.points.len(), 27);
        assert_eq!(mesh.cells.len(), 8);
        assert_eq!(mesh.cells[0].points, &[0, 1, 2, 3, 4, 5, 6, 7]);
        assert_eq!(mesh.cells[7].points, &[6, 11, 17, 14, 20, 23, 26, 24]);
        for point in &mesh.points {
            let mut radius = 0.0;
            for i in 0..2 {
                assert!(point.coords[i] >= 0.0);
                assert!(point.coords[i] <= block.args_ring.rmax);
                radius += point.coords[i] * point.coords[i];
            }
            radius = f64::sqrt(radius);
            if [0, 3, 4, 7, 13, 15, 18, 21, 25].contains(&point.id) {
                approx_eq(radius, block.args_ring.rmin, 1e-15);
            }
            if [1, 2, 5, 6, 12, 14, 19, 20, 24].contains(&point.id) {
                approx_eq(radius, (block.args_ring.rmin + block.args_ring.rmax) / 2.0, 1e-17);
            }
            if [8, 9, 10, 11, 16, 17, 22, 23, 26].contains(&point.id) {
                approx_eq(radius, block.args_ring.rmax, 1e-17);
            }
            if [0, 1, 2, 3, 8, 9, 12, 13, 16].contains(&point.id) {
                approx_eq(point.coords[2], block.args_ring.zmin, 1e-15);
            }
            if [4, 5, 6, 7, 10, 11, 14, 15, 17].contains(&point.id) {
                approx_eq(
                    point.coords[2],
                    (block.args_ring.zmin + block.args_ring.zmax) / 2.0,
                    1e-15,
                );
            }
            if [18, 19, 20, 21, 22, 23, 24, 25, 26].contains(&point.id) {
                approx_eq(point.coords[2], block.args_ring.zmax, 1e-15);
            }
        }
        if SAVE_FIGURE {
            let mut fig = Figure::new();
            fig.figure_size = Some((600.0, 600.0));
            mesh.draw(
                Some(fig),
                "/tmp/gemlab/test_block_transform_into_ring_3d.svg",
                |_, _| {},
            )
            .unwrap();
        }
    }

    #[test]
    fn transform_into_ring_works_3d_hex32() {
        let mut block = Block::new_cube(1.0);
        block.set_ndiv(&[2, 2, 2]).unwrap();
        block
            .set_transform_into_ring(Some(ArgsRing {
                amin: 0.0,
                amax: PI / 2.0,
                rmin: 3.0,
                rmax: 8.0,
                zmin: 0.0,
                zmax: 2.0,
            }))
            .unwrap();
        let mesh = block.subdivide(GeoKind::Hex32).unwrap();
        mesh.check_all().unwrap();
        for point in &mesh.points {
            let mut radius = 0.0;
            for i in 0..2 {
                assert!(point.coords[i] >= 0.0);
                assert!(point.coords[i] <= block.args_ring.rmax);
                radius += point.coords[i] * point.coords[i];
            }
            radius = f64::sqrt(radius);
            if [0, 15, 14, 3, 61, 60, 53, 24, 25, 4, 7, 67, 66, 55, 123, 122, 71, 94].contains(&point.id) {
                approx_eq(radius, block.args_ring.rmin, 1e-15);
            }
            if [90, 91, 118, 62, 56, 63, 57, 119, 5, 27, 1].contains(&point.id) {
                approx_eq(radius, (block.args_ring.rmin + block.args_ring.rmax) / 2.0, 1e-17);
            }
            if [32, 48, 49, 72, 130, 78, 109, 44, 39, 38].contains(&point.id) {
                approx_eq(radius, block.args_ring.rmax, 1e-17);
            }
        }
        if SAVE_FIGURE {
            let mut fig = Figure::new();
            fig.figure_size = Some((600.0, 600.0));
            mesh.draw(
                Some(fig),
                "/tmp/gemlab/test_block_transform_into_ring_3d_hex32.svg",
                |_, _| {},
            )
            .unwrap();
        }
    }

    #[test]
    fn constraints_2d_handles_errors() {
        #[rustfmt::skip]
        let mut block = Block::new(&[
            [0.0, 0.0],
            [2.0, 0.0],
            [2.0, 2.0],
            [0.0, 2.0],
        ]).unwrap();
        block
            .set_edge_constraint(0, Some(Constraint2D::Circle(0.0, 0.0, 1.0)))
            .unwrap();
        assert_eq!(
            block.subdivide(GeoKind::Qua4).err(),
            Some("cannot apply constraint because a point is at the center of the circle")
        );
    }

    #[test]
    fn constraints_2d_works() {
        // block does not touch constraint
        #[rustfmt::skip]
        let mut block = Block::new(&[
            [1.2, 0.0],
            [2.0, 0.0],
            [0.0, 2.0],
            [0.0, 1.2],
        ]).unwrap();
        block.set_ndiv(&[2, 2]).unwrap();

        // circle pushes point
        let ct = Constraint2D::Circle(0.0, 0.0, 1.0);
        block.set_edge_constraint(3, Some(ct)).unwrap();
        let mesh = block.subdivide(GeoKind::Qua4).unwrap();
        for p in [0, 3, 7] {
            let d = point_point_distance(&mesh.points[p].coords, &[0.0, 0.0]).unwrap();
            approx_eq(d, 1.0, 1e-15);
        }
        if SAVE_FIGURE {
            draw(
                mesh,
                &block,
                true,
                "/tmp/gemlab/test_block_constraints_2d_qua4_1.svg",
                |_| {},
            );
        }

        // circle pulls point
        let ct = Constraint2D::Circle(0.0, 0.0, 0.5);
        block.set_edge_constraint(3, Some(ct)).unwrap();
        let mesh = block.subdivide(GeoKind::Qua4).unwrap();
        for p in [0, 3, 7] {
            let d = point_point_distance(&mesh.points[p].coords, &[0.0, 0.0]).unwrap();
            approx_eq(d, 0.5, 1e-15);
        }
        if SAVE_FIGURE {
            draw(
                mesh,
                &block,
                true,
                "/tmp/gemlab/test_block_constraints_2d_qua4_2.svg",
                |_| {},
            );
        }

        // block touches constraint (also we need to move mid nodes)
        #[rustfmt::skip]
        let mut block = Block::new(&[
            [1.0, 0.0],
            [2.0, 0.0],
            [0.0, 2.0],
            [0.0, 1.0],
        ]).unwrap();
        block.set_ndiv(&[2, 2]).unwrap();

        // define circle at negative space with larger radius
        // but touching the points (1.0,0.0) and (0.0,1.0)
        let r = 3.0;
        let beta = f64::asin(SQRT_2 / (2.0 * r));
        let theta = PI / 4.0 - beta;
        let xc = -r * f64::sin(theta);
        let ct = Constraint2D::Circle(xc, xc, r);
        block.set_edge_constraint(3, Some(ct)).unwrap();
        let mesh = block.subdivide(GeoKind::Qua8).unwrap();
        // side 3
        for p in [0, 3, 7, 14, 17] {
            let d = point_point_distance(&mesh.points[p].coords, &[xc, xc]).unwrap();
            approx_eq(d, r, 1e-15);
        }
        // middle nodes
        for (a, mid, b) in [(3, 6, 2)] {
            let xmid = (mesh.points[a].coords[0] + mesh.points[b].coords[0]) / 2.0;
            let ymid = (mesh.points[a].coords[1] + mesh.points[b].coords[1]) / 2.0;
            vec_approx_eq(&mesh.points[mid].coords, &[xmid, ymid], 1e-15);
        }
        if SAVE_FIGURE {
            draw(
                mesh,
                &block,
                true,
                "/tmp/gemlab/test_block_constraints_2d_qua8.svg",
                |_| {},
            );
        }
    }

    #[test]
    fn constraints_2d_multiple_works() {
        let half_l = 2.0;
        #[rustfmt::skip]
        let mut block = Block::new(&[
            [-half_l, -half_l],
            [ half_l, -half_l],
            [ half_l,  half_l],
            [-half_l,  half_l],
        ]).unwrap();
        block.set_ndiv(&[2, 2]).unwrap();
        let r = 8.0;
        let theta = f64::asin(half_l / r);
        let cen_minus = -half_l - r * f64::cos(theta);
        let cen_plus = half_l + r * f64::cos(theta);
        block
            .set_edge_constraint(0, Some(Constraint2D::Circle(0.0, cen_minus, r)))
            .unwrap();
        block
            .set_edge_constraint(1, Some(Constraint2D::Circle(cen_plus, 0.0, r)))
            .unwrap();
        block
            .set_edge_constraint(2, Some(Constraint2D::Circle(0.0, cen_plus, r)))
            .unwrap();
        block
            .set_edge_constraint(3, Some(Constraint2D::Circle(cen_minus, 0.0, r)))
            .unwrap();
        let mesh = block.subdivide(GeoKind::Qua8).unwrap();
        mesh.check_all().unwrap();
        // side 0
        for p in [0, 4, 1, 10, 8] {
            let d = point_point_distance(&mesh.points[p].coords, &[0.0, cen_minus]).unwrap();
            approx_eq(d, r, 1e-15);
        }
        // side 1
        for p in [8, 11, 9, 19, 18] {
            let d = point_point_distance(&mesh.points[p].coords, &[cen_plus, 0.0]).unwrap();
            approx_eq(d, r, 1e-15);
        }
        // side 2
        for p in [14, 16, 13, 20, 18] {
            let d = point_point_distance(&mesh.points[p].coords, &[0.0, cen_plus]).unwrap();
            approx_eq(d, r, 1e-15);
        }
        // side 3
        for p in [0, 7, 3, 17, 14] {
            let d = point_point_distance(&mesh.points[p].coords, &[cen_minus, 0.0]).unwrap();
            approx_eq(d, r, 1e-15);
        }
        // center vertical line
        for p in [5, 2, 15] {
            assert_eq!(mesh.points[p].coords[0], 0.0);
        }
        // horizontal vertical line
        for p in [6, 2, 12] {
            assert_eq!(mesh.points[p].coords[1], 0.0);
        }
        // middle nodes
        for (a, mid, b) in [(3, 6, 2), (2, 12, 9), (1, 5, 2), (2, 15, 13)] {
            let xmid = (mesh.points[a].coords[0] + mesh.points[b].coords[0]) / 2.0;
            let ymid = (mesh.points[a].coords[1] + mesh.points[b].coords[1]) / 2.0;
            vec_approx_eq(&mesh.points[mid].coords, &[xmid, ymid], 1e-15);
        }
        if SAVE_FIGURE {
            draw(
                mesh,
                &block,
                false,
                "/tmp/gemlab/test_block_constraints_2d_multiple.svg",
                |_| {},
            );
        }
    }

    #[test]
    fn constraints_3d_works() {
        let half_l = 2.0;
        #[rustfmt::skip]
        let mut block = Block::new(&[
            [-half_l, -half_l, -half_l],
            [ half_l, -half_l, -half_l],
            [ half_l,  half_l, -half_l],
            [-half_l,  half_l, -half_l],
            [-half_l, -half_l,  half_l],
            [ half_l, -half_l,  half_l],
            [ half_l,  half_l,  half_l],
            [-half_l,  half_l,  half_l],
        ]).unwrap();
        let r = 5.0;
        let theta = f64::asin(half_l / r);
        let cen_minus = -half_l - r * f64::cos(theta);
        let cen_plus = half_l + r * f64::cos(theta);
        block
            .set_face_constraint(0, Some(Constraint3D::CylinderZ(cen_minus, 0.0, r)))
            .unwrap();
        block
            .set_face_constraint(1, Some(Constraint3D::CylinderZ(cen_plus, 0.0, r)))
            .unwrap();
        block
            .set_face_constraint(2, Some(Constraint3D::CylinderZ(0.0, cen_minus, r)))
            .unwrap();
        block
            .set_face_constraint(3, Some(Constraint3D::CylinderZ(0.0, cen_plus, r)))
            .unwrap();
        let mesh = block.subdivide(GeoKind::Hex20).unwrap();
        mesh.check_all().unwrap();
        // corner points
        for p in [0, 20, 33, 44] {
            assert_eq!(mesh.points[p].coords[2], -half_l);
        }
        for p in [51, 63, 77, 71] {
            assert_eq!(mesh.points[p].coords[2], half_l);
        }
        // face 0
        for p in [
            33, 38, 3, 11, 0, // z-min
            43, 19, 16, //
            35, 41, 7, 15, 4, // z-mid
            76, 62, 59, //
            71, 74, 54, 58, 51, // z-max
        ] {
            let x = &[mesh.points[p].coords[0], mesh.points[p].coords[1]];
            let d = point_point_distance(x, &[cen_minus, 0.0]).unwrap();
            approx_eq(d, r, 1e-15);
        }
        // face 1
        for p in [
            20, 25, 21, 46, 44, // z-min
            30, 31, 50, //
            22, 28, 23, 48, 45, // z-mid
            68, 69, 80, //
            63, 66, 64, 78, 77, // z-max
        ] {
            let x = &[mesh.points[p].coords[0], mesh.points[p].coords[1]];
            let d = point_point_distance(x, &[cen_plus, 0.0]).unwrap();
            approx_eq(d, r, 1e-15);
        }
        // face 2
        for p in [
            0, 8, 1, 24, 20, // z-min
            16, 17, 30, //
            4, 12, 5, 27, 22, // z-mid
            59, 60, 68, //
            51, 55, 52, 65, 68, // z-max
        ] {
            let x = &[mesh.points[p].coords[0], mesh.points[p].coords[1]];
            let d = point_point_distance(x, &[0.0, cen_minus]).unwrap();
            approx_eq(d, r, 1e-15);
        }
        // face 3
        for p in [
            44, 47, 32, 37, 33, // z-min
            50, 42, 43, //
            45, 49, 34, 40, 35, // z-mid
            80, 75, 76, //
            77, 79, 70, 73, 71, // z-max
        ] {
            let x = &[mesh.points[p].coords[0], mesh.points[p].coords[1]];
            let d = point_point_distance(x, &[0.0, cen_plus]).unwrap();
            approx_eq(d, r, 1e-15);
        }
        // vertical line at center
        for p in [2, 18, 6, 61, 53] {
            assert_eq!(mesh.points[p].coords[0], 0.0);
            assert_eq!(mesh.points[p].coords[1], 0.0);
        }
        // middle nodes
        for (a, mid, b) in [
            (2, 10, 3),
            (2, 26, 21),
            (2, 9, 1),
            (2, 36, 32), // z-min
            (6, 13, 5),
            (6, 39, 34),
            (6, 14, 7),
            (6, 29, 23), // z-mid
            (53, 56, 52),
            (53, 72, 70),
            (53, 57, 54),
            (53, 67, 64), // z-max
        ] {
            let xmid = (mesh.points[a].coords[0] + mesh.points[b].coords[0]) / 2.0;
            let ymid = (mesh.points[a].coords[1] + mesh.points[b].coords[1]) / 2.0;
            let zmid = (mesh.points[a].coords[2] + mesh.points[b].coords[2]) / 2.0;
            vec_approx_eq(&mesh.points[mid].coords, &[xmid, ymid, zmid], 1e-15);
        }
        if SAVE_FIGURE {
            draw(
                mesh,
                &block,
                false,
                "/tmp/gemlab/test_block_constraints_3d.svg",
                |plot| {
                    let mut surf = Surface::new();
                    const NP: usize = 81;
                    surf.set_solid_color("#ff000020");
                    surf.draw_cylinder(&[cen_minus, 0.0, -half_l], &[cen_minus, 0.0, half_l], r, 5, NP)
                        .unwrap();
                    surf.set_solid_color("#00ff0020");
                    surf.draw_cylinder(&[cen_plus, 0.0, -half_l], &[cen_plus, 0.0, half_l], r, 5, NP)
                        .unwrap();
                    surf.set_solid_color("#0000ff20");
                    surf.draw_cylinder(&[0.0, cen_minus, -half_l], &[0.0, cen_minus, half_l], r, 5, NP)
                        .unwrap();
                    surf.set_solid_color("#ff00ff20");
                    surf.draw_cylinder(&[0.0, cen_plus, -half_l], &[0.0, cen_plus, half_l], r, 5, NP)
                        .unwrap();
                    plot.add(&surf);
                    plot.set_range_3d(-half_l, half_l, -half_l, half_l, -half_l, half_l);
                },
            );
        }
    }

    #[test]
    fn constraints_3d_works_cylinder_xy() {
        let half_l = 2.0;
        #[rustfmt::skip]
        let mut block = Block::new(&[
            [-half_l, -half_l, -half_l],
            [ half_l, -half_l, -half_l],
            [ half_l,  half_l, -half_l],
            [-half_l,  half_l, -half_l],
            [-half_l, -half_l,  half_l],
            [ half_l, -half_l,  half_l],
            [ half_l,  half_l,  half_l],
            [-half_l,  half_l,  half_l],
        ]).unwrap();
        let r = 5.0;
        let theta = f64::asin(half_l / r);
        let cen_minus = -half_l - r * f64::cos(theta);
        let cen_plus = half_l + r * f64::cos(theta);
        block
            .set_face_constraint(4, Some(Constraint3D::CylinderX(0.0, cen_minus, r)))
            .unwrap();
        block
            .set_face_constraint(5, Some(Constraint3D::CylinderY(0.0, cen_plus, r)))
            .unwrap();
        let mesh = block.subdivide(GeoKind::Hex20).unwrap();
        mesh.check_all().unwrap();
        // corner points
        for p in [0, 20, 33, 44] {
            assert_eq!(mesh.points[p].coords[2], -half_l);
        }
        for p in [51, 63, 77, 71] {
            assert_eq!(mesh.points[p].coords[2], half_l);
        }
        // face 4
        for p in [
            0, 11, 3, 38, 33, // x-min
            8, 10, 37, //
            1, 9, 2, 36, 32, // x-mid
            24, 26, 47, //
            20, 25, 21, 46, 44, // x-max
        ] {
            let x = &[mesh.points[p].coords[1], mesh.points[p].coords[2]];
            let d = point_point_distance(x, &[0.0, cen_minus]).unwrap();
            approx_eq(d, r, 1e-15);
        }
        // face 3
        for p in [
            51, 55, 52, 65, 63, // y-min
            58, 56, 66, //
            54, 57, 53, 67, 64, // y-mid
            74, 72, 78, //
            71, 73, 70, 79, 77, // y-max
        ] {
            let x = &[mesh.points[p].coords[0], mesh.points[p].coords[2]];
            let d = point_point_distance(x, &[0.0, cen_plus]).unwrap();
            approx_eq(d, r, 1e-15);
        }
        // vertical line at center
        for p in [2, 18, 6, 61, 53] {
            assert_eq!(mesh.points[p].coords[0], 0.0);
            assert_eq!(mesh.points[p].coords[1], 0.0);
        }
        // middle nodes
        for (a, mid, b) in [
            (3, 19, 7),
            (2, 18, 6),
            (21, 31, 23), // z-min
            (5, 60, 52),
            (6, 61, 53),
            (34, 75, 70), // z-max
        ] {
            let xmid = (mesh.points[a].coords[0] + mesh.points[b].coords[0]) / 2.0;
            let ymid = (mesh.points[a].coords[1] + mesh.points[b].coords[1]) / 2.0;
            let zmid = (mesh.points[a].coords[2] + mesh.points[b].coords[2]) / 2.0;
            vec_approx_eq(&mesh.points[mid].coords, &[xmid, ymid, zmid], 1e-15);
        }
        if SAVE_FIGURE {
            draw(
                mesh,
                &block,
                false,
                "/tmp/gemlab/test_block_constraints_3d_cylinder_xy.svg",
                |plot| {
                    let mut surf = Surface::new();
                    const NP: usize = 81;
                    surf.set_solid_color("#ff000020");
                    surf.draw_cylinder(&[-half_l, 0.0, cen_minus], &[half_l, 0.0, cen_minus], r, 5, NP)
                        .unwrap();
                    surf.set_solid_color("#00ff0020");
                    surf.draw_cylinder(&[0.0, -half_l, cen_plus], &[0.0, half_l, cen_plus], r, 5, NP)
                        .unwrap();
                    plot.add(&surf);
                    plot.set_range_3d(-half_l, half_l, -half_l, half_l, -half_l, half_l);
                },
            );
        }
    }

    #[test]
    fn constraints_2d_handles_imprecision() {
        let mut block = Block::new_square(6.0);
        block.set_ndiv(&[3, 3]).unwrap();
        let ct = Constraint2D::Circle(0.0, 0.0, 6.0);
        block.set_edge_constraint(1, Some(ct.clone())).unwrap();
        block.set_edge_constraint(2, Some(ct.clone())).unwrap();
        let mesh = block.subdivide(GeoKind::Qua4).unwrap();
        for p in [6, 7, 11, 15] {
            let d = point_point_distance(&mesh.points[p].coords, &[0.0, 0.0]).unwrap();
            approx_eq(d, 6.0, 1e-15);
        }
        for p in [13, 12, 14] {
            let d = point_point_distance(&mesh.points[p].coords, &[0.0, 0.0]).unwrap();
            approx_eq(d, 6.0, 1e-15);
        }
        if SAVE_FIGURE {
            draw(
                mesh,
                &block,
                true,
                "/tmp/gemlab/test_block_constraints_imprecision_2d.svg",
                |_| {},
            );
        }
    }

    #[test]
    fn constraints_handles_interior_nodes_qua8() {
        let radius = 6.0;
        let m = radius / 2.0;
        let n = radius / SQRT_2;
        let p = 1.15 * m / SQRT_2;
        let mut block = Block::new(&[[m, 0.0], [radius, 0.0], [n, n], [p, p]]).unwrap();
        block.set_ndiv(&[3, 3]).unwrap();
        let ct = Constraint2D::Circle(0.0, 0.0, radius);
        block.set_edge_constraint(1, Some(ct)).unwrap();
        let mesh = block.subdivide(GeoKind::Qua8).unwrap();
        for p in [13, 16, 14, 27, 26, 38, 37] {
            let d = point_point_distance(&mesh.points[p].coords, &[0.0, 0.0]).unwrap();
            approx_eq(d, radius, 1e-15);
        }
        for (a, mid, b) in [(23, 28, 26), (9, 17, 14)] {
            let xmid = (mesh.points[a].coords[0] + mesh.points[b].coords[0]) / 2.0;
            let ymid = (mesh.points[a].coords[1] + mesh.points[b].coords[1]) / 2.0;
            vec_approx_eq(&mesh.points[mid].coords, &[xmid, ymid], 1e-15);
        }
        if SAVE_FIGURE {
            draw(
                mesh,
                &block,
                true,
                "/tmp/gemlab/test_block_constraints_interior_nodes_qua8.svg",
                |_| {},
            );
        }
    }

    #[test]
    fn constraints_handles_interior_nodes_qua9() {
        let radius = 6.0;
        let m = radius / 2.0;
        let n = radius / SQRT_2;
        let p = 1.15 * m / SQRT_2;
        let mut block = Block::new(&[[m, 0.0], [radius, 0.0], [n, n], [p, p]]).unwrap();
        block.set_ndiv(&[3, 3]).unwrap();
        let ct = Constraint2D::Circle(0.0, 0.0, radius);
        block.set_edge_constraint(1, Some(ct)).unwrap();
        let mesh = block.subdivide(GeoKind::Qua9).unwrap();
        for p in [15, 18, 16, 32, 31, 46, 45] {
            let d = point_point_distance(&mesh.points[p].coords, &[0.0, 0.0]).unwrap();
            approx_eq(d, radius, 1e-15);
        }
        for (a, mid, b) in [(42, 48, 46), (27, 33, 31), (28, 34, 32), (10, 19, 16), (12, 20, 18)] {
            let xmid = (mesh.points[a].coords[0] + mesh.points[b].coords[0]) / 2.0;
            let ymid = (mesh.points[a].coords[1] + mesh.points[b].coords[1]) / 2.0;
            vec_approx_eq(&mesh.points[mid].coords, &[xmid, ymid], 1e-15);
        }
        if SAVE_FIGURE {
            draw(
                mesh,
                &block,
                true,
                "/tmp/gemlab/test_block_constraints_interior_nodes_qua9.svg",
                |_| {},
            );
        }
    }

    #[test]
    fn constraints_handles_interior_nodes_qua16() {
        let radius = 6.0;
        let m = radius / 2.0;
        let n = radius / SQRT_2;
        let p = 1.15 * m / SQRT_2;
        let mut block = Block::new(&[[m, 0.0], [radius, 0.0], [n, n], [p, p]]).unwrap();
        block.set_ndiv(&[3, 3]).unwrap();
        let ct = Constraint2D::Circle(0.0, 0.0, radius);
        block.set_edge_constraint(1, Some(ct)).unwrap();
        let mesh = block.subdivide(GeoKind::Qua16).unwrap();
        for (a, c, d, b) in [
            (85, 99, 98, 94),
            (83, 96, 97, 92),
            (52, 65, 63, 61),
            (55, 69, 68, 64),
            (53, 66, 67, 62),
            (17, 35, 32, 29),
            (22, 39, 38, 34),
            (19, 36, 37, 31),
        ] {
            let lx = mesh.points[b].coords[0] - mesh.points[a].coords[0];
            let ly = mesh.points[b].coords[1] - mesh.points[a].coords[1];
            let xc = mesh.points[a].coords[0] + lx / 3.0;
            let yc = mesh.points[a].coords[1] + ly / 3.0;
            let xd = mesh.points[a].coords[0] + 2.0 * lx / 3.0;
            let yd = mesh.points[a].coords[1] + 2.0 * ly / 3.0;
            vec_approx_eq(&mesh.points[c].coords, &[xc, yc], 1e-14);
            vec_approx_eq(&mesh.points[d].coords, &[xd, yd], 1e-15);
        }
        if SAVE_FIGURE {
            draw(
                mesh,
                &block,
                true,
                "/tmp/gemlab/test_block_constraints_interior_nodes_qua16.svg",
                |_| {},
            );
        }
    }

    #[test]
    fn constraints_handles_interior_nodes_qua17() {
        let radius = 6.0;
        let m = radius / 2.0;
        let n = radius / SQRT_2;
        let p = 1.15 * m / SQRT_2;
        let mut block = Block::new(&[[m, 0.0], [radius, 0.0], [n, n], [p, p]]).unwrap();
        block.set_ndiv(&[3, 3]).unwrap();
        let ct = Constraint2D::Circle(0.0, 0.0, radius);
        block.set_edge_constraint(1, Some(ct)).unwrap();
        let mesh = block.subdivide(GeoKind::Qua17).unwrap();
        for (a, mid, b) in [(82, 92, 90), (54, 64, 62), (20, 34, 32)] {
            let xmid = (mesh.points[a].coords[0] + mesh.points[b].coords[0]) / 2.0;
            let ymid = (mesh.points[a].coords[1] + mesh.points[b].coords[1]) / 2.0;
            vec_approx_eq(&mesh.points[mid].coords, &[xmid, ymid], 1e-15);
        }
        if SAVE_FIGURE {
            draw(
                mesh,
                &block,
                true,
                "/tmp/gemlab/test_block_constraints_interior_nodes_qua17.svg",
                |_| {},
            );
        }
    }

    #[test]
    fn constraints_handles_interior_nodes_hex20() {
        let radius = 6.0;
        let m = radius / 2.0;
        let n = radius / SQRT_2;
        let p = 1.15 * m / SQRT_2;
        #[rustfmt::skip]
        let mut block = Block::new(&[
            [     m, 0.0, 0.0],
            [radius, 0.0, 0.0],
            [     n,   n, 0.0],
            [     p,   p, 0.0],
            [     m, 0.0, 1.0],
            [radius, 0.0, 1.0],
            [     n,   n, 1.0],
            [     p,   p, 1.0],
        ]).unwrap();
        block.set_ndiv(&[2, 2, 1]).unwrap();
        let ct = Constraint3D::CylinderZ(0.0, 0.0, radius);
        block.set_face_constraint(1, Some(ct)).unwrap();
        let mesh = block.subdivide(GeoKind::Hex20).unwrap();
        for p in [
            20, 25, 21, 46, 44, // z-min
            30, 31, 50, // z-mid
            22, 28, 23, 48, 45, // z-max
        ] {
            let d = point_point_distance(&mesh.points[p].coords, &[0.0, 0.0, mesh.points[p].coords[2]]).unwrap();
            approx_eq(d, radius, 1e-15);
        }
        for (a, mid, b) in [(2, 26, 21), (6, 29, 23)] {
            let xmid = (mesh.points[a].coords[0] + mesh.points[b].coords[0]) / 2.0;
            let ymid = (mesh.points[a].coords[1] + mesh.points[b].coords[1]) / 2.0;
            let zmid = (mesh.points[a].coords[2] + mesh.points[b].coords[2]) / 2.0;
            vec_approx_eq(&mesh.points[mid].coords, &[xmid, ymid, zmid], 1e-15);
        }
        if SAVE_FIGURE {
            draw(
                mesh,
                &block,
                true,
                "/tmp/gemlab/test_block_constraints_interior_nodes_hex20.svg",
                |_| {},
            );
        }
    }
}
