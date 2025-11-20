use std::collections::{HashMap, VecDeque};

use super::bvh::BoundingVolumeHierarchy as BVH;
use super::primitives::{triangle_intersection, AxisAlignedBoundingBox};
use super::traits::{Bounded, Surface};
use super::AffineTransform;
use super::{Ray, Vec3D};
use std::rc::Rc;

/// TriangleMesh is defined by a array of vertices. Each triangle is defined by three indices into
/// that array, one for each corner of the triangle.
/// TriangleMesh itself does not implement Surface or Bounded. Instead, a BVH over the mesh can be
/// built using TriangleRef primitives, which reference individual triangles within the mesh.
#[derive(Debug, Clone)]
pub struct TriangleMesh {
    pub vertices: Vec<Vec3D>,
    pub triangles: Vec<[usize; 3]>,
}

/// TriangleRef is a reference to a single triangle within a TriangleMesh.
/// It implements Surface and Bounded by looking up the relevant data from the mesh.
#[derive(Debug, Clone)]
pub struct TriangleRef {
    mesh: Rc<TriangleMesh>,
    index: usize,
}

impl Surface for TriangleRef {
    fn hit(&self, ray: &Ray, t_min: f32, t_max: f32) -> Option<(f32, Vec3D)> {
        let tri = self.mesh.triangles[self.index];
        let v0 = self.mesh.vertices[tri[0]];
        let v1 = self.mesh.vertices[tri[1]];
        let v2 = self.mesh.vertices[tri[2]];

        triangle_intersection(v0, v1, v2, ray, t_min, t_max)
    }
}

impl Bounded for TriangleRef {
    #[inline]
    fn bounding_box(&self) -> AxisAlignedBoundingBox {
        let tri = self.mesh.triangles[self.index];
        AxisAlignedBoundingBox::empty()
            .grow(self.mesh.vertices[tri[0]])
            .grow(self.mesh.vertices[tri[1]])
            .grow(self.mesh.vertices[tri[2]])
    }

    #[inline]
    fn centroid_axis(&self, axis: usize) -> f32 {
        let tri = self.mesh.triangles[self.index];
        (self.mesh.vertices[tri[0]][axis]
            + self.mesh.vertices[tri[1]][axis]
            + self.mesh.vertices[tri[2]][axis])
            / 3.0
    }

    #[inline]
    fn centroid(&self) -> Vec3D {
        let tri = self.mesh.triangles[self.index];
        (self.mesh.vertices[tri[0]] + self.mesh.vertices[tri[1]] + self.mesh.vertices[tri[2]]) / 3.0
    }
}

impl TriangleMesh {
    #[inline]
    fn get_edges(&self, triangle: usize) -> [(usize, usize); 3] {
        let vertices = &self.triangles[triangle];
        [0, 1, 2].map(|i| (vertices[i], vertices[(i + 1) % 3]))
    }

    pub fn build_bvh(mesh: &Rc<TriangleMesh>) -> Rc<BVH<TriangleRef>> {
        let prims: Vec<TriangleRef> = (0..mesh.triangles.len())
            .map(|i| TriangleRef {
                mesh: Rc::clone(mesh),
                index: i,
            })
            .collect();
        Rc::new(BVH::build(prims))
    }

    pub fn apply_transform(&mut self, transform: AffineTransform) {
        self.vertices
            .iter_mut()
            .for_each(|v| *v = transform.transform_point(*v));
    }

    pub fn area_centroid(&self) -> Vec3D {
        // Average centroid weighted by area.
        let mut total_area = 0.0;
        let mut area_centroid = Vec3D::ZERO;
        for triangle in &self.triangles {
            let vertices = triangle.map(|i| self.vertices[i]);
            let area = area_weighted_normal(vertices).norm();
            let centroid = vertices.iter().sum::<Vec3D>() / 3.0;
            total_area += area;
            area_centroid += area * centroid;
        }
        return area_centroid / total_area;
    }

    pub fn volume_centroid(&self) -> Vec3D {
        // Average centroid weighted by volume.
        let mut total_volume = 0.0;
        let mut volume_centroid = Vec3D::ZERO;
        // For each triangle we compute the signed volume weighted centroid of the
        // tetrahedron defined by adding the origin point (0, 0, 0) to the triangle.
        for triangle in &self.triangles {
            let vertices = triangle.map(|i| self.vertices[i]);
            let centroid = vertices.iter().sum::<Vec3D>() / 4.0;
            let signed_volume =
                Vec3D::dot(vertices[0], Vec3D::cross(vertices[1], vertices[2])) / 6.0;
            total_volume += signed_volume;
            volume_centroid += signed_volume * centroid;
        }
        return volume_centroid / total_volume;
    }

    /// Attempts to fix the triangle mesh to ensure a globally consistent normal vector field
    /// (i.e., consistent face orientation). For closed surfaces, the normal vector field
    /// points outwards. This is guaranteed by enforcing the following rules:
    ///
    /// 1. Mesh is a manifold, meaning that each edge is shared by at most two triangles:
    ///    - One triangle for boundary edges.
    ///    - Two triangles for internal edges.
    ///
    /// 2. The vertex ordering of neighboring triangles ensures that they traverse
    ///    their shared edge in opposite directions.
    ///
    /// 3. The signed volume of a closed surface as determined by the divergence theorem
    ///    should be positive.
    ///
    /// Returns Ok(()) if mesh was successfully fixed. Can fail in one of two ways. Either the
    /// mesh does not define a manifold or the mesh defines a non-orientable surface.
    pub fn try_fix_mesh(&mut self) -> Result<(), &'static str> {
        // HashMap mapping undirected edges to triangles that share that edge.
        let mut edge_map: HashMap<(usize, usize), Vec<usize>> =
            HashMap::with_capacity(3 * self.triangles.len());
        for triangle in 0..self.triangles.len() {
            for edge in self.get_edges(triangle) {
                let edge = to_undirected(edge);

                let entry = edge_map.entry(edge).or_insert(Vec::new());
                entry.push(triangle);
                if entry.len() > 2 {
                    // More than two triangles share this edge.
                    return Err("Mesh does not define a manifold.");
                }
            }
        }

        let mut is_closed = true;
        let mut queue = VecDeque::new();
        let mut checked = vec![false; self.triangles.len()];
        // Arbitrarily set first triangle as correct orientation
        queue.push_back(0);
        checked[0] = true;
        // Using BFS, fix triangles such that all edges are consistent with this orientation.
        while let Some(triangle) = queue.pop_front() {
            for edge in self.get_edges(triangle) {
                // Each undirected edge need only be checked once.
                if let Some(tris) = edge_map.remove(&to_undirected(edge)) {
                    let neighbour = tris.into_iter().find(|other| other != &triangle);
                    if neighbour.is_none() {
                        // Edge is a boundary => no neighbour to check
                        is_closed = false;
                        continue;
                    }
                    let neighbour = neighbour.unwrap();

                    // If neighbour shares the flipped edge, then orientation is consistent
                    let are_consistent = self.get_edges(neighbour).contains(&(edge.1, edge.0));
                    match (are_consistent, checked[neighbour]) {
                        (true, true) => {
                            // Neighbour already checked and correctly oriented
                        }
                        (false, true) => {
                            // Neighbour needs to be flipped but is already checked => contradiction
                            return Err("Mesh is not orientable");
                        }
                        (true, false) => {
                            // Add neighbour to queue
                            queue.push_back(neighbour);
                            checked[neighbour] = true;
                        }
                        (false, false) => {
                            // Flip orientation of neighbour then add it to queue
                            self.triangles[neighbour].swap(0, 1);
                            queue.push_back(neighbour);
                            checked[neighbour] = true;
                        }
                    }
                }
            }
        }

        // At this point the mesh should have a consistent orientation.
        if is_closed {
            self.fix_global_orientation();
        }

        return Ok(());
    }

    /// This works by observing that by the divergence theorem the signed volume
    /// of a closed surface can be computed as
    ///     ∫ (x * n_x) dA = ±V
    /// where x is the x coordinate and n_x is the x-coordinate of the unit normal
    /// vector. For a triangulated surface, the surface integral can be evaluated
    /// as a sum over triangles. Over a single triangle the integral evaluates to
    ///     n_x * A * <x>
    /// where A is the area of the triangle and <x> the average x coordinate.
    /// Notice that n_x * A is actually just the x-coordinate of the area weighted
    /// normal vector. Final simplification is done by the fact we only care about
    /// the sign, and so can discard scale factors.
    fn fix_global_orientation(&mut self) {
        // Compute signed volume assuming consistent normal vector orientation.
        let mut volume = 0.0;
        for triangle in &self.triangles {
            let vertices = triangle.map(|i| self.vertices[i]);
            let nx = area_weighted_normal(vertices).x;
            let avg_x = vertices.iter().map(|v| v.x).sum::<f32>();
            volume += nx * avg_x;
        }
        // Flip global orientation
        if volume.is_sign_negative() {
            for triangle in self.triangles.iter_mut() {
                triangle.swap(0, 1);
            }
        }
    }
}

#[inline]
fn area_weighted_normal(triangle: [Vec3D; 3]) -> Vec3D {
    let e1 = triangle[1] - triangle[0];
    let e2 = triangle[2] - triangle[0];
    return 0.5 * Vec3D::cross(e1, e2);
}

#[inline]
fn to_undirected(edge: (usize, usize)) -> (usize, usize) {
    (edge.0.min(edge.1), edge.0.max(edge.1))
}

#[cfg(test)]
mod tests {
    use std::iter::zip;

    use super::*;

    fn mesh_eq(mesh1: &TriangleMesh, mesh2: &TriangleMesh) -> bool {
        for (tri1, tri2) in zip(&mesh1.triangles, &mesh2.triangles) {
            let mut are_equal = false;
            for i in 0..3 {
                let cycle = [tri1[i], tri1[(i + 1) % 3], tri1[(i + 2) % 3]];
                are_equal ^= &cycle == tri2;
            }
            if !are_equal {
                return false;
            }
        }
        return true;
    }

    #[test]
    fn fix_mesh() {
        // use super::super::stl::read_stl;
        // let mut mesh: TriangleMesh = read_stl("../src/meshes/Baby_Yoda.stl").unwrap();

        // Properly oriented mesh
        let original = TriangleMesh {
            vertices: vec![
                Vec3D::new(1.0, 1.0, 1.0),
                Vec3D::new(1.0, -1.0, -1.0),
                Vec3D::new(-1.0, 1.0, -1.0),
                Vec3D::new(-1.0, -1.0, 1.0),
            ],
            triangles: vec![[0, 1, 2], [1, 3, 2], [0, 3, 1], [0, 2, 3]],
        };

        // Make non-consistent clone. By flipping the orientation of the first
        // triangle, we also test that the global orientation is fixed.
        let mut mesh = TriangleMesh {
            vertices: original.vertices.clone(),
            triangles: original.triangles.clone(),
        };
        mesh.triangles[0].swap(0, 1);
        mesh.triangles[1].swap(2, 1);

        assert!(!mesh_eq(&original, &mesh));

        // Try fix
        let result = mesh.try_fix_mesh();
        assert!(result.is_ok());
        assert!(mesh_eq(&original, &mesh));
    }
}
