use super::primitives::AxisAlignedBoundingBox;
use super::traits::{Bounded, Surface};
use super::{Point3D, Ray, Vec3D};
use std::rc::Rc;

/// TriangleMesh is defined by a array of vertices. Each triangle is defined by three indices into
/// that array, one for each corner of the triangle. TriangleMesh does not implement Surface, and
/// so is basically just a datastructure used to hold triangle data.
pub struct TriangleMesh {
    pub vertices: Vec<Point3D>,
    pub triangles: Vec<[usize; 3]>,
}

impl TriangleMesh {
    pub fn triangles(self: &Rc<Self>) -> Vec<Triangle> {
        return (0..self.triangles.len())
            .map(move |index| Triangle {
                mesh: Rc::clone(self),
                index,
            })
            .collect();
    }
}

/// Triangle's do not live by themselves, but refer to a TriangleMesh and themselves only store
/// their corresponding index into that TriangleMesh's triangle array.
pub struct Triangle {
    mesh: Rc<TriangleMesh>,
    index: usize,
}

impl Triangle {
    pub fn vertices(&self) -> [&Point3D; 3] {
        let indices = self.mesh.triangles[self.index];
        return indices.map(|index| &self.mesh.vertices[index]);
    }
}

impl Surface for Triangle {
    fn hit(&self, ray: &Ray, t_min: f32, t_max: f32) -> Option<(f32, Vec3D)> {
        // Möller-Trumbore intersection algorithm. More or less copied from Wikipedia
        let vertices = self.vertices();
        let e1 = vertices[1] - vertices[0];
        let e2 = vertices[2] - vertices[0];
    
        let ray_cross_e2 = Vec3D::cross(ray.direction, e2);
        let det = Vec3D::dot(e1, ray_cross_e2);
    
        if det > -f32::EPSILON && det < f32::EPSILON {
            return None; // Ray is parallel to triangle.
        }
    
        let inv_det = 1.0 / det;
        let s = ray.origin - vertices[0];
        let u = inv_det * Vec3D::dot(s, ray_cross_e2);
        if u < 0.0 || u > 1.0 {
            // Ray misses triangle in u-dir
            return None;
        }
    
        let s_cross_e1 = Vec3D::cross(s, e1);
        let v = inv_det * Vec3D::dot(ray.direction, s_cross_e1);
        if v < 0.0 || u + v > 1.0 {
            // Ray misses triangle in v-dir
            return None;
        }
        
        let t = inv_det * Vec3D::dot(e2, s_cross_e1);
        if t > t_min && t < t_max {
            // Hit!
            let normal = Vec3D::cross(e1, e2).normalise();
            return Some((t, normal));
        }
        else {
            // Ray hits outside of given interval
            return None;
        }
    }
}

impl Bounded for Triangle {
    fn bounding_box(&self) -> AxisAlignedBoundingBox {
        let mut upper = Vec3D::ZERO;
        let mut lower = Vec3D::ZERO;
        let vertices = self.vertices();
        for i in 0..3 {
            upper[i] = vertices
                .iter()
                .fold(f32::NEG_INFINITY, |max, p| max.max(p[i]));
            lower[i] = vertices.iter().fold(f32::INFINITY, |min, p| min.min(p[i]));
        }
        return AxisAlignedBoundingBox { upper, lower };
    }
}
