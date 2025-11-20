use crate::geometry::AxisAlignedBoundingBox;
use crate::Ray;
use crate::Vec3D;

pub trait Surface {
    fn hit(&self, ray: &Ray, t_min: f32, t_max: f32) -> Option<(f32, Vec3D)>;
    /// Should be overridden if a more efficient implementation is possible
    fn hit_bool(&self, ray: &Ray, t_min: f32, t_max: f32) -> bool {
        self.hit(ray, t_min, t_max).is_some()
    }
}

pub trait Bounded {
    fn bounding_box(&self) -> AxisAlignedBoundingBox;
    // Should be overridden if a more efficient/precise implementation is possible
    fn centroid(&self) -> Vec3D {
        self.bounding_box().centroid()
    }
    fn centroid_axis(&self, axis: usize) -> f32 {
        self.centroid()[axis]
    }
}
