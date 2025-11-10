use super::primitives::AxisAlignedBoundingBox;
use super::{Ray, Vec3D};

pub trait Surface {
    fn hit(&self, ray: &Ray, t_min: f32, t_max: f32) -> Option<(f32, Vec3D)>;
    fn hit_bool(&self, ray: &Ray, t_min: f32, t_max: f32) -> bool {
        self.hit(ray, t_min, t_max).is_some()
    }
}

pub trait Bounded {
    fn bounding_box(&self) -> AxisAlignedBoundingBox;
}
