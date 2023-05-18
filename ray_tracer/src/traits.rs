use super::{Ray, Vec3D};
use super::primitives::AxisAlignedBoundingBox;

pub trait Surface {
    fn hit(&self, ray: &Ray, t_min: f32, t_max: f32) -> Option<(f32, Vec3D)>;
}

pub trait Bounded {
    fn bounding_box(&self) -> AxisAlignedBoundingBox;
}

pub trait BoundedSurface: Bounded + Surface {}