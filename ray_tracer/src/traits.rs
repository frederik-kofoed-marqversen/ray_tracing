use fastrand::Rng;

use super::primitives::AxisAlignedBoundingBox;
use super::{Ray, Vec3D};

pub trait Surface {
    fn hit(&self, ray: &Ray, t_min: f32, t_max: f32) -> Option<(f32, Vec3D)>;
}

pub trait Bounded {
    fn bounding_box(&self) -> AxisAlignedBoundingBox;
}

pub trait Light {
    fn sample_light(&self, ref_point: Vec3D, rng: &mut Rng) -> Option<(Vec3D, Vec3D, f32)>;
}
