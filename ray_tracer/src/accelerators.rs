use super::primitives::AxisAlignedBoundingBox;
use super::traits::*;
use super::{Ray, Vec3D};

pub type BVH = BoundingVolumeHierarchy;
pub struct BoundingVolumeHierarchy {
    pub node1: Box<dyn BoundedSurface>,
    pub node2: Box<dyn BoundedSurface>,
    bounding_box: AxisAlignedBoundingBox,
}

impl BoundingVolumeHierarchy {
    pub fn build(mut objects: Vec<Box<dyn BoundedSurface>>) -> Box<dyn BoundedSurface> {
        match objects.len() {
            0 => panic!("zero objects given!"),
            1 => return objects.pop().unwrap(),
            2 => {
                let bounding_box = AxisAlignedBoundingBox::combine(
                    &objects[0].bounding_box(),
                    &objects[1].bounding_box(),
                );
                let result: Box<dyn BoundedSurface> = Box::new(Self {
                    node1: objects.pop().unwrap(),
                    node2: objects.pop().unwrap(),
                    bounding_box,
                });
                return result;
            },
            _ => {
                let (left, right) = center_partition(objects);
                let node1 = Self::build(left);
                let node2 = Self::build(right);
                let bounding_box =
                    AxisAlignedBoundingBox::combine(&node1.bounding_box(), &node2.bounding_box());
                let result: Box<dyn BoundedSurface> = Box::new(Self {
                    node1,
                    node2,
                    bounding_box,
                });
                return result;
            }
        }
    }
}

impl Bounded for BoundingVolumeHierarchy {
    fn bounding_box(&self) -> AxisAlignedBoundingBox {
        self.bounding_box.clone()
    }
}

impl Surface for BoundingVolumeHierarchy {
    fn hit(&self, ray: &Ray, t_min: f32, t_max: f32) -> Option<(f32, Vec3D)> {
        let mut node1 = &self.node1;
        let mut node2 = &self.node2;

        // First check if bounding box of either node is hit.
        let hit1 = node1.bounding_box().hit(ray, t_min, t_max).map(|(t, _)| t);
        let hit2 = node2.bounding_box().hit(ray, t_min, t_max).map(|(t, _)| t);
        let (t1, t2) = match (hit1, hit2) {
            (None, None) => return None,
            (Some(_), None) => return node1.hit(ray, t_min, t_max),
            (None, Some(_)) => return node2.hit(ray, t_min, t_max),
            (Some(t1), Some(t2)) => (t1, t2),
        };

        // Both bounding boxes were hit.
        // Now check for hits inside each node, starting with the closer of the two.
        if t1 > t2 {
            std::mem::swap(&mut node1, &mut node2);
        }
        if let Some(hit) = node1.hit(ray, t_min, t_max) {
            return Some(hit);
        } else {
            return node2.hit(ray, t_min, t_max);
        }
    }
}

fn center_partition(
    objects: Vec<Box<dyn BoundedSurface>>,
) -> (Vec<Box<dyn BoundedSurface>>, Vec<Box<dyn BoundedSurface>>) {
    let mut upper = Vec3D::fill(f32::NEG_INFINITY);
    let mut lower = Vec3D::fill(f32::INFINITY);
    for obj in objects[1..].iter() {
        let centroid = obj.bounding_box().centroid();
        upper = Vec3D::max(upper, centroid);
        lower = Vec3D::min(lower, centroid);
    }

    let widths = upper - lower;
    let mut max_width = 0.0;
    let mut axis = 0;
    for i in 0..3 {
        if widths[i] > max_width {
            max_width = widths[i];
            axis = i;
        }
    }

    let center = (upper + lower)[axis] / 2.0;
    return objects
        .into_iter()
        .partition(|obj| obj.bounding_box().centroid()[axis] < center);
}
