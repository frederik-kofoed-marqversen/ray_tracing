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
        let n = objects.len();
        if n == 0 {
            panic!("zero objects given!");
        }
        if n == 1 {
            return objects.pop().unwrap()
        }
        if n == 2 {
            let bounding_box = AxisAlignedBoundingBox::combine(&objects[0].bounding_box(), &objects[1].bounding_box());
            let result: Box<dyn BoundedSurface> = Box::new(Self {node1:objects.pop().unwrap(), node2:objects.pop().unwrap(), bounding_box});
            return result
        }

        let (left, right) = center_partition(objects);

        let node1 = Self::build(left);
        let node2 = Self::build(right);
        let bounding_box = AxisAlignedBoundingBox::combine(&node1.bounding_box(), &node2.bounding_box());
        let result: Box<dyn BoundedSurface> = Box::new(Self {node1, node2, bounding_box});
        return result
    }
}

impl BoundedSurface for BoundingVolumeHierarchy {}

impl Bounded for BoundingVolumeHierarchy {
    fn bounding_box(&self) -> AxisAlignedBoundingBox {
        self.bounding_box.clone()
    }
}

impl Surface for BoundingVolumeHierarchy {
    fn hit(&self, ray: &Ray, t_min: f32, mut t_max: f32) -> Option<(f32, Vec3D)> {
        if self.bounding_box.hit(ray, t_min, t_max).is_none() {
            // ray will miss both left and right => no hit
            return None
        }

        let hit1 = self.node1.hit(ray, t_min, t_max);        
        if hit1.is_some() {
            // only hit right if closer than hit left
            t_max = hit1.unwrap().0;
        }
        match self.node2.hit(ray, t_min, t_max) {
            None => return hit1,
            hit2 => return hit2
        }
    }
}


fn center_partition(objects: Vec<Box<dyn BoundedSurface>>) -> (Vec<Box<dyn BoundedSurface>>, Vec<Box<dyn BoundedSurface>>) {
    let mut upper = objects[0].bounding_box().centroid();
    let mut lower = upper.clone();
    for obj in objects[1..].iter() {
        let centroid = obj.bounding_box().centroid();
        for i in 0..3 {
            upper[i] = upper[i].max(centroid[i]);
            lower[i] = lower[i].min(centroid[i]);
        }
    }
    let widths = &upper - &lower;
    
    let mut max_width = 0.0;
    let mut axis = 0;
    for i in 0..3 {
        if widths[i] > max_width {
            max_width = widths[i];
            axis = i;
        }
    }
    
    let center = (upper + lower)[axis] / 2.0;
    return objects.into_iter().partition(|obj| obj.bounding_box().centroid()[axis] < center)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn build_bvh() {
        BVH::build(Vec::new());
    }
}