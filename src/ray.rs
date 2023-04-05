use super::{Point3D, Vec3D};

#[derive(Debug, Copy, Clone)]
pub struct Ray {
    pub origin: Point3D,
    pub direction: Vec3D,
}

impl Ray {
    pub fn new(origin: Point3D, direction: Vec3D) -> Self {
        Self{origin, direction}
    }
    
    pub fn at(&self, t: f64) -> Point3D {
        self.origin + self.direction * t
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ray_at() {
        let point = Point3D::new(1.0, 1.0, 1.0);
        let vec = Vec3D::new(0.0, 1.0, 2.0);
        let ray = Ray::new(point, vec);
        let point = ray.at(2.0);
        assert_eq!(point.x, 1.0);
        assert_eq!(point.y, 3.0);
        assert_eq!(point.z, 5.0);
    }
}