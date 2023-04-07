use super::{Vec3D, Point3D, Ray};

pub trait Surface {
    fn hit(&self, ray: &Ray, t_min: f64, t_max: f64) -> Option<f64>;

    fn normal(&self, point: Point3D) -> Vec3D;
}

pub struct Sphere {
    center: Point3D,
    radius: f64,
}

impl Sphere {
    pub fn new(center: Point3D, radius: f64) -> Self {
        Self{center, radius}
    }
}

impl Surface for Sphere {
    fn hit(&self, ray: &Ray, t_min: f64, t_max: f64) -> Option<f64> {
        let r_oc: Vec3D = ray.origin - self.center;
        let a = ray.direction.norm_squared();
        let half_b = Vec3D::dot(&ray.direction, &r_oc);
        let c = r_oc.norm_squared() - self.radius * self.radius;
        let d = half_b*half_b - a*c;

        if d < 0.0 {return None}
        let sqrt_d = d.sqrt();
        
        let mut t = -(half_b + sqrt_d) / a;
        if t < t_max && t > t_min {return Some(t)}
        
        t = (-half_b + sqrt_d) / a;
        if t < t_max && t > t_min {return Some(t)}

        return None
    }

    #[inline]
    fn normal(&self, point: Point3D) -> Vec3D {
        return (point - self.center) / self.radius
    }
}

pub struct Environment {
    pub surfaces: Vec<Box<dyn Surface>>,
}

impl Environment {
    pub fn new() -> Self {
        Self{surfaces: Vec::new()}
    }

    pub fn add(&mut self, surface: Box<dyn Surface>) {
        self.surfaces.push(surface);
    }

    pub fn hit(&self, ray: &Ray, t_min: f64, t_max: f64) -> Option<(f64, &Box<dyn Surface>)> {
        let mut hit = false;
        let mut closest = &self.surfaces[0];
        let mut distance = f64::INFINITY;

        for surface in self.surfaces.iter() {
            match surface.hit(ray, t_min, t_max) {
                Some(t) => {
                    hit = true;
                    if t < distance {
                        closest = surface;
                        distance = t;
                    }
                },
                None => {}
            }
        }

        if hit {
            return Some((distance, closest))
        } else {
            return None
        }
    }
}