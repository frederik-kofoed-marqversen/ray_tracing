use super::{Vec3D, Point3D, Colour, Ray};

#[derive(Default)]
pub struct Material {
    pub albedo: Colour,
    pub specular_colour: Colour,
    pub reflectiveness: f64,
    pub roughness: f64,
    pub emission: Colour,
    pub intensity: f64,
}

impl Material {
    #[inline]
    pub fn diffuse(albedo: Colour) -> Self {
        Self {
            albedo,
            specular_colour: Colour::zero(),
            reflectiveness: 0.0,
            roughness: 0.0,
            emission: Colour::zero(),
            intensity: 0.0,
        }
    }

    #[inline]
    pub fn light_source(emission: Colour, intensity: f64) -> Self {
        Self {
            albedo: Colour::zero(),
            specular_colour: Colour::zero(),
            reflectiveness: 0.0,
            roughness: 0.0,
            emission,
            intensity,
        }
    }

    #[inline]
    pub fn mirror() -> Self {
        Self {
            albedo: Colour::zero(),
            specular_colour: Colour::ones(),
            reflectiveness: 1.0,
            roughness: 0.0,
            emission: Colour::zero(),
            intensity: 0.0,
        }
    }

    #[inline]
    pub fn metal(specular_colour: Colour, roughness: f64) -> Self {
        Self {
            albedo: Colour::zero(),
            specular_colour,
            reflectiveness: 1.0,
            roughness,
            emission: Colour::zero(),
            intensity: 0.0,
        }
    }
}

pub trait Surface {
    fn hit(&self, ray: &Ray, t_min: f64, t_max: f64) -> Option<(f64, Vec3D)>;
}

pub struct Sphere {
    center: Point3D,
    radius: f64,
}

impl Sphere {
    pub fn new(center: Point3D, radius: f64) -> Self {
        Self{center, radius}
    }

    #[inline]
    fn normal(&self, point: Point3D) -> Vec3D {
        return (point - self.center) / self.radius
    }
}

impl Surface for Sphere {
    fn hit(&self, ray: &Ray, t_min: f64, t_max: f64) -> Option<(f64, Vec3D)> {
        let r_oc: Vec3D = ray.origin - self.center;
        let a = ray.direction.norm_squared();
        let half_b = Vec3D::dot(&ray.direction, &r_oc);
        let c = r_oc.norm_squared() - self.radius * self.radius;
        let d = half_b*half_b - a*c;

        if d < 0.0 {return None}
        let sqrt_d = d.sqrt();
        
        let mut t = -(half_b + sqrt_d) / a;
        if t < t_max && t > t_min {
            return Some((t, self.normal(ray.at(t))))
        }
        
        t = (-half_b + sqrt_d) / a;
        if t < t_max && t > t_min {
            return Some((t, self.normal(ray.at(t))))
        }

        return None
    }
}

pub struct Plane {
    d: f64,  // distance from origo to plane along the normal direction
    normal: Vec3D,
}

impl Plane {
    pub fn new(p: Point3D, mut normal: Vec3D) -> Self {
        normal.normalise();
        let d = -Vec3D::dot(&p, &normal);
        Self {d, normal}
    }
}

impl Surface for Plane {
    fn hit(&self, ray: &Ray, t_min: f64, t_max: f64) -> Option<(f64, Vec3D)> {
        let denom = Vec3D::dot(&self.normal, &ray.direction);
        if denom == 0.0 {return None}
        
        let t = - (Vec3D::dot(&self.normal, &ray.origin) + self.d) / denom;
        if t < t_max && t > t_min {
            return Some((t, self.normal))
        } else {
            return None
        }
    }
}

pub struct Triangle {
    p1: Point3D,
    p2: Point3D,
    p3: Point3D,
    d: f64,  // distance from origo to plane along the normal direction
    normal: Vec3D,
}

impl Triangle {
    pub fn new(p1: Point3D, p2: Point3D, p3: Point3D) -> Self {
        let v12 = &p2 - &p1;
        let v23 = &p3 - &p2;
        let normal = *Vec3D::cross(&v12, &v23).normalise();
        let d = -Vec3D::dot(&p1, &normal);
        Self {p1, p2, p3, d, normal}
    }
}

impl Surface for Triangle {
    fn hit(&self, ray: &Ray, t_min: f64, t_max: f64) -> Option<(f64, Vec3D)> {
        let denom = Vec3D::dot(&self.normal, &ray.direction);
        if denom == 0.0 {return None}
        
        let t = - (Vec3D::dot(&self.normal, &ray.origin) + self.d) / denom;
        if t > t_max || t < t_min {
            return None
        }

        let e1 = self.p2 - self.p1;
        let e2 = self.p3 - self.p2;
        let e3 = self.p1 - self.p3;

        let p = ray.at(t);
        let c1 = p - self.p1;
        let c2 = p - self.p2;
        let c3 = p - self.p3;

        if  Vec3D::dot(&self.normal, &Vec3D::cross(&e1, &c1)) < 0.0 ||
            Vec3D::dot(&self.normal, &Vec3D::cross(&e2, &c2)) < 0.0 ||
            Vec3D::dot(&self.normal, &Vec3D::cross(&e3, &c3)) < 0.0
        {
            return None
        } else {
            return Some((t, self.normal))
        }
    }
}

pub struct Object {
    pub surface: Box<dyn Surface>,
    pub material: Material
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hit_sphere() {
        let sphere = Sphere::new(Point3D::new(1.0, 1.0, 1.0), 1.0,);
        let ray = Ray::new(
            Point3D::zero(),
            Vec3D::new(1.0, 1.0, 1.0) / f64::sqrt(3.0)
        );
        let (t, normal) = sphere.hit(&ray, 0.0, f64::INFINITY).unwrap();
        assert!((f64::sqrt(3.0) - 1.0 - t).abs() < 1e-15);
        assert!((normal.norm() - 1.0).abs() < 1e-15);
    }

    #[test]
    fn hit_sphere_tangent() {
        let sphere = Sphere::new(Point3D::new(1.0, 1.0, 1.0), 1.0);
        let ray = Ray::new(
            Point3D::zero(),
            Vec3D::new(1.0, 1.0, 0.0) / f64::sqrt(2.0)
        );
        let (t, normal) = sphere.hit(&ray, 0.0, f64::INFINITY).unwrap();
        assert!((f64::sqrt(2.0) - t).abs() < 1e-15);
        assert!((normal.norm() - 1.0).abs() < 1e-15);
    }

    #[test]
    fn miss_sphere() {
        let sphere = Sphere::new(Point3D::new(1.0, 1.0, 1.0), 1.0);
        let ray = Ray::new(
            Point3D::zero(),
            Vec3D::e3()
        );
        let result = sphere.hit(&ray, 0.0, f64::INFINITY);
        assert_eq!(result, None);
    }
}