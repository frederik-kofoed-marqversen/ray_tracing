use super::{Vec3D, Ray, Object, Colour, Rng};

pub struct Environment;
impl Environment {
    pub fn colour(&self, ray: &Ray) -> Colour {
        return Colour::zero();
        /* let t = 0.5*(ray.direction.z/ray.direction.norm() + 1.0);
        let colour = (1.0-t)*Colour::ones() + t*Colour::new(0.5, 0.7, 1.0);
        return colour; */
    }
}

pub struct Scene {
    pub objects: Vec<Object>,
    pub environment: Environment,
}

impl Scene {
    pub fn new() -> Self {
        Self{objects: Vec::new(), environment: Environment}
    }

    pub fn add(&mut self, object: Object) {
        self.objects.push(object);
    }

    pub fn hit(&self, ray: &Ray, t_min: f64, t_max: f64) -> Option<(f64, Vec3D, &Object)> {
        let mut result: Option<(f64, Vec3D, &Object)> = None;
        let mut distance = f64::INFINITY;
        for object in self.objects.iter() {
            match object.surface.hit(ray, t_min, t_max) {
                Some((t, normal)) => {  // the ray has hit something
                    if t < distance {  // the hit object is closer
                        distance = t;
                        result = Some((t, normal, object));
                    }
                },
                None => {}
            }
        }
        return result
    }

    pub fn ray_colour(&self, mut ray: Ray, depth: u32, rng: &mut impl Rng) -> Colour {
        if depth == 0 {return Colour::zero()}
        match self.hit(&ray, 0.001, f64::INFINITY) {
            Some((t, normal, obj)) => {  // an object has been hit
                ray.origin = ray.at(t);
                let emitted_light = obj.material.emission * obj.material.intensity;
                let mut scattered_direction = &normal + Vec3D::random(1.0, rng);
                if scattered_direction.almost_zero(1e-12) {
                    scattered_direction = normal.clone();
                }

                if rng.gen::<f64>() < obj.material.reflectiveness { // specular reflection
                    let reflected_direction = ray.direction - 2.0 * Vec3D::dot(&ray.direction, &normal) * normal;
                    ray.direction = Vec3D::interpolate(
                        &reflected_direction,
                        &scattered_direction, 
                        obj.material.roughness
                    );
                    return emitted_light + obj.material.specular_colour * self.ray_colour(ray, depth-1, rng);
                } else {  // diffuse scattering
                    ray.direction = scattered_direction;
                    return emitted_light + obj.material.albedo * self.ray_colour(ray, depth-1, rng);
                }
            },
            None => {
                return Colour::zero();
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::{Sphere, Point3D, Material};

    #[test]
    fn hit() {
        let mut world = Scene::new();
        world.add(Object{
            surface: Box::new(Sphere::new(Point3D::new(1.0, 1.0, 1.0), 1.0)),
            material: Material::default()
        });
        world.add(Object{
            surface: Box::new(Sphere::new(Point3D::e3(), 0.5)),
            material: Material::default()
        });
        let ray = Ray::new(
            Point3D::zero(),
            Vec3D::e3()
        );
        let (t, _, _) = world.hit(&ray, 0.0, f64::INFINITY).unwrap();
        assert_eq!(t, 0.5);
    }
}