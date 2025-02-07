extern crate fastrand;
use super::materials::Colour;
use super::{Camera, Object, Ray, Vec3D};
use fastrand::Rng;
use std::io;
use std::io::Write;

#[inline]
pub fn clamp(val: f32, min: f32, max: f32) -> f32 {
    if val < min {
        return min;
    }
    if val > max {
        return max;
    }
    return val;
}

fn write_pixel(
    lock: &mut io::StdoutLock,
    pixel_colour: Colour,
) -> io::Result<()> {
    writeln!(
        lock,
        "{} {} {}",
        (256.0 * clamp(pixel_colour.x.sqrt(), 0.0, 0.999)) as u8,
        (256.0 * clamp(pixel_colour.y.sqrt(), 0.0, 0.999)) as u8,
        (256.0 * clamp(pixel_colour.z.sqrt(), 0.0, 0.999)) as u8,
    )
}

fn sky_colour(ray: &Ray) -> Colour {
    let t = 0.5 * (ray.direction.z / ray.direction.norm() + 1.0);
    let colour = (1.0 - t) * Colour::ones() + t * Colour::new(0.5, 0.7, 1.0);
    return colour;
}

fn reflection_coefficient(eta: f32, cos_theta: f32) -> f32 {
    let sin_theta = f32::sqrt(1.0 - cos_theta * cos_theta);
    if eta * sin_theta > 1.0 {
        return 1.0;
    }
    // Schlick's approximation
    let r0 = f32::powi((eta - 1.0) / (eta + 1.0), 2);
    return r0 + (1.0 - r0) * f32::powi(1.0 - cos_theta, 5);
}

fn refract(ray_direction: &Vec3D, normal: &Vec3D, eta: f32, cos_theta: f32) -> Vec3D {
    let refracted_parallel = eta * (ray_direction + cos_theta * normal);
    let refracted_orthogonal = -f32::sqrt(1.0 - refracted_parallel.norm_squared()) * normal;
    return refracted_parallel + refracted_orthogonal;
}

fn reflect(ray_direction: &Vec3D, normal: &Vec3D) -> Vec3D {
    ray_direction - 2.0 * Vec3D::dot(&ray_direction, &normal) * normal
}

fn scatter(ray_direction: &Vec3D, normal: &Vec3D, roughness: f32, rng: &mut Rng) -> Vec3D {
    // Compute scattering direction using true Lambertian distribution.
    // This is done by adding a random unit vector from a spherically symmetric distribution
    // to the normal vector.
    // roughness=0 => reflection, roughness=1 => diffuse

    let mut scattered_direction = normal + Vec3D::random_rejection(rng);
    // let mut scattered_direction = &normal + Vec3D::random_direct(rng);  // alternative

    // Remove almost zero cases to avoid rounding errors.
    // Direction zero is the same as the normal direction.
    if scattered_direction.almost_zero(1e-12) {
        scattered_direction = normal.clone();
    }

    // Add scattering direction to ray direction according to roughness
    // roughness=0 => reflection, roughness=1 => diffuse
    let result = Vec3D::interpolate(ray_direction, &scattered_direction, roughness);

    return result / result.norm();
}

pub struct Engine {
    objects: Vec<Object>,
    camera: Camera,
}

impl Engine {
    pub fn new(objects: Vec<Object>, camera: Camera) -> Self {
        return Self { objects, camera };
    }

    fn hit(&self, ray: &Ray, t_min: f32, t_max: f32) -> Option<(f32, Vec3D, &Object)> {
        let mut result: Option<(f32, Vec3D, &Object)> = None;
        let mut distance = f32::INFINITY;
        for object in self.objects.iter() {
            match object.surface.hit(ray, t_min, t_max) {
                Some((t, normal)) => {
                    // the ray has hit something
                    if t < distance {
                        // the hit object is closer
                        distance = t;
                        result = Some((t, normal, object));
                    }
                }
                None => {}
            }
        }
        return result;
    }

    fn ray_colour(&self, ray: Ray, depth: u32, rng: &mut Rng) -> Colour {
        if depth == 0 {
            return Colour::zero();
        }
        match self.hit(&ray, 0.001, f32::INFINITY) {
            Some((t, mut normal, obj)) => {
                // an object has been hit
                // Compute new emitted ray
                let intersection = ray.at(t);

                // Compute angle of incidence and normal vector pointing up from surface in direction of where
                // the ray came from. Also compute eta = n1/n2 where n1 is the refractive index of material
                // where the ray comes from and n2 that of the material that ray goes into.
                let mut cos_theta = -Vec3D::dot(&ray.direction, &normal); // cosine of angle of incidence
                let mut eta = obj.material.refractive_index;
                if cos_theta < 0.0 {
                    cos_theta = -cos_theta;
                    normal = -normal;
                } else {
                    eta = eta.and_then(|n| Some(1.0 / n));
                }

                // Compute reflection coefficient and compute emitted ray by reflection or refraction
                let reflection_coefficient = if let Some(eta) = eta {
                    reflection_coefficient(eta, cos_theta)
                } else {
                    1.0
                };
                let mut emitted_direction;
                if rng.f32() < reflection_coefficient {
                    // Reflect
                    emitted_direction = reflect(&ray.direction, &normal);
                } else {
                    // Refract
                    emitted_direction =
                        refract(&ray.direction, &normal, eta.unwrap_or(1.0), cos_theta);
                    // Refracted rays scatter along the negative of the normal vector
                    normal = -normal;
                };

                // Add scattering from roughness
                emitted_direction =
                    scatter(&emitted_direction, &normal, obj.material.roughness, rng);

                // Compute colour of emitted ray
                let emitted_ray = Ray {
                    origin: intersection,
                    direction: emitted_direction,
                };
                let mut colour = obj.material.albedo * self.ray_colour(emitted_ray, depth - 1, rng);

                // Add emitted light from the object
                if let Some((emission, intensity)) = obj.material.emission {
                    colour = colour + emission * intensity
                }

                return colour;
            }
            None => {
                // return sky_colour(&ray);
                return Colour::zero();
            }
        }
    }

    pub fn render(
        &self,
        aspect_ratio: f32,
        image_width: usize,
        samples_per_pixel: u32,
        ray_depth: u32,
    ) -> io::Result<()> {
        let image_height: usize = (image_width as f32 / aspect_ratio) as usize;
        let mut rng = fastrand::Rng::new();

        let stdout = io::stdout();
        let mut lock = stdout.lock();
        let stderr = io::stderr();
        let mut err_lock = stderr.lock();

        write!(lock, "P3\n{image_width} {image_height}\n255\n")?;
        for j in (0..image_height).rev() {
            writeln!(err_lock, "Scanlines remaining: {j}")?;
            for i in 0..image_width {
                let mut colour = Colour::zero();
                for _ in 0..samples_per_pixel {
                    let u = (i as f32 + rng.f32()) / (image_width - 1) as f32;
                    let v = (j as f32 + rng.f32()) / (image_height - 1) as f32;
                    let ray = self.camera.ray(u, v);
                    colour += self.ray_colour(ray, ray_depth, &mut rng);
                }
                write_pixel(&mut lock, colour / samples_per_pixel as f32)?;
            }
        }
        writeln!(err_lock, "Done.")?;

        Ok(())
    }
}
