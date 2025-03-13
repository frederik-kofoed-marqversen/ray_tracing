extern crate fastrand;
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

#[inline]
fn write_pixel(lock: &mut io::StdoutLock, pixel_colour: Vec3D) -> io::Result<()> {
    writeln!(
        lock,
        "{} {} {}",
        (256.0 * clamp(pixel_colour.x.sqrt(), 0.0, 0.999)) as u8,
        (256.0 * clamp(pixel_colour.y.sqrt(), 0.0, 0.999)) as u8,
        (256.0 * clamp(pixel_colour.z.sqrt(), 0.0, 0.999)) as u8,
    )
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
            match object.hit(ray, t_min, t_max) {
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

    fn ray_colour(&self, ray: Ray, depth: u32, rng: &mut Rng) -> Vec3D {
        if depth == 0 {
            return Vec3D::ZERO;
        }
        match self.hit(&ray, 0.001, f32::INFINITY) {
            Some((t, normal, obj)) => {
                // an object has been hit
                // Compute new emitted ray
                let intersection = ray.at(t);

                let emitted_dir = obj.material.sample_bsdf(-ray.direction, normal, rng);

                // Compute colour of emitted ray
                let emitted_ray = Ray::new(intersection, emitted_dir);
                let mut colour = Vec3D::mul_elemwise(
                    obj.material.albedo,
                    self.ray_colour(emitted_ray, depth - 1, rng),
                );

                // Add emitted light from the object
                if let Some((emission, intensity)) = obj.material.emission {
                    colour = colour + emission * intensity
                }

                return colour;
            }
            None => {
                // return sky_colour(&ray);
                return Vec3D::ZERO;
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
                let mut colour = Vec3D::ZERO;
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
