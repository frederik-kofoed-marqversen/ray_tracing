extern crate fastrand;
use super::bsdf::{to_local_coords, to_world_coords};
use super::{Camera, Light, Object, Ray, Vec3D};
use fastrand::Rng;
use std::io;
use std::io::Write;
use std::rc::Rc;

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
    lights: Vec<Rc<dyn Light>>,
    camera: Camera,
}

impl Engine {
    pub fn new(objects: Vec<Object>, lights: Vec<Rc<dyn Light>>, camera: Camera) -> Self {
        return Self {
            objects,
            lights,
            camera,
        };
    }

    fn hit(&self, ray: &Ray, t_min: f32, t_max: f32) -> Option<(f32, Vec3D, &Object)> {
        let mut result: Option<(f32, Vec3D, &Object)> = None;
        let mut distance = f32::INFINITY;
        for object in self.objects.iter() {
            if let Some((t, normal)) = object.hit(ray, t_min, t_max) {
                if t < distance {
                    // the hit object is closer
                    distance = t;
                    result = Some((t, normal, object));
                }
            }
        }
        return result;
    }

    fn hit_light(&self, ray: &Ray, t_min: f32, t_max: f32) -> Option<(f32, Vec3D, Vec3D)> {
        let mut result = None;
        let mut distance = f32::INFINITY;
        for light in self.lights.iter() {
            if let Some((t, normal, radiance)) = light.hit(ray, t_min, t_max) {
                if t < distance {
                    // the hit object is closer
                    distance = t;
                    result = Some((t, normal, radiance));
                }
            }
        }
        return result;
    }

    fn sample_lights_uniform(&self, rng: &mut Rng) -> (&Rc<dyn Light>, f32) {
        let n = self.lights.len();
        let light_index = rng.usize(0..n);
        let light = &self.lights[light_index];
        let p = 1.0 / n as f32;

        (light, p)
    }

    fn ray_colour(&self, mut ray: Ray, mut depth: u32, rng: &mut Rng) -> Vec3D {
        const T_MIN: f32 = 0.001;
        const T_MAX: f32 = f32::INFINITY;
        const THRESHOLD: f32 = 0.001;

        let mut beta = Vec3D::ONES;
        let mut result = Vec3D::ZERO;

        while beta.x + beta.y + beta.z > THRESHOLD && depth > 0 {
            depth -= 1;

            // intersect ray with scene
            let hit = self.hit(&ray, T_MIN, T_MAX);
            let t_max = hit.map_or(T_MAX, |(t, _, _)| t);
            // check if ray hits a light source before object
            if let Some((_, _, radiance)) = self.hit_light(&ray, T_MIN, t_max) {
                result += Vec3D::mul_elemwise(beta, radiance);
                break; // lights are assumed to not reflect light
            }
            if hit.is_none() {
                return result;
            }
            let (t, normal, object) = hit.unwrap();
            let intersect = ray.at(t);
            let dir_out = to_local_coords(-ray.direction, normal);

            // direct light sampling
            let (light, p) = self.sample_lights_uniform(rng);
            let light_sample = light.sample(intersect, rng);
            let shadow_ray = Ray {
                origin: intersect,
                direction: light_sample.direction,
            };
            let dir_in = to_local_coords(shadow_ray.direction, normal);
            let bsdf = &object.material.eval(dir_in, dir_out) * dir_in.z.abs();
            if !bsdf.is_zero() && self.hit(&shadow_ray, T_MIN, T_MAX).is_none() {
                // light is visible and samplable by BSDF
                result +=
                    Vec3D::mul_elemwise(Vec3D::mul_elemwise(beta, light_sample.spectrum), bsdf)
                        / (light_sample.pdf * p);
            };

            // indirect light sampling
            let bsdf_sample = object.material.sample(dir_out, rng);
            beta = Vec3D::mul_elemwise(beta, bsdf_sample.spectrum) * bsdf_sample.dir_in.z.abs()
                / bsdf_sample.pdf;

            // update ray for next iteration
            ray = Ray {
                origin: intersect,
                direction: to_world_coords(bsdf_sample.dir_in, normal),
            };
        }

        return result;
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
            if j % 10 == 0 {
                writeln!(err_lock, "Scanlines remaining: {j}")?;
            }
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
