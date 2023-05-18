extern crate fastrand;
use std::io;
use std::io::Write;
use fastrand::Rng;
use super::{Vec3D, Ray, Object, Camera};
use super::materials::Colour;

#[inline]
pub fn clamp(val: f32, min: f32, max: f32) -> f32 {
    if val < min {return min}
    if val > max {return max}
    return val
}

fn write_pixel(lock: &mut io::StdoutLock, mut pixel_colour: Colour, num_samples: u32) -> io::Result<()> {
    pixel_colour *= 1.0 / num_samples as f32;
    
    writeln!(lock, "{} {} {}", 
        (256.0 * clamp(pixel_colour.x.sqrt(), 0.0, 0.999)) as u8,
        (256.0 * clamp(pixel_colour.y.sqrt(), 0.0, 0.999)) as u8,
        (256.0 * clamp(pixel_colour.z.sqrt(), 0.0, 0.999)) as u8,
    )
}

fn sky_colour(ray: &Ray) -> Colour {
    let t = 0.5*(ray.direction.z/ray.direction.norm() + 1.0);
    let colour = (1.0-t)*Colour::ones() + t*Colour::new(0.5, 0.7, 1.0);
    return colour;
}

pub struct Engine {
    objects: Vec<Object>,
    camera: Camera,
}

impl Engine {
    pub fn new(objects: Vec<Object>, camera: Camera) -> Self {
        return Self {objects, camera}
    }

    fn hit(&self, ray: &Ray, t_min: f32, t_max: f32) -> Option<(f32, Vec3D, &Object)> {
        let mut result: Option<(f32, Vec3D, &Object)> = None;
        let mut distance = f32::INFINITY;
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

    fn ray_colour(&self, mut ray: Ray, depth: u32, rng: &mut Rng) -> Colour {
        if depth == 0 {return Colour::zero()}
        match self.hit(&ray, 0.001, f32::INFINITY) {
            Some((t, normal, obj)) => {  // an object has been hit
                ray.origin = ray.at(t);
                let emitted_light = obj.material.emission * obj.material.intensity;
                
                let reflected_direction = ray.direction - 2.0 * Vec3D::dot(&ray.direction, &normal) * normal;
                let mut scattered_direction = &normal + Vec3D::random_rejection(rng);
                // let mut scattered_direction = &normal + Vec3D::random_direct(rng);  // alternative
                if scattered_direction.almost_zero(1e-12) {
                    scattered_direction = normal.clone();
                }
                
                ray.direction = Vec3D::interpolate(
                    &reflected_direction,
                    &scattered_direction, 
                    obj.material.roughness
                ); // roughness=0 => reflection, roughness=1 => diffuse

                return emitted_light + obj.material.albedo * self.ray_colour(ray, depth-1, rng);
            },
            None => {
                // return sky_colour(&ray);
                return Colour::zero();
            }
        }
    }

    pub fn render(&self, aspect_ratio: f32, image_width: usize, samples_per_pixel: u32, ray_depth: u32) -> io::Result<()> {
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
                    let u = (i as f32 + rng.f32()) / (image_width-1) as f32;
                    let v = (j as f32 + rng.f32()) / (image_height-1) as f32;
                    let ray = self.camera.ray(u, v);
                    colour += self.ray_colour(ray, ray_depth, &mut rng);
                }
                write_pixel(&mut lock, colour, samples_per_pixel)?;
            }
        }
        writeln!(err_lock, "Done.")?;
    
        Ok(())
    }
}