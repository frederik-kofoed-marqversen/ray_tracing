extern crate fastrand;
use super::bsdf::{to_local_coords, to_world_coords};
use super::{Camera, Light, Object, Ray, Vec3D};
use fastrand::Rng;
use std::io;
use std::io::Write;
use std::rc::Rc;

#[inline]
fn write_pixel(lock: &mut io::StdoutLock, pixel_colour: Vec3D) -> io::Result<()> {
    writeln!(
        lock,
        "{} {} {}",
        (256.0 * pixel_colour.x.sqrt().clamp(0.0, 0.999)) as u8,
        (256.0 * pixel_colour.y.sqrt().clamp(0.0, 0.999)) as u8,
        (256.0 * pixel_colour.z.sqrt().clamp(0.0, 0.999)) as u8,
    )
}

#[inline]
fn mis_weight(power: i32, pa: f32, pb: f32) -> f32 {
    // power heuristic
    let a = pa.powi(power);
    let b = pb.powi(power);
    if a + b == 0.0 {
        0.0
    } else {
        a / (a + b)
    }
}

pub struct Engine {
    objects: Vec<Object>,
    lights: Vec<Rc<dyn Light>>,
    camera: Camera,
}

impl Engine {
    #[inline]
    pub fn new(objects: Vec<Object>, lights: Vec<Rc<dyn Light>>, camera: Camera) -> Self {
        return Self {
            objects,
            lights,
            camera,
        };
    }

    #[inline]
    fn hit(&self, ray: &Ray, t_min: f32, mut t_max: f32) -> Option<(f32, Vec3D, &Object)> {
        let mut result: Option<(f32, Vec3D, &Object)> = None;
        for object in self.objects.iter() {
            if let Some((t, normal)) = object.hit(ray, t_min, t_max) {
                t_max = t;
                result = Some((t, normal, object));
            }
        }
        return result;
    }

    #[inline]
    fn hit_bool(&self, ray: &Ray, t_min: f32, t_max: f32) -> bool {
        self.objects
            .iter()
            .any(|obj| obj.hit_bool(ray, t_min, t_max))
    }

    #[inline]
    fn hit_light(
        &self,
        ray: &Ray,
        t_min: f32,
        mut t_max: f32,
    ) -> Option<(f32, Vec3D, &Rc<dyn Light>)> {
        let mut result = None;
        for light in self.lights.iter() {
            if let Some((t, _normal, radiance)) = light.hit(ray, t_min, t_max) {
                t_max = t;
                result = Some((t, radiance, light));
            }
        }
        return result;
    }

    #[inline]
    fn sample_lights_uniform(&self, rng: &mut Rng) -> (&Rc<dyn Light>, f32) {
        let n = self.lights.len();
        let light_index = rng.usize(0..n);
        let light = &self.lights[light_index];
        let p = 1.0 / n as f32;

        (light, p)
    }

    #[inline]
    fn lights_pdf(&self, _light: &Rc<dyn Light>) -> f32 {
        1.0 / self.lights.len() as f32
    }

    pub fn ray_colour(&self, mut ray: Ray, mut depth: u32, rng: &mut fastrand::Rng) -> Vec3D {
        const T_MIN: f32 = 0.001;
        const T_MAX: f32 = f32::INFINITY;

        let mut radiance = Vec3D::ZERO; // accumulated radiance
        let mut beta = Vec3D::ONES; // path throughput
        // let mut eta_scaling = 1.0;  // scaling factor for refractive index changes along the path

        // `prev_pdf` stores the pdf of the sampling method that produced `ray`
        // (when the current ray was generated). For the camera primary ray we
        // initialize it to 1.0 so MIS with emission treats camera as delta-like.
        let mut prev_pdf: f32 = 1.0;
        let mut prev_specular = true;

        while depth > 0 {
            depth -= 1;

            // Intersect with scene geometry
            let hit = self.hit(&ray, T_MIN, T_MAX);
            let t_hit = hit.map_or(T_MAX, |(t, _, _)| t);

            // Intersect with lights and add contribution
            if let Some((_, emitted, light)) = self.hit_light(&ray, T_MIN, t_hit) {
                if prev_specular {
                    radiance += Vec3D::mul_elemwise(beta, emitted);
                } else {
                    let p_light = self.lights_pdf(light) * light.pdf(&ray);
                    let w = mis_weight(2, prev_pdf, p_light);
                    radiance += w * Vec3D::mul_elemwise(beta, emitted);
                }
                break; // light terminates path (assumed black body with no reflection)
            }

            // No light hit: if no object hit either, ray escapes scene
            let (t, normal, object) = match hit {
                Some(h) => h,
                None => break,
            };

            let interaction_point = ray.at(t);
            let wo_local = to_local_coords(-ray.direction, normal);
            prev_specular = object.material.is_specular();

            // If surface is non-specular, do explicit direct-light sampling
            if !prev_specular {
                let (light, select_light_pdf) = self.sample_lights_uniform(rng);
                let light_sample = light.sample(interaction_point, rng);

                // Ensure sampled light is visible
                let shadow_ray = Ray {
                    origin: interaction_point,
                    direction: light_sample.direction,
                };
                if !self.hit_bool(&shadow_ray, T_MIN, T_MAX) {
                    let wi_local = to_local_coords(shadow_ray.direction, normal);

                    // pdf of sampling this direction via light-sampling strategy
                    let p_light = light_sample.pdf * select_light_pdf;

                    // bsdf value * cos theta
                    let bsdf_val = object.material.eval(wi_local, wo_local);
                    let cos = wi_local.z.abs();
                    let contribution = Vec3D::mul_elemwise(light_sample.spectrum, bsdf_val * cos);

                    // pdf of sampling the same direction via BSDF sampling
                    let p_bsdf = object.material.pdf(wi_local, wo_local);

                    let w = mis_weight(2, p_light, p_bsdf);
                    if p_light > 0.0 {
                        radiance += w * Vec3D::mul_elemwise(beta, contribution) / p_light;
                    }
                }
            }

            // Sample BSDF to get new direction (indirect lighting)
            let bsdf_sample = object.material.sample(wo_local, rng);
            let bsdf_val_cos = bsdf_sample.spectrum * bsdf_sample.dir_in.z.abs();

            // Update throughput and prepare next ray
            beta = Vec3D::mul_elemwise(beta, bsdf_val_cos) / bsdf_sample.pdf;
            // eta_scaling /= bsdf_sample.eta_rel_it().powi(2); // only upon transmission
            prev_pdf = bsdf_sample.pdf;
            prev_specular = object.material.is_specular();

            ray = Ray {
                origin: interaction_point,
                direction: to_world_coords(bsdf_sample.dir_in, normal),
            };

            // Russian roulette for path termination
            let q = 1.0 - beta.max_elem(); // termination probability (removed redundant .max(0.0))
            // let eta_corrected_beta = beta * eta_scaling;
            // let q = 1.0 - eta_corrected_beta.max_elem();
            if rng.f32() < q {
                // terminate path
                break;
            }
            beta /= 1.0 - q;
        }

        return radiance;
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
        let progress_interval = image_height / 10;
        for j in (0..image_height).rev() {
            if j % progress_interval == 0 {
                let progress = 100 - (j * 100) / image_height;
                writeln!(err_lock, "Progress: {progress}%")?;
            }
            for i in 0..image_width {
                let mut colour = Vec3D::ZERO;
                for _ in 0..samples_per_pixel {
                    let u = (i as f32 + rng.f32()) / (image_width - 1) as f32;
                    let v = (j as f32 + rng.f32()) / (image_height - 1) as f32;
                    let ray = self.camera.ray(u, v);
                    colour += self.ray_colour(ray, ray_depth, &mut rng);
                }
                if colour.x.is_nan() || colour.y.is_nan() || colour.z.is_nan() {
                    colour = Vec3D::new(255.0, 255.0, 0.0);
                }
                write_pixel(&mut lock, colour / samples_per_pixel as f32)?;
            }
        }
        writeln!(err_lock, "Done.")?;

        Ok(())
    }
}
