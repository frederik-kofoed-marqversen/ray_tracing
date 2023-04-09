#[macro_use]
extern crate impl_ops;
extern crate rand;
extern crate rand_distr;
mod vec3d;
mod ray;
mod surfaces;
mod scene;
mod camera;

use std::io;
use std::io::Write;
use rand::Rng;
use rand_distr::StandardNormal;

use vec3d::Vec3D;
use ray::Ray;
use surfaces::*;
use scene::Scene;

type Colour = Vec3D;
type Point3D = Vec3D;

#[inline]
pub fn clamp(val: f64, min: f64, max: f64) -> f64 {
    if val < min {return min}
    if val > max {return max}
    return val
}

fn write_pixel(lock: &mut io::StdoutLock, mut pixel_colour: Colour, num_samples: u32) -> io::Result<()> {
    pixel_colour *= 1.0 / num_samples as f64;
    
    writeln!(lock, "{} {} {}", 
        (256.0 * clamp(pixel_colour.x.sqrt(), 0.0, 0.999)) as u8,
        (256.0 * clamp(pixel_colour.y.sqrt(), 0.0, 0.999)) as u8,
        (256.0 * clamp(pixel_colour.z.sqrt(), 0.0, 0.999)) as u8,
    )
}

fn main() -> io::Result<()> {
    const ASPECT_RATIO: f64 = 16.0 / 9.0;
    const IMAGE_WIDTH: usize = 800;
    const IMAGE_HEIGHT: usize = (IMAGE_WIDTH as f64 / ASPECT_RATIO) as usize;
    const SAMPLES_PER_PIXEL: u32 = 512;
    const RAY_DEPTH: u32 = 8;

    // World setup
    let mut world = Scene::new();
    
    // Balls
    world.add(Object{
        surface: Box::new(Sphere::new(Point3D::new(-2.0, 0.0, 0.5), 0.5)),
        material: Material::diffuse(Colour::e1())
    });
    world.add(Object{
        surface: Box::new(Sphere::new(Point3D::new(-0.5, 0.0, 0.8), 0.75)),
        // material: Material::diffuse(Colour::e2())
        material: Material::metal(Colour::new(1.0, 0.5, 0.5), 0.5) 
    });
    world.add(Object{
        surface: Box::new(Sphere::new(Point3D::new(1.5, 0.0, 1.0), 1.0)),
        material: Material::diffuse(Colour::e3())
    });
    world.add(Object{
        surface: Box::new(Sphere::new(Point3D::new(2.5, -2.5, 0.5), 0.5)),
        material: Material { 
            albedo: Colour::e2(),
            specular_colour: Colour::zero(),
            reflectiveness: 0.0,
            roughness: 0.0,
            emission: Colour::e2(), 
            intensity: 1.0
        }
    });
    world.add(Object{
        surface: Box::new(Triangle::new(
            Point3D::new(-4.5, 2.0, 0.2),
            Point3D::new(-2.0, 3.0, 0.2),
            Point3D::new(-3.0, 2.5, 3.2),
        )),
        material: Material::mirror()
    });

    // Globe
    world.add(Object{
        surface: Box::new(Sphere::new(Point3D::new(0.0, 0.0, -100.0), 100.0)),
        material: Material::diffuse(Colour::ones()*0.5)
    });

    // Light
    world.add(Object{
        surface: Box::new(Sphere::new(Point3D::new(0.0, 110.0, 0.0), 100.0)),
        material: Material::light_source(Colour::ones(), 4.0)
    });

    // Camera
    let cam = camera::Camera::new(
        Point3D::new(-2.5, -4.0, 3.0),
        Vec3D::new(0.5, 1.0, -0.5),
        ASPECT_RATIO
    );

    // Render image
    let mut rng = rand::thread_rng();

    let stdout = io::stdout();
    let mut lock = stdout.lock();
    let stderr = io::stderr();
    let mut err_lock = stderr.lock();

    write!(lock, "P3\n{IMAGE_WIDTH} {IMAGE_HEIGHT}\n255\n")?;
    for j in (0..IMAGE_HEIGHT).rev() {
        writeln!(err_lock, "Scanlines remaining: {j}")?;
        for i in 0..IMAGE_WIDTH {
            let mut colour = Colour::zero();
            for _ in 0..SAMPLES_PER_PIXEL {
                let u = (i as f64 + rng.gen::<f64>()) / (IMAGE_WIDTH-1) as f64;
                let v = (j as f64 + rng.gen::<f64>()) / (IMAGE_HEIGHT-1) as f64;
                let ray = cam.ray(u, v);
                colour += world.ray_colour(ray, RAY_DEPTH, &mut rng);
            }
            write_pixel(&mut lock, colour, SAMPLES_PER_PIXEL)?;
        }
    }
    writeln!(err_lock, "Done.")?;

    Ok(())
}
