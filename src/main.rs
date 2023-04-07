#[macro_use]
extern crate impl_ops;
mod vec3d;
mod ray;
mod surfaces;

use std::io;
use std::io::Write;

use vec3d::Vec3D;
use ray::Ray;
use surfaces::{Sphere, Environment};

type Colour = Vec3D;
type Point3D = Vec3D;

fn write_pixel(lock: &mut io::StdoutLock, pixel_colour: Colour) -> io::Result<()> {
    writeln!(lock, "{} {} {}", 
        (pixel_colour.x * 255.999) as u8,
        (pixel_colour.y * 255.999) as u8,
        (pixel_colour.z * 255.999) as u8,
    )
}

fn ray_colour(ray: &Ray, world: &Environment) -> Colour {
    match world.hit(ray, 0.0, f64::INFINITY) {
        Some((t, surface)) => {
            let normal = surface.normal(ray.at(t));
            return 0.5 * Colour::new(normal.x+1.0, normal.z+1.0, -normal.y+1.0)
        },
        None => {
            let t = 0.5 * (ray.direction.z / ray.direction.norm() + 1.0);
            return (1.0 - t) * Colour::new(1.0, 1.0, 1.0) + t * Colour::new(0.5, 0.7, 1.0)
        }
    }
}

fn main() -> io::Result<()> {
    // Output image
    const ASPECT_RATIO: f64 = 16.0 / 9.0;
    const IMAGE_WIDTH: usize = 400;
    const IMAGE_HEIGHT: usize = (IMAGE_WIDTH as f64 / ASPECT_RATIO) as usize;

    // Camera
    let window_height: f64 = 2.0;
    let window_width: f64 = window_height * ASPECT_RATIO;
    let focal_length: f64 = 1.0;

    let origin = Point3D::zero();
    let horizontal = Vec3D::e1() * window_width;
    let vertical = Vec3D::e3() * window_height;
    let forward = Vec3D::e2() * focal_length;
    let lower_left_corner: Point3D = &origin + &forward - &horizontal / 2.0 - &vertical / 2.0;

    // Render image
    let stdout = io::stdout();
    let mut lock = stdout.lock();
    let stderr = io::stderr();
    let mut err_lock = stderr.lock();

    write!(lock, "P3\n{IMAGE_WIDTH} {IMAGE_HEIGHT}\n255\n")?;
    
    let mut world = Environment::new();
    world.add(Box::new(Sphere::new(Point3D::new(0.0, 1.0, 0.0), 0.5)));
    world.add(Box::new(Sphere::new(Point3D::new(0.0, 1.0, -100.5), 100.0)));

    for j in (0..IMAGE_HEIGHT).rev() {
        writeln!(err_lock, "Scanlines remaining: {j}")?;
        
        let v = j as f64 / (IMAGE_HEIGHT-1) as f64;
        for i in 0..IMAGE_WIDTH {
            let u = i as f64 / (IMAGE_WIDTH-1) as f64;
            let direction: Vec3D = &lower_left_corner + u*&horizontal + v*&vertical - &origin;
            let ray = Ray::new(origin.clone(), direction);
            let colour = ray_colour(&ray, &world);
            write_pixel(&mut lock, colour)?;
        }
    }
    writeln!(err_lock, "Done.")?;

    Ok(())
}
