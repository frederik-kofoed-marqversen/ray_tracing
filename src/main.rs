#![feature(trait_upcasting)]

extern crate ray_tracer;
use ray_tracer::materials::*;
use ray_tracer::primitives::*;
use ray_tracer::subdivision_surface::loop_subdivide;
use ray_tracer::traits::BoundedSurface;
use ray_tracer::triangle_mesh::*;
use ray_tracer::Camera;
use ray_tracer::{Object, Point3D, Vec3D};

const ASPECT_RATIO: f32 = 16.0 / 9.0;
const IMAGE_WIDTH: usize = 400;
const SAMPLES_PER_PIXEL: u32 = 100;
const RAY_DEPTH: u32 = 4;

fn main() -> std::io::Result<()> {
    let (scene, camera) = scene_tetrahedron();
    // Render image
    let eng = ray_tracer::engine::Engine::new(scene, camera);
    return eng.render(ASPECT_RATIO, IMAGE_WIDTH, SAMPLES_PER_PIXEL, RAY_DEPTH);
}

#[allow(dead_code)]
fn scene_tetrahedron() -> (Vec<Object>, Camera) {
    // Regular tetrahedron centred at (0, 0, 0)
    let mesh = TriangleMesh {
        vertices: vec![
            Vec3D::new(1.0, 1.0, 1.0),
            Vec3D::new(1.0, -1.0, -1.0),
            Vec3D::new(-1.0, 1.0, -1.0),
            Vec3D::new(-1.0, -1.0, 1.0),
        ],
        triangles: vec![[0, 1, 2], [1, 3, 2], [0, 3, 1], [0, 2, 3]],
    };

    let mesh = loop_subdivide(&mesh, 1);
    let (_mesh, triangles) = mesh.to_triangles();
    let boxed_triangles: Vec<Box<dyn BoundedSurface>> = triangles
        .into_iter()
        .map(|obj| Box::new(obj) as Box<dyn BoundedSurface>)
        .collect();
    // let bvh = ray_tracer::accelerators::BVH::build(boxed_triangles);
    // let tetrahedron = Object {
    //     surface: bvh,
    //     material: Material::diffuse(Colour::new(1.0, 0.5, 0.5)),
    // };
    let mut triangles = boxed_triangles.into_iter().map(|obj| Object {
        surface: obj,
        material: Material::diffuse(Colour::new(1.0, 0.5, 0.5)),
    }).collect();

    // Ground to cast shadows onto
    let ground = Object {
        surface: Box::new(Plane::new(-1.5 * Vec3D::Z, Vec3D::Z)),
        material: Material::diffuse(Colour::ONES * 0.8),
    };

    // Sun for lighting
    let light = Object {
        surface: Box::new(Sphere::new(2000.0 * Point3D::new(1.0, 0.0, 1.0), 1000.0)),
        material: Material::light_source(Colour::ONES, 6.0),
    };

    // Camera
    let angle = 0.0 * std::f32::consts::PI;
    let loc = Vec3D::new_polar(3.0, 0.4 * std::f32::consts::PI, angle);
    let mut camera = ray_tracer::Camera::new(loc, Vec3D::ONES, ASPECT_RATIO);
    camera.look_at(Point3D::ZERO);

    let mut scene = vec![ground, light];
    scene.append(&mut triangles);
    return (scene, camera);
}

#[allow(dead_code)]
fn scene_yoda() -> (Vec<Object>, Camera) {
    // Baby Yoda (Grogu)
    let mesh = ray_tracer::stl::read_stl("./meshes/Baby_Yoda.stl").unwrap();
    let centre = mesh.vertices.iter().sum::<Vec3D>() / mesh.vertices.len() as f32;
    let ground_level = mesh
        .vertices
        .iter()
        .fold(f32::INFINITY, |min, p| min.min(p.z));

    // Unfortunately the .stl file is not good and does not define a manifold
    // mesh.try_fix_mesh();
    // let mesh = loop_subdivide(&mesh, 1);

    let (_, triangles) = mesh.to_triangles();
    let boxed_triangles = triangles
        .into_iter()
        .map(|obj| Box::new(obj) as Box<dyn BoundedSurface>)
        .collect();
    let bvh = ray_tracer::accelerators::BVH::build(boxed_triangles);
    let baby_yoda = Object {
        surface: bvh,
        // material: Material::dielectric(Colour::new(0.7, 0.7, 0.9), 1.5),
        material: Material::diffuse(Colour::fill(0.1)),
    };

    // Ground level to cast shadows onto
    let ground = Object {
        surface: Box::new(Plane::new(ground_level * Vec3D::Z, Vec3D::Z)),
        material: Material::diffuse(Colour::ONES * 0.8),
    };

    // Set up camera
    let angle = 0.4 * std::f32::consts::PI;
    let dir = Vec3D::new_polar(1.0, std::f32::consts::PI / 2.0, angle);
    let mut camera = ray_tracer::Camera::new(
        centre + 90.0 * dir + 40.0 * Vec3D::Z,
        Vec3D::ONES,
        ASPECT_RATIO,
    );
    camera.look_at(centre);

    let up = Vec3D::Z;
    let left = Vec3D::cross(Vec3D::Z, camera.direction);
    let mut away = camera.direction.clone();
    away.z = 0.0;

    // Set up a sun to light up the scene
    let light = Object {
        surface: Box::new(Sphere::new(
            centre + 2000.0 * up + 2000.0 * left - 1000.0 * away,
            1000.0,
        )),
        material: Material::light_source(Colour::ONES, 6.0),
    };

    // An incandescent ball
    let ball = Object {
        surface: Box::new(Sphere::new(Vec3D::new(-4.0, 6.0, 59.0), 5.0)),
        material: Material {
            albedo: Colour::new(0.8, 0.2, 0.2),
            roughness: 0.7,
            emission: Some((Colour::new(1.0, 0.03, 0.03), 1.0)),
            refractive_index: None,
        },
    };

    let scene = vec![baby_yoda, ground, light, ball];

    return (scene, camera);
}

#[allow(dead_code)]
fn scene_primitives() -> (Vec<Object>, Camera) {
    // World setup
    let mut scene = Vec::new();

    // Small objects
    scene.push(Object {
        surface: Box::new(Sphere::new(Point3D::new(-2.0, 0.0, 0.5), 0.5)),
        material: Material::diffuse(Colour::new(1.0, 0.01, 0.01)),
    });
    scene.push(Object {
        surface: Box::new(Sphere::new(Point3D::new(-0.5, 0.0, 0.8), 0.75)),
        // material: Material::diffuse(Colour::e2())
        material: Material::metal(Colour::new(1.0, 0.5, 0.5), 0.5),
    });
    scene.push(Object {
        surface: Box::new(Sphere::new(Point3D::new(1.5, 0.0, 1.0), 1.0)),
        material: Material::diffuse(Colour::new(0.01, 0.01, 1.0)),
    });
    scene.push(Object {
        surface: Box::new(Sphere::new(Point3D::new(2.5, -2.5, 0.5), 0.5)),
        material: Material::light_source(Colour::new(0.01, 1.0, 0.01), 2.0),
    });
    scene.push(Object {
        surface: Box::new(Triangle::new_lonely(
            Point3D::new(-4.5, 2.0, 0.2),
            Point3D::new(-2.0, 3.0, 0.2),
            Point3D::new(-3.0, 2.5, 3.2),
        )),
        material: Material::mirror(),
    });
    scene.push(Object {
        surface: Box::new(Sphere::new(Point3D::new(0.0, -1.6, 0.8), 0.75)),
        material: Material::dielectric(Colour::new(1.0, 1.0, 1.0), 1.5),
    });

    // Globe
    scene.push(Object {
        surface: Box::new(Sphere::new(Point3D::new(0.0, 0.0, -100.0), 100.0)),
        material: Material::diffuse(Colour::ONES * 0.5),
    });

    // Light
    scene.push(Object {
        surface: Box::new(Sphere::new(Point3D::new(0.0, 110.0, 0.0), 100.0)),
        material: Material::light_source(Colour::ONES, 3.0),
    });

    // Camera
    let camera = ray_tracer::Camera::new(
        Point3D::new(-2.5, -4.0, 3.0),
        Vec3D::new(0.5, 1.0, -0.5),
        ASPECT_RATIO,
    );

    return (scene, camera);
}
