use std::rc::Rc;

extern crate ray_tracer;

use ray_tracer::bvh::BoundingVolumeHierarchy as BVH;
use ray_tracer::materials::*;
use ray_tracer::primitives::*;
use ray_tracer::subdivision_surface::loop_subdivide;
use ray_tracer::triangle_mesh::{Triangle, TriangleMesh};
use ray_tracer::Camera;
use ray_tracer::{Object, Point3D, Vec3D};

const ASPECT_RATIO: f32 = 16.0 / 9.0;
const IMAGE_WIDTH: usize = 400;
const SAMPLES_PER_PIXEL: u32 = 200;
const RAY_DEPTH: u32 = 6;

fn main() -> std::io::Result<()> {
    let (scene, camera) = scene_tetrahedron();
    // Render image
    let eng = ray_tracer::engine::Engine::new(scene, camera);
    return eng.render(ASPECT_RATIO, IMAGE_WIDTH, SAMPLES_PER_PIXEL, RAY_DEPTH);
}

#[allow(dead_code)]
fn scene_tetrahedron() -> (Vec<Object>, Camera) {
    let tetrahedron = TriangleMesh {
        vertices: vec![
            Vec3D::new(1.0, 1.0, 1.0),
            Vec3D::new(1.0, -1.0, -1.0),
            Vec3D::new(-1.0, 1.0, -1.0),
            Vec3D::new(-1.0, -1.0, 1.0),
        ],
        triangles: vec![[0, 1, 2], [1, 3, 2], [0, 3, 1], [0, 2, 3]],
    };
    let tetrahedron = loop_subdivide(&tetrahedron, 3);
    let tetrahedron_material = Material::diffuse(Colour::new(1.0, 0.5, 0.5));

    // Ground to cast shadows onto
    let ground = Object {
        surface: Rc::new(Plane::new(-1.2 * Vec3D::Z, Vec3D::Z)),
        material: Material::diffuse(Colour::ONES * 0.8),
    };

    // Sun for lighting
    let light = Object {
        surface: Rc::new(Sphere::new(2000.0 * Point3D::new(1.0, 0.0, 1.0), 1000.0)),
        material: Material::light_source(Colour::ONES, 6.0),
    };

    // Camera
    let angle = -0.2 * std::f32::consts::PI;
    let loc = Vec3D::new_polar(1.0, 0.3 * std::f32::consts::PI, angle);
    let mut camera = ray_tracer::Camera::new(loc, Vec3D::ONES, ASPECT_RATIO);
    camera.look_at(Vec3D::ZERO);

    // Prepare scene for rendering
    let (_, triangles) = tetrahedron.to_triangles();
    let bvh = BVH::build(triangles);
    let tetrahedron = Object {
        surface: Rc::new(bvh),
        material: tetrahedron_material,
    };
    let scene = vec![tetrahedron, ground, light];

    return (scene, camera);
}

#[allow(dead_code)]
fn scene_teapot() -> (Vec<Object>, Camera) {
    let teapot = ray_tracer::stl::read_stl("./meshes/Utah_teapot.stl").unwrap();
    let teapot = loop_subdivide(&teapot, 1);
    let teapot_material = Material::diffuse(Colour::new(1.0, 0.5, 0.5));

    // Ground to cast shadows onto
    let ground = teapot
        .vertices
        .iter()
        .map(|v| v.z)
        .fold(f32::INFINITY, |a, b| a.min(b));
    let ground = Object {
        surface: Rc::new(Plane::new(Vec3D::Z * ground, Vec3D::Z)),
        material: Material::diffuse(Colour::ONES * 0.8),
    };

    // Sun for lighting
    let light = Object {
        surface: Rc::new(Sphere::new(2000.0 * Point3D::new(1.0, 0.0, 1.0), 1000.0)),
        material: Material::light_source(Colour::ONES, 6.0),
    };

    // Camera
    let angle = 0.5 * std::f32::consts::PI;
    let loc = Vec3D::new_polar(18.0, 0.3 * std::f32::consts::PI, angle);
    let mut camera = ray_tracer::Camera::new(loc, Vec3D::ONES, ASPECT_RATIO);
    camera.look_at(teapot.volume_centroid());

    // Prepare scene for rendering
    let (_, triangles) = teapot.to_triangles();
    let bvh = BVH::build(triangles);
    let teapot = Object {
        surface: Rc::new(bvh),
        material: teapot_material,
    };
    let scene = vec![teapot, ground, light];
    return (scene, camera);
}

#[allow(dead_code)]
fn scene_baby_yoda() -> (Vec<Object>, Camera) {
    // Baby Yoda (Grogu)
    let grogu = ray_tracer::stl::read_stl("./meshes/Baby_Yoda.stl").unwrap();
    // Unfortunately the .stl file is not good and does not define a manifold
    // mesh.try_fix_mesh();
    // let grogu = loop_subdivide(&grogu, 1);
    // let grogu_material = Material::diffuse(Colour::fill(0.1));
    let grogu_material = Material::dielectric(Colour::new(0.7, 0.7, 0.9), 1.5);

    // Ground to cast shadows onto
    let ground = grogu
        .vertices
        .iter()
        .fold(f32::INFINITY, |min, v| min.min(v.z));
    let ground = Object {
        surface: Rc::new(Plane::new(Vec3D::Z * ground, Vec3D::Z)),
        material: Material::diffuse(Colour::ONES * 0.8),
    };

    // Set up camera
    let centre = grogu.volume_centroid();
    let angle = 0.4 * std::f32::consts::PI;
    let dir = Vec3D::new_polar(1.0, std::f32::consts::PI / 2.0, angle);
    let mut camera = ray_tracer::Camera::new(
        centre + 90.0 * dir + 40.0 * Vec3D::Z,
        Vec3D::ONES,
        ASPECT_RATIO,
    );
    camera.look_at(centre);

    // Some directions
    let up = Vec3D::Z;
    let left = Vec3D::cross(Vec3D::Z, camera.direction);
    let mut away = camera.direction.clone();
    away.z = 0.0;

    // Set up a sun to light up the scene
    let light = Object {
        surface: Rc::new(Sphere::new(
            centre + 2000.0 * up + 2000.0 * left - 1000.0 * away,
            1000.0,
        )),
        material: Material::light_source(Colour::ONES, 6.0),
    };

    // An incandescent ball
    let ball = Object {
        surface: Rc::new(Sphere::new(Vec3D::new(-4.0, 6.0, 59.0), 5.0)),
        material: Material {
            albedo: Colour::new(0.8, 0.2, 0.2),
            roughness: 0.7,
            emission: Some((Colour::new(1.0, 0.03, 0.03), 1.0)),
            refractive_index: None,
        },
    };

    // Prepare scene for rendering
    let (_, triangles) = grogu.to_triangles();
    let bvh = BVH::build(triangles);
    let grogu = Object {
        surface: Rc::new(bvh),
        material: grogu_material,
    };
    let scene = vec![grogu, ground, light, ball];

    return (scene, camera);
}

#[allow(dead_code)]
fn scene_primitives() -> (Vec<Object>, Camera) {
    // World setup
    let mut scene = Vec::new();

    // Small objects
    scene.push(Object {
        surface: Rc::new(Sphere::new(Point3D::new(-2.0, 0.0, 0.5), 0.5)),
        material: Material::diffuse(Colour::new(1.0, 0.01, 0.01)),
    });
    scene.push(Object {
        surface: Rc::new(Sphere::new(Point3D::new(-0.5, 0.0, 0.8), 0.75)),
        // material: Material::diffuse(Colour::e2())
        material: Material::metal(Colour::new(1.0, 0.5, 0.5), 0.5),
    });
    scene.push(Object {
        surface: Rc::new(Sphere::new(Point3D::new(1.5, 0.0, 1.0), 1.0)),
        material: Material::diffuse(Colour::new(0.01, 0.01, 1.0)),
    });
    scene.push(Object {
        surface: Rc::new(Sphere::new(Point3D::new(2.5, -2.5, 0.5), 0.5)),
        material: Material::light_source(Colour::new(0.01, 1.0, 0.01), 2.0),
    });
    scene.push(Object {
        surface: Rc::new(Triangle::new_lonely(
            Point3D::new(-4.5, 2.0, 0.2),
            Point3D::new(-2.0, 3.0, 0.2),
            Point3D::new(-3.0, 2.5, 3.2),
        )),
        material: Material::mirror(),
    });
    scene.push(Object {
        surface: Rc::new(Sphere::new(Point3D::new(0.0, -1.6, 0.8), 0.75)),
        material: Material::dielectric(Colour::new(1.0, 1.0, 1.0), 1.5),
    });

    // Globe
    scene.push(Object {
        surface: Rc::new(Sphere::new(Point3D::new(0.0, 0.0, -100.0), 100.0)),
        material: Material::diffuse(Colour::ONES * 0.5),
    });

    // Light
    scene.push(Object {
        surface: Rc::new(Sphere::new(Point3D::new(0.0, 110.0, 0.0), 100.0)),
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
