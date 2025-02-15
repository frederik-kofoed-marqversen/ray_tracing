use core::f32;
use std::collections::HashMap;

use super::triangle_mesh::TriangleMesh;
use super::{Point3D, Vec3D};

/// This is an implementation of Loop Subdivision Surface. The algorithm is laid out in
/// Hoppe et al. 1994. "Piecewise smooth surface reconstruction". In Proceedings of
/// SIGGRAPH ’94, Computer Graphics Proceedings, Annual Conference Series, Orlando,
/// Florida, 295–302.
///
/// To properly understand the article i was aided by:
/// https://pbr-book.org/3ed-2018/Shapes/Subdivision_Surfaces#LoopSubdiv::beta
///
/// This is a minimal implementation. As discribed in the paper, much more is possible like
/// computing vertex normals, and tangent vectors to make a UV mesh.
pub fn loop_subdivide(mesh: &TriangleMesh, divisions: usize) -> TriangleMesh {
    let mut sd_surface = SDSurface::from_triangle_mesh(mesh);
    for _ in 0..divisions {
        sd_surface = sd_surface.subdivide();
    }
    sd_surface.push_to_limit_surface();
    return sd_surface.to_triangle_mesh();
}

/// Next index in cyclic ordering on three elements
#[inline]
fn next(i: usize) -> usize {
    (i + 1) % 3
}

/// Previous index in cyclic ordering on three elements
#[inline]
fn prev(i: usize) -> usize {
    (i + 2) % 3
}

#[derive(Debug)]
struct SDVertex {
    point: Vec3D,
    reference_face: usize, // index of any incident face
}

#[derive(Hash, PartialEq, Eq, Clone)]
struct SDEdge {
    v1: usize,
    v2: usize,
}

impl SDEdge {
    fn new(v1: usize, v2: usize) -> Self {
        let mut vertices = [v1, v2];
        vertices.sort();
        return Self {
            v1: vertices[0],
            v2: vertices[1],
        };
    }
}

/// Represents a triangular face
///
/// # Ordering of Neighbours
/// - The `neighbours` are ordered according to the positive cyclic ordering of vertices:
///   - If `neighbours[i]` is `Some(j)`, then the `j`th face shares the edge given by
///     (`vertices[i]`, `vertices[next(i)]`).
///   - If `neighbours[i]` is `None`, that edge is a **boundary edge**.
#[derive(Debug)]
struct SDFace {
    vertices: [usize; 3],           // index of vertices in positive order
    neighbours: [Option<usize>; 3], // indices of neighbouring faces.
}

impl SDFace {
    // Get the index into `vertices` of a given vertex.
    fn vertex_index(&self, vertex: usize) -> usize {
        if let Some(index) = self.vertices.iter().position(|other| other == &vertex) {
            return index;
        } else {
            // We should never ask for a vertex that is not part of a face
            panic!("Logic error!")
        }
    }
}

#[derive(Debug)]
struct SDSurface {
    vertices: Vec<SDVertex>,
    faces: Vec<SDFace>,
}

impl SDSurface {
    fn from_triangle_mesh(mesh: &TriangleMesh) -> Self {
        let mut vertices = Vec::with_capacity(mesh.vertices.len());
        let mut faces = Vec::with_capacity(mesh.triangles.len());

        // Initialise vec of vertices with their point in 3D space. Remaining
        // fields will be computed and correctly set later.
        for vertex in &mesh.vertices {
            vertices.push(SDVertex {
                point: vertex.clone(),
                reference_face: 0,
            });
        }

        // Initialise vec of triangle faces by setting the indices into the
        // vec of vertices. Neighbours will be set later
        for triangle in &mesh.triangles {
            faces.push(SDFace {
                vertices: triangle.clone(),
                neighbours: [None; 3],
            })
        }

        // Construct SDSurface object
        let mut sd_surface = Self { vertices, faces };

        // As reference face for each vertex set the index of any neighbouring face.
        for (f_index, face) in sd_surface.faces.iter().enumerate() {
            for v_index in face.vertices {
                sd_surface.vertices[v_index].reference_face = f_index;
            }
        }

        // Compute face neighbours
        // An edge is given by two vertex indices in positive direction as defined by the mesh.
        // When encountering an edge, it is stored together with its corresponding face, and
        // the index of that edge into the neighbours field of that face.
        // When we are done, `edges` contains all edges that are boundary edges.
        let mut edges: HashMap<SDEdge, (usize, usize)> = HashMap::new();
        let faces = &mut sd_surface.faces;
        for f_index in 0..faces.len() {
            for e_index in 0..3 {
                let edge = SDEdge::new(
                    faces[f_index].vertices[e_index],
                    faces[f_index].vertices[next(e_index)],
                );
                // Assuming consistent orientation of triangles over the entire mesh. Thus,
                // at most two faces share an edge.
                if let Some((neighbour, neighbour_e_index)) = edges.remove(&edge) {
                    // Edge has been seen before and so is removed and neighbouring faces are updated.
                    faces[f_index].neighbours[e_index] = Some(neighbour);
                    faces[neighbour].neighbours[neighbour_e_index] = Some(f_index);
                } else {
                    // Edge has not been seen before.
                    edges.insert(edge, (f_index, e_index));
                }
            }
        }

        return sd_surface;
    }

    fn subdivide(&self) -> Self {
        // Vec for storing the new vertices. There is approximately F ≈ 2*V in a triangle mesh.
        // The new mesh will have 4 * F faces, and so approximately 2 * F vertices!
        let mut new_vertices: Vec<SDVertex> = Vec::with_capacity(2 * self.faces.len());

        // Calculate updated positions for all existing vertices and add them to the list.
        // Ordering of already existing vertices will be preserved.
        for (v_index, vertex) in self.vertices.iter().enumerate() {
            let (neighbours, on_boundary) = self.vertex_neighbours(v_index);
            let valence = neighbours.len();

            let new_point;
            if !on_boundary {
                let n = valence as f32;
                let a =
                    5.0 / 8.0 - (3.0 + 2.0 * f32::cos(2.0 * f32::consts::PI / n)).powi(2) / 64.0;
                let alpha = n * (1.0 - a) / a;
                let neighbour_sum: Point3D =
                    neighbours.iter().map(|&i| self.vertices[i].point).sum();
                new_point = (alpha * vertex.point + neighbour_sum) / (alpha + n);
            } else {
                let neighbour_sum = self.vertices[neighbours[0]].point
                    + self.vertices[neighbours[valence - 1]].point;
                new_point = (6.0 * vertex.point + neighbour_sum) / 8.0;
            }

            // An existing vertex which is the i'th vertex of a face will always be incident to
            // the i'th child face. The children of the j'th fase begin at new index 4*j
            let i = self.faces[vertex.reference_face].vertex_index(v_index);
            let new_reference_face = 4 * vertex.reference_face + i;

            // Init new updated vertex
            new_vertices.push(SDVertex {
                point: new_point,
                reference_face: new_reference_face,
            });
        }

        // Compute new additional vertices which all lie on an existing edge. Notice that
        // positions of new vertices are based on old positions of already existing vertices.
        // `edges` is a map from edges to vertex index in `new_vertices`.
        let mut edges: HashMap<SDEdge, usize> = HashMap::new();
        for (f_index, face) in self.faces.iter().enumerate() {
            for i in 0..3 {
                let edge = SDEdge::new(face.vertices[i], face.vertices[next(i)]);
                if edges.contains_key(&edge) {
                    // Vertex has already been created
                    continue;
                }

                // Record that vertex on this edge will be created
                edges.insert(edge.clone(), new_vertices.len());

                // Added vertices are always incident to the central 4th child of a face.
                let reference_face = 4 * f_index + 3;
                // Compute new vertex and add to vertex list
                if let Some(f2_index) = face.neighbours[i] {
                    // Edge is NOT a boundary edge => created vertex is NOT on the boundary
                    let face2 = &self.faces[f2_index];
                    let point = 3.0 / 8.0 * self.vertices[edge.v1].point
                        + 3.0 / 8.0 * self.vertices[edge.v2].point
                        + 1.0 / 8.0 * self.vertices[face.vertices[prev(i)]].point
                        + 1.0 / 8.0 * self.vertices[face2.vertices[prev(i)]].point;
                    new_vertices.push(SDVertex {
                        point,
                        reference_face,
                    });
                } else {
                    // Edge is a boundary edge => created vertex is on the boundary
                    let point = 0.5 * (self.vertices[edge.v1].point + self.vertices[edge.v2].point);
                    new_vertices.push(SDVertex {
                        point,
                        reference_face,
                    });
                }
            }
        }

        // Initialize new faces. Since each face has 4 children the
        // children of face `i` will be at indices 4*i..4*(i+1)
        let mut new_faces: Vec<SDFace> = (0..4 * self.faces.len())
            .map(|_| SDFace {
                vertices: [0; 3],
                neighbours: [None; 3],
            })
            .collect();

        // Compute new topology (face data). Figure 3.36 defines the ordering of these things.
        // Loop over each of the parent faces and make sure that its' children are properly updated
        for (f_index, face) in self.faces.iter().enumerate() {
            let child = |i: usize| -> usize { 4 * f_index + i };
            // Loop over the three corner vertices.
            for i in 0..3 {
                // The corner vertex at index `i` belongs to `child(i)`
                new_faces[child(i)].vertices[i] = face.vertices[i];
                // For each old vertex a new vertex has been added along the edge (`i`, `next(i)`)
                // This vertex belongs to three children of the current face
                let edge = SDEdge::new(face.vertices[i], face.vertices[next(i)]);
                let vertex = *edges.get(&edge).unwrap();
                new_faces[child(i)].vertices[next(i)] = vertex;
                new_faces[child(next(i))].vertices[i] = vertex;
                new_faces[child(3)].vertices[i] = vertex;

                // Update the `neighbours` field of `child(i)`
                // Central child `child(3)` is a neighbour
                new_faces[child(i)].neighbours[next(i)] = Some(child(3));
                new_faces[child(3)].neighbours[prev(i)] = Some(child(i));

                // There are two potential neighbouring parents. These might not exist due to borders.
                // If a parent exist, they share the corner vertex.
                if let Some(f2_index) = face.neighbours[i] {
                    let corner_vertex = face.vertices[i];
                    let f2_i = self.faces[f2_index].vertex_index(corner_vertex);
                    let f2_child = 4 * f2_index + f2_i;
                    new_faces[child(i)].neighbours[i] = Some(f2_child);
                }
                if let Some(f2_index) = face.neighbours[prev(i)] {
                    let corner_vertex = face.vertices[i];
                    let f2_i = self.faces[f2_index].vertex_index(corner_vertex);
                    let f2_child = 4 * f2_index + f2_i;
                    new_faces[child(i)].neighbours[prev(i)] = Some(f2_child);
                }
            }
        }

        // Return new mesh
        return SDSurface {
            vertices: new_vertices,
            faces: new_faces,
        };
    }

    fn push_to_limit_surface(&mut self) {
        // Save updated points in new vector since update depends on old points.
        let mut limit_points = Vec::with_capacity(self.vertices.len());
        for vertex in 0..self.vertices.len() {
            let (neighbours, on_boundary) = self.vertex_neighbours(vertex);
            let valence = neighbours.len();
            let point = self.vertices[vertex].point;
            if !on_boundary {
                let n = valence as f32;
                let a =
                    5.0 / 8.0 - (3.0 + 2.0 * f32::cos(2.0 * f32::consts::PI / n)).powi(2) / 64.0;
                let omega = 3.0 * n / (8.0 * a);

                let neighbour_sum: Point3D =
                    neighbours.iter().map(|&i| self.vertices[i].point).sum();
                let new_point = (omega * point + neighbour_sum) / (omega + n);
                limit_points.push(new_point);
            } else {
                // For simplicity we always use regular crease vertex weights.
                let neighbour_sum = self.vertices[neighbours[0]].point
                    + self.vertices[neighbours[valence - 1]].point;
                let new_point = (3.0 * point + neighbour_sum) / 5.0;
                limit_points.push(new_point);
            }
        }
        // Set vertex points to their new updated points
        for vertex in 0..self.vertices.len() {
            self.vertices[vertex].point = limit_points[vertex];
        }
    }

    fn to_triangle_mesh(&self) -> TriangleMesh {
        let vertices = self.vertices.iter().map(|v| v.point).collect();
        let triangles = self.faces.iter().map(|f| f.vertices).collect();
        return TriangleMesh {
            vertices,
            triangles,
        };
    }

    #[inline]
    fn next_vertex(&self, vertex: usize, face: usize) -> usize {
        let face = &self.faces[face];
        face.vertices[next(face.vertex_index(vertex))]
    }

    #[inline]
    fn prev_vertex(&self, vertex: usize, face: usize) -> usize {
        let face = &self.faces[face];
        face.vertices[prev(face.vertex_index(vertex))]
    }

    // Given a vertex `i` of a face the next face is defined as that which shares the edge
    // (`i`, `next(i)`) which is simply the face `neighbours[i]`.
    #[inline]
    fn next_face(&self, vertex: usize, face: usize) -> Option<usize> {
        let face = &self.faces[face];
        face.neighbours[face.vertex_index(vertex)]
    }

    #[inline]
    fn prev_face(&self, vertex: usize, face: usize) -> Option<usize> {
        let face = &self.faces[face];
        face.neighbours[prev(face.vertex_index(vertex))]
    }

    /// Computes the neighbouring vertices for a given vertex along with a bool
    /// which is true if the given vertex is a boundary vertex and false if not.
    fn vertex_neighbours(&self, vertex: usize) -> (Vec<usize>, bool) {
        let start_face = self.vertices[vertex].reference_face;
        let mut neighbours = Vec::new();

        // Loop through incident faces until return to start or hit boundary.
        // For each face, store the next vertex along that face.
        let mut face = start_face;
        neighbours.push(self.next_vertex(vertex, face));
        while let Some(next_face) = self.next_face(vertex, face) {
            if next_face == start_face {
                // Has returned to starting point, and vertex is thus not a boundary point
                return (neighbours, false);
            } else {
                face = next_face;
                neighbours.push(self.next_vertex(vertex, next_face));
            }
        }

        // A boundary must have been hit.
        // Reverse the order of neighbours so that boundary vertices will be at either end
        neighbours.reverse();
        // Start again but go in opposite direction.
        face = start_face;
        neighbours.push(self.prev_vertex(vertex, face));
        while let Some(prev_face) = self.prev_face(vertex, face) {
            face = prev_face;
            neighbours.push(self.prev_vertex(vertex, prev_face));
        }

        return (neighbours, true);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn vertex_neighbours() {
        let mesh = TriangleMesh {
            vertices: vec![Vec3D::default(); 5],
            triangles: vec![[0, 1, 4], [1, 2, 4], [2, 3, 4], [3, 0, 4]],
        };
        let sds = SDSurface::from_triangle_mesh(&mesh);

        // Check neighbours of boundary vertices
        let neighbours = [vec![1, 4, 3], vec![2, 4, 0], vec![3, 4, 1], vec![0, 4, 2]];
        for (i, neighbours) in neighbours.into_iter().enumerate() {
            assert_eq!(sds.vertex_neighbours(i), (neighbours, true));
        }

        // Check neighbours of internal vertex
        let (mut neighbours, boundary) = sds.vertex_neighbours(4);
        neighbours.sort();
        assert_eq!((neighbours, boundary), (vec![0, 1, 2, 3], false));
    }

    #[test]
    fn subdivision() {
        let mesh = TriangleMesh {
            vertices: vec![Vec3D::default(); 3],
            triangles: vec![[0, 1, 2]],
        };

        let sds = SDSurface::from_triangle_mesh(&mesh);
        let sds = sds.subdivide();

        // Check vertices
        for i in 0..3 {
            assert_eq!(sds.vertices[i].reference_face, i);
        }
        for i in 3..6 {
            assert_eq!(sds.vertices[i].reference_face, 3);
        }

        // Check faces
        let checks = [
            ([0, 3, 5], [None, Some(3), None]),
            ([3, 1, 4], [None, None, Some(3)]),
            ([5, 4, 2], [Some(3), None, None]),
            ([3, 4, 5], [Some(1), Some(2), Some(0)]),
        ];
        for i in 0..4 {
            let face = &sds.faces[i];
            assert_eq!((face.vertices, face.neighbours), checks[i]);
        }
    }

    #[test]
    fn topology() {
        use super::super::stl::read_stl;
        use std::collections::HashSet;

        let mesh = read_stl("../src/meshes/Baby_Yoda.stl").unwrap();
        let sds = SDSurface::from_triangle_mesh(&mesh);

        let faces: Vec<(usize, &SDFace)> = sds
            .faces
            .iter()
            .enumerate()
            .filter(|(_i, face)| face.vertices.contains(&194))
            .collect();

        let vertices: HashSet<&usize> = faces
            .iter()
            .flat_map(|(_i, face)| face.vertices.iter())
            .collect();

        dbg!(faces);
        dbg!(vertices);
    }
}
