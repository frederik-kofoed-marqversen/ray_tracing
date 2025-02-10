use std::collections::HashMap;

use super::{Point3D, Vec3D};

#[derive(Clone)]
pub struct TriangleMesh {
    vertices: Vec<Point3D>,
    triangles: Vec<[usize; 3]>,
}

/// This is an implementation of Loop Subdivision Surface. The basic ideas are from
/// https://pbr-book.org/3ed-2018/Shapes/Subdivision_Surfaces#LoopSubdiv::beta
/// However, major restructuring and optimisations have been done. They also describe 
/// how to compute vertex normals. This has not been implemented here as of yet.
pub fn loop_subdivide(mesh: &TriangleMesh, divisions: usize) -> TriangleMesh {
    let mut sd_surface = SDSurface::from_triangle_mesh(mesh);

    for _ in 0..divisions {
        sd_surface.subdivide();
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

struct SDVertex {
    point: Vec3D,
    reference_face: usize, // index of any incident face
}

/// Represents a triangular face
///
/// # Ordering of Neighbours
/// - The `neighbours` are ordered according to the cyclic ordering of vertices:
///   - If `neighbours[i]` is `Some(j)`, then the `j`th face shares the edge given by
///     (`vertices[i]`, `vertices[next(i)]`).
///   - If `neighbours[i]` is `None`, that edge is a **boundary edge**.
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
        let mut edges: HashMap<(usize, usize), (usize, usize)> = HashMap::new();
        let faces = &mut sd_surface.faces;
        for f_index in 0..faces.len() {
            for e_index in 0..3 {
                let edge = (
                    faces[f_index].vertices[e_index],
                    faces[f_index].vertices[next(e_index)],
                );
                // Assuming consistent orientation of triangles over the entire mesh. Thus,
                // when two faces share an edge, the edge is given in opposite direction.
                if let Some((neighbour, neighbour_e_index)) = edges.remove(&(edge.1, edge.0)) {
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

    fn subdivide(&mut self) {
        // Initialize new faces. Since each face has 4 children the
        // children of face `i` will be at indices 4*i..4*(i+1)
        let mut new_faces: Vec<SDFace> = (0..4 * self.faces.len())
            .map(|_| SDFace {
                vertices: [0; 3],
                neighbours: [None; 3],
            })
            .collect();

        // Update positions of existing vertices
        for v_index in 0..self.vertices.len() {
            let vertex = &self.vertices[v_index];
            let (neighbours, on_boundary) = self.vertex_neighbours(v_index);
            let valence = neighbours.len();

            let new_point;
            if !on_boundary {
                let beta = if valence == 6 {
                    1.0 / 16.0
                } else if valence == 3 {
                    3.0 / 16.0
                } else {
                    3.0 / (8 * valence) as f32
                };

                let neighbour_sum: Point3D =
                    neighbours.iter().map(|&i| self.vertices[i].point).sum();
                new_point = (1.0 - valence as f32 * beta) * vertex.point + beta * neighbour_sum;
            } else {
                let beta = 1.0 / 8.0;
                let neighbour_sum = self.vertices[neighbours[0]].point
                    + self.vertices[neighbours[valence - 1]].point;
                new_point = (1.0 - 2.0 * beta) * vertex.point + beta * neighbour_sum;
            }

            // An existing vertex which is the i'th vertex of a face will always be incident to
            // the i'th child face. The children of the j'th fase begin at new index 4*j
            let i = self.faces[vertex.reference_face].vertex_index(v_index);
            let new_reference_face = 4 * vertex.reference_face + i;

            // Update vertex. The on_boundary field is preserved.
            let vertex = &mut self.vertices[v_index];
            vertex.reference_face = new_reference_face;
            vertex.point = new_point
        }

        // Compute new additional vertices
        // `edges` is a map from edges to vertex index in `new_vertices`.
        let mut edges: HashMap<(usize, usize), usize> = HashMap::new();
        for f_index in 0..self.faces.len() {
            let face = &self.faces[f_index];
            for i in 0..3 {
                // Edges will be ordered pairs.
                let mut edge = [face.vertices[i], face.vertices[next(i)]];
                edge.sort();
                let edge = (edge[0], edge[1]);
                if edges.contains_key(&edge) {
                    // Vertex has already been created
                    continue;
                }

                // Record that vertex on this edge will be created
                edges.insert(edge, self.vertices.len());

                // Added vertices are always incident to the central 4th child of a face.
                let reference_face = 4 * (f_index + 1);
                // Compute new vertex and add to vertex list
                if let Some(f2_index) = face.neighbours[i] {
                    // Edge is NOT a boundary edge => created vertex is NOT on the boundary
                    let face2 = &self.faces[f2_index];
                    let point = 3.0 / 8.0 * self.vertices[edge.0].point
                        + 3.0 / 8.0 * self.vertices[edge.1].point
                        + 1.0 / 8.0 * self.vertices[face.vertices[prev(i)]].point
                        + 1.0 / 8.0 * self.vertices[face2.vertices[prev(i)]].point;
                    self.vertices.push(SDVertex {
                        point,
                        reference_face,
                    });
                } else {
                    // Edge is a boundary edge => created vertex is on the boundary
                    let point = 0.5 * (self.vertices[edge.0].point + self.vertices[edge.1].point);
                    self.vertices.push(SDVertex {
                        point,
                        reference_face,
                    });
                }
            }
        }

        // Compute new topology (face data)
        // Loop over each of the parent faces and make sure that its' children are properly updated
        for (f_index, face) in self.faces.iter().enumerate() {
            let child = |i: usize| -> usize { 4 * f_index + i };
            // Loop over the three corner vertices.
            for i in 0..3 {
                // The corner vertex at index `i` belongs to `child(i)`
                new_faces[child(i)].vertices[i] = face.vertices[i];
                // For each old vertex a new vertex has been added along the edge (`i`, `next(i)`)
                // This vertex is belongs to three children of the current face
                let vertex = *edges
                    .get(&(face.vertices[i], face.vertices[next(i)]))
                    .unwrap();
                new_faces[child(i)].vertices[next(i)] = vertex;
                new_faces[child(next(i))].vertices[i] = vertex;
                new_faces[child(3)].vertices[i] = vertex;

                // Update the `neighbours` field of `child(i)`
                // Central child `child(3)` is a neighbour
                new_faces[child(3)].neighbours[prev(i)] = Some(child(i));
                new_faces[child(i)].neighbours[next(i)] = Some(child(3));

                // There are two potential neighbouring parents. These might not exist due to borders.
                // If a parent exist, they share the corner vertex. The index of a corner vertex in a face
                // is also the index of the child in that corner.
                if let Some(f2_index) = face.neighbours[i] {
                    let corner_vertex = face.vertices[i];
                    let f2_i = &self.faces[f2_index].vertex_index(corner_vertex);
                    let f2_child = 4 * f2_index + f2_i;
                    new_faces[child(i)].neighbours[i] = Some(f2_child);
                }
                if let Some(f2_index) = face.neighbours[prev(i)] {
                    let corner_vertex = face.vertices[i];
                    let f2_i = &self.faces[f2_index].vertex_index(corner_vertex);
                    let f2_child = 4 * f2_index + f2_i;
                    new_faces[child(i)].neighbours[prev(i)] = Some(f2_child);
                }
            }
        }

        // Replace old topology with new
        self.faces = new_faces;
    }

    fn push_to_limit_surface(&mut self) {
        for vertex in 0..self.vertices.len() {
            let (neighbours, on_boundary) = self.vertex_neighbours(vertex);
            let valence = neighbours.len();
            let point = self.vertices[vertex].point;
            if !on_boundary {
                let beta = if valence == 3 {
                    3.0 / 16.0
                } else {
                    3.0 / (8.0 * valence as f32)
                };
                let beta = 1.0 / (valence as f32 + 3.0 / (8.0 * beta));

                let neighbour_sum: Point3D =
                    neighbours.iter().map(|&i| self.vertices[i].point).sum();
                self.vertices[vertex].point =
                    (1.0 - valence as f32 * beta) * point + beta * neighbour_sum;
            } else {
                let beta = 1.0 / 5.0;
                let neighbour_sum = self.vertices[neighbours[0]].point
                    + self.vertices[neighbours[valence - 1]].point;
                self.vertices[vertex].point = (1.0 - 2.0 * beta) * point + beta * neighbour_sum;
            }
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

    fn next_vertex(&self, vertex: usize, face: usize) -> usize {
        let face = &self.faces[face];
        face.vertices[next(face.vertex_index(vertex))]
    }

    fn prev_vertex(&self, vertex: usize, face: usize) -> usize {
        let face = &self.faces[face];
        face.vertices[prev(face.vertex_index(vertex))]
    }

    // Given a vertex `i` of a face the next face is defined as that which shares the edge
    // (`i`, `next(i)`) which is simply the face `neighbours[i]`.
    fn next_face(&self, vertex: usize, face: usize) -> Option<usize> {
        let face = &self.faces[face];
        face.neighbours[face.vertex_index(vertex)]
    }

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
                // Has returned to starting point, and thus not a boundary point
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
        while let Some(next_face) = self.prev_face(vertex, face) {
            face = next_face;
            neighbours.push(self.prev_vertex(vertex, next_face));
        }

        return (neighbours, true);
    }
}
