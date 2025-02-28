use super::primitives::AxisAlignedBoundingBox as AABB;
use super::traits::{Bounded, Surface};
use super::triangle_mesh::Triangle;
use super::{Ray, Vec3D};

impl Triangle {
    #[inline]
    fn centroid_axis(&self, axis: usize) -> f32 {
        self.vertices().iter().map(|p| p[axis]).sum::<f32>() / 3.0
    }

    #[inline]
    fn centroid(&self) -> Vec3D {
        self.vertices().iter().sum::<Vec3D>() / 3.0
    }
}

/// Nodes have a bounding box, a counter for the number of primitives it contains,
/// and an index. If `num_prim == 0` then the Node is not a leaf node and `index` is
/// the index of the first child node. The second child node always has index equal
/// to `index + 1`. If on the other hand `num_prim > 0`, then the node contains
/// primitives and is a leaf node. In that case `index` is the index of the first
/// primitive.
#[derive(Debug)]
pub struct Node {
    bounding_box: AABB,
    index: u32,
    num_prim: u32,
}

impl Node {
    fn is_leaf(&self) -> bool {
        self.num_prim > 0
    }
}

/// Bounding Volume Hierarchy contains a reference to a Vec of primitives (for now
/// only triangles), and a Vec of nodes. The `indices` is a list of indices into
/// `triangles` such that if a node points to some `index`, then the corresponding
/// triangle is `triangles[indices[index]]`.
pub struct BoundingVolumeHierarchy {
    triangles: Vec<Triangle>,
    indices: Vec<usize>,
    nodes: Vec<Node>,
}

impl BoundingVolumeHierarchy {
    pub fn build(triangles: Vec<Triangle>, fast: bool) -> Self {
        // Initialise BVH.
        // The number of nodes is upper bounded by 2n-1 since the maximal number of
        // leafs is n, and for each layer of m nodes, there are a maximum of 2m
        // children. Thus, the maximum number of nodes is 1 + 2 +... + n = 2n-1.
        let n = triangles.len();
        assert!(2 * n - 1 < u32::MAX as usize);
        let mut bvh = Self {
            triangles,
            indices: (0..n).collect(),
            nodes: Vec::with_capacity(2 * n - 1),
        };

        // Add root node
        let root_node = bvh.new_node(0, n as u32);
        bvh.nodes.push(root_node);
        // Recursively subdivide the root node
        bvh.subdivide_recursive(0, fast);
        bvh.nodes.shrink_to_fit();

        return bvh;
    }

    #[inline]
    fn get_triangles(&self, first_prim: u32, num_prim: u32) -> impl Iterator<Item = &Triangle> {
        (first_prim..first_prim + num_prim).map(|i| &self.triangles[self.indices[i as usize]])
    }

    #[inline]
    fn new_node(&self, first_prim: u32, num_prim: u32) -> Node {
        let bounding_box = self
            .get_triangles(first_prim, num_prim)
            .map(|tri| tri.bounding_box())
            .fold(AABB::empty(), |a, b| AABB::combine(&a, &b));

        return Node {
            bounding_box,
            index: first_prim,
            num_prim,
        };
    }

    fn subdivide_recursive(&mut self, node_index: usize, fast: bool) {
        let partition_point;
        match if fast {
            self.midpoint_partition_point(node_index)
        } else {
            self.sah_partition_point(node_index)
        } {
            Some(i) => partition_point = i,
            None => return,
        }

        let first_prim = self.nodes[node_index].index;
        let left_count = partition_point - first_prim;
        let right_count = self.nodes[node_index].num_prim - left_count;
        if left_count == 0 || right_count == 0 {
            // One of the partitions are empty and
            // so the node cannot be subdivided.
            return;
        };

        // Update current node since it is being subdivided
        let left_child_index = self.nodes.len();
        let node = &mut self.nodes[node_index];
        node.index = left_child_index as u32;
        node.num_prim = 0;

        // Spawn children
        let left_node = self.new_node(first_prim, left_count);
        let right_node = self.new_node(partition_point, right_count);
        self.nodes.push(left_node);
        self.nodes.push(right_node);

        // Recursively subdivide child nodes
        self.subdivide_recursive(left_child_index, fast);
        self.subdivide_recursive(left_child_index + 1, fast);
    }

    fn sah_partition_point(&mut self, node_index: usize) -> Option<u32> {
        let node = &self.nodes[node_index];
        let first_prim = node.index;
        let num_prim = node.num_prim;

        // Determine split axis and split point using Surface Area Heuristic
        let mut axis = 0;
        let mut split_point = 0.0;
        let mut cost = f32::INFINITY;

        for test_axis in 0..3 {
            for triangle in self.get_triangles(first_prim, num_prim) {
                let test_point = triangle.centroid_axis(test_axis);
                let test_cost = self.evaluate_sah(node, test_axis, test_point);
                if test_cost < cost {
                    axis = test_axis;
                    split_point = test_point;
                    cost = test_cost;
                }
            }
        }

        let current_cost = num_prim as f32 * node.bounding_box.half_area();
        if cost >= current_cost {
            // Partitioning does not improve the BVH
            return None;
        }

        let i = self.partition_sort(first_prim, num_prim, axis, split_point);
        return Some(i);
    }

    fn evaluate_sah(&self, node: &Node, axis: usize, pos: f32) -> f32 {
        let mut left = AABB::empty();
        let mut right = AABB::empty();
        let mut left_count = 0;
        let mut right_count = 0;

        for triangle in self.get_triangles(node.index, node.num_prim) {
            if triangle.centroid_axis(axis) < pos {
                left_count += 1;
                left = AABB::combine(&left, &triangle.bounding_box());
            } else {
                right_count += 1;
                right = AABB::combine(&right, &triangle.bounding_box());
            }
        }

        let cost = left_count as f32 * left.half_area() + right_count as f32 * right.half_area();
        return if cost > 0.0 { cost } else { f32::INFINITY };
    }

    fn midpoint_partition_point(&mut self, node_index: usize) -> Option<u32> {
        let node = &self.nodes[node_index];
        let first_prim = node.index;
        let num_prim = node.num_prim;

        // Do not partition nodes with less than two primitives
        if num_prim <= 2 {
            return None;
        }

        // Construct centroid bounding box
        let centroid_aabb = self
            .get_triangles(first_prim, num_prim)
            .fold(AABB::empty(), |res, tri| res.grow(tri.centroid()));

        // Find largest axis and cut it in half
        let widths = centroid_aabb.upper - centroid_aabb.lower;
        let (axis, width) = (0..3)
            .map(|i| (i, widths[i]))
            .max_by(|(_, a), (_, b)| a.total_cmp(b))
            .unwrap();
        let split_point = centroid_aabb.lower[axis] + 0.5 * width;

        let i = self.partition_sort(first_prim, num_prim, axis, split_point);
        return Some(i);
    }

    fn partition_sort(
        &mut self,
        first_prim: u32,
        num_prim: u32,
        axis: usize,
        split_point: f32,
    ) -> u32 {
        // Move triangles contained in this node to match the partitioning
        let mut i = first_prim;
        let mut j = first_prim + num_prim - 1;
        while i <= j {
            if self.triangles[self.indices[i as usize]].centroid_axis(axis) < split_point {
                i += 1;
            } else {
                self.indices.swap(i as usize, j as usize);
                j -= 1;
            }
        }
        // Return the partition point (beginning of second partition)
        return i;
    }
}

impl Surface for BoundingVolumeHierarchy {
    fn hit(&self, ray: &Ray, t_min: f32, mut t_max: f32) -> Option<(f32, Vec3D)> {
        let mut hit = None;

        let mut queue: Vec<&Node> = vec![&self.nodes[0]];
        while let Some(node) = queue.pop() {
            if node.is_leaf() {
                // Intersect triangles contained in node.
                for triangle in self.get_triangles(node.index, node.num_prim) {
                    if let Some(new_hit) = triangle.hit(ray, t_min, t_max) {
                        t_max = new_hit.0;
                        hit = Some(new_hit);
                    }
                }
                continue;
            }

            // Intersect ray with both children
            let mut child1 = &self.nodes[node.index as usize];
            let mut child2 = &self.nodes[node.index as usize + 1];
            let hit1 = child1.bounding_box.hit(ray, t_min, t_max).map(|(t, _)| t);
            let hit2 = child2.bounding_box.hit(ray, t_min, t_max).map(|(t, _)| t);

            // Check if either child is hit and add them to the queue.
            match (hit1, hit2) {
                (None, None) => continue,
                (Some(_), None) => {
                    queue.push(child1);
                    continue;
                }
                (None, Some(_)) => {
                    queue.push(child2);
                    continue;
                }
                (Some(t1), Some(t2)) => {
                    // Both bounding boxes were hit => Add both to the queue, the closest
                    // is added last so that it is checked first.
                    if t1 > t2 {
                        std::mem::swap(&mut child1, &mut child2);
                    }
                    queue.push(child2);
                    queue.push(child1);
                }
            };
        }

        return hit;
    }
}
