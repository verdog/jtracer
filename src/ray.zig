const std = @import("std");

const Tuple = @import("tuple.zig").Tuple;
const Vector = @import("tuple.zig").Vector;
const Point = @import("tuple.zig").Point;

const Sphere = @import("sphere.zig").Sphere;

pub const Ray = struct {
    origin: Tuple,
    direction: Tuple,

    pub fn init(origin: Tuple, direction: Tuple) This {
        std.debug.assert(origin.isPoint());
        std.debug.assert(direction.isVector());

        return .{
            .origin = origin,
            .direction = direction,
        };
    }

    pub fn position(self: This, t: f64) Tuple {
        return self.origin.plus(self.direction.scaled(t));
    }

    const This = @This();
};

// TODO make this interface more general
pub fn intersect(sphere: Sphere, ray: Ray) []f64 {
    // TODO replace this
    var gpa_impl = std.heap.GeneralPurposeAllocator(.{}){};
    const gpa = gpa_impl.allocator();

    // for now, spheres are centered at the origin
    _ = sphere;
    const sphere_to_ray = ray.origin.minus(Point.init(0, 0, 0));

    const a = ray.direction.dot(ray.direction);
    const b = 2 * ray.direction.dot(sphere_to_ray);
    const c = sphere_to_ray.dot(sphere_to_ray) - 1;

    const discriminant = b * b - 4 * a * c;

    if (discriminant < 0) return &.{}; // XXX danger?

    var intersections = gpa.alloc(f64, 2) catch unreachable;

    intersections[0] = (-b - @sqrt(discriminant)) / (2 * a);
    intersections[1] = (-b + @sqrt(discriminant)) / (2 * a);

    return intersections;
}

const expect = std.testing.expect;

test "Ceating and querying a ray" {
    const origin = Point.init(1, 2, 3);
    const direction = Vector.init(4, 5, 6);

    const ray = Ray.init(origin, direction);

    try expect(ray.origin.equals(origin));
    try expect(ray.direction.equals(direction));
}

test "Computing a point from a distance" {
    const ray = Ray.init(Point.init(2, 3, 4), Vector.init(1, 0, 0));

    try expect(ray.position(0).equals(Point.init(2, 3, 4)));
    try expect(ray.position(1).equals(Point.init(3, 3, 4)));
    try expect(ray.position(-1).equals(Point.init(1, 3, 4)));
    try expect(ray.position(2.5).equals(Point.init(4.5, 3, 4)));
}

test "A ray intersects a sphere at two points" {
    const ray = Ray.init(Point.init(0, 0, -5), Vector.init(0, 0, 1));
    const sphere = Sphere.init();

    const xs = intersect(sphere, ray);

    try expect(xs.len == 2);
    try expect(xs[0] == 4.0);
    try expect(xs[1] == 6.0);
}

test "A ray intersects a sphere at a tangent" {
    const ray = Ray.init(Point.init(0, 1, -5), Vector.init(0, 0, 1));
    const sphere = Sphere.init();

    const xs = intersect(sphere, ray);

    try expect(xs.len == 2);
    try expect(xs[0] == 5.0);
    try expect(xs[1] == 5.0);
}

test "A ray misses a sphere" {
    const ray = Ray.init(Point.init(0, 2, -5), Vector.init(0, 0, 1));
    const sphere = Sphere.init();

    const xs = intersect(sphere, ray);

    try expect(xs.len == 0);
}

test "A ray originates inside a sphere" {
    const ray = Ray.init(Point.init(0, 0, 0), Vector.init(0, 0, 1));
    const sphere = Sphere.init();

    const xs = intersect(sphere, ray);

    try expect(xs.len == 2);
    try expect(xs[0] == -1.0);
    try expect(xs[1] == 1.0);
}

test "A sphere behind a ray" {
    const ray = Ray.init(Point.init(0, 0, 5), Vector.init(0, 0, 1));
    const sphere = Sphere.init();

    const xs = intersect(sphere, ray);

    try expect(xs.len == 2);
    try expect(xs[0] == -6.0);
    try expect(xs[1] == -4.0);
}
