//! intersection algorithms

const std = @import("std");

const Tuple = @import("tuple.zig").Tuple;
const Vector = @import("tuple.zig").Vector;
const Point = @import("tuple.zig").Point;
const Ray = @import("ray.zig").Ray;

const VolPtr = @import("world.zig").VolumePtr;
const Sphere = @import("sphere.zig").Sphere;

const trans = @import("transform.zig");
const mtx = @import("matrix.zig");

pub const Intersection = struct {
    t: f64,
    vptr: VolPtr,

    pub fn init(t: f64, vptr: VolPtr) This {
        return .{ .t = t, .vptr = vptr };
    }

    pub fn equals(self: This, other: This) bool {
        return self.t == other.t and std.meta.eql(self.vptr, other.vptr);
    }

    const This = @This();
};

pub const Intersections = struct {
    ixs: std.ArrayList(Intersection),

    pub fn init(alctr: std.mem.Allocator, inits: anytype) This {
        const InitsT = @TypeOf(inits);
        const inits_type_info = @typeInfo(InitsT);
        if (inits_type_info != .Struct) {
            @compileError("Init with .{} syntax");
        }
        const fields_info = inits_type_info.Struct.fields;

        var xs = This{
            .ixs = std.ArrayList(Intersection).init(alctr),
        };

        inline for (fields_info) |field| {
            xs.ixs.append(@field(inits, field.name)) catch unreachable;
        }

        return xs;
    }

    pub fn deinit(self: This) void {
        self.ixs.deinit();
    }

    pub fn hit(self: This) ?Intersection {
        // TODO eventually this should turn into a full-on sort
        var min_positive_t: ?f64 = null;
        var min_positive_i: usize = self.ixs.items.len;

        for (self.ixs.items) |x, i| {
            if (x.t >= 0 and x.t < min_positive_t orelse std.math.floatMax(f64)) {
                min_positive_t = x.t;
                min_positive_i = i;
            }
        }

        if (min_positive_i < self.ixs.items.len)
            return self.ixs.items[min_positive_i]
        else
            return null;
    }

    const This = @This();
};

// TODO make this interface more general
pub fn intersect(vptr: VolPtr, sphere: Sphere, ray: Ray, alctr: std.mem.Allocator) Intersections {
    const td_ray = ray.transformed(sphere.transform.inverted() catch unreachable);
    const sphere_to_ray = td_ray.origin.minus(Point.init(0, 0, 0));

    const a = td_ray.direction.dot(td_ray.direction);
    const b = 2 * td_ray.direction.dot(sphere_to_ray);
    const c = sphere_to_ray.dot(sphere_to_ray) - 1;

    const discriminant = b * b - 4 * a * c;

    if (discriminant < 0) return Intersections.init(alctr, .{});

    var intersections = Intersections.init(alctr, .{
        Intersection.init((-b - @sqrt(discriminant)) / (2 * a), vptr),
        Intersection.init((-b + @sqrt(discriminant)) / (2 * a), vptr),
    });

    return intersections;
}

const expect = std.testing.expect;

test "An intersection encapsulates t and object" {
    const sphereidx = VolPtr{ .sphere_idx = 0 };
    const t: f64 = 3.5;

    const intx = Intersection.init(t, sphereidx);

    try expect(intx.t == 3.5);
    try expect(intx.vptr.sphere_idx == sphereidx.sphere_idx);
}

test "Aggregating intersections" {
    const alctr = std.testing.allocator;
    const sphereidx = VolPtr{ .sphere_idx = 0 };
    const int1 = Intersection.init(1.0, sphereidx);
    const int2 = Intersection.init(2.0, sphereidx);

    const xs = Intersections.init(alctr, .{ int1, int2 });
    defer xs.deinit();

    try expect(xs.ixs.items.len == 2);
    try expect(xs.ixs.items[0].t == 1.0);
    try expect(xs.ixs.items[1].t == 2.0);
}

test "A ray intersects a sphere at two points" {
    const alctr = std.testing.allocator;
    const ray = Ray.init(Point.init(0, 0, -5), Vector.init(0, 0, 1));
    const sphere = Sphere.init();
    const vptr = VolPtr{ .sphere_idx = 0 };

    const xs = intersect(vptr, sphere, ray, alctr);
    defer xs.deinit();

    try expect(xs.ixs.items.len == 2);
    try expect(xs.ixs.items[0].t == 4.0);
    try expect(xs.ixs.items[1].t == 6.0);
}

test "A ray intersects a sphere at a tangent" {
    const alctr = std.testing.allocator;
    const ray = Ray.init(Point.init(0, 1, -5), Vector.init(0, 0, 1));
    const sphere = Sphere.init();
    const vptr = VolPtr{ .sphere_idx = 0 };

    const xs = intersect(vptr, sphere, ray, alctr);
    defer xs.deinit();

    try expect(xs.ixs.items.len == 2);
    try expect(xs.ixs.items[0].t == 5.0);
    try expect(xs.ixs.items[1].t == 5.0);
}

test "A ray misses a sphere" {
    const alctr = std.testing.allocator;
    const ray = Ray.init(Point.init(0, 2, -5), Vector.init(0, 0, 1));
    const sphere = Sphere.init();
    const vptr = VolPtr{ .sphere_idx = 0 };

    const xs = intersect(vptr, sphere, ray, alctr);

    try expect(xs.ixs.items.len == 0);
}

test "A ray originates inside a sphere" {
    const alctr = std.testing.allocator;
    const ray = Ray.init(Point.init(0, 0, 0), Vector.init(0, 0, 1));
    const sphere = Sphere.init();
    const vptr = VolPtr{ .sphere_idx = 0 };

    const xs = intersect(vptr, sphere, ray, alctr);
    defer xs.deinit();

    try expect(xs.ixs.items.len == 2);
    try expect(xs.ixs.items[0].t == -1.0);
    try expect(xs.ixs.items[1].t == 1.0);
}

test "A sphere behind a ray" {
    const alctr = std.testing.allocator;
    const ray = Ray.init(Point.init(0, 0, 5), Vector.init(0, 0, 1));
    const sphere = Sphere.init();
    const vptr = VolPtr{ .sphere_idx = 0 };

    const xs = intersect(vptr, sphere, ray, alctr);
    defer xs.deinit();

    try expect(xs.ixs.items.len == 2);
    try expect(xs.ixs.items[0].t == -6.0);
    try expect(xs.ixs.items[1].t == -4.0);
}

test "Intersect a scaled sphere with a ray" {
    const alctr = std.testing.allocator;
    const ray = Ray.init(Point.init(0, 0, -5), Vector.init(0, 0, 1));
    var s = Sphere.init();
    const vptr = VolPtr{ .sphere_idx = 7 };

    s.transform = trans.makeScaling(2, 2, 2);

    const xs = intersect(vptr, s, ray, alctr);
    defer xs.deinit();

    try expect(xs.ixs.items.len == 2);
    try expect(xs.ixs.items[0].t == 3.0);
    try expect(xs.ixs.items[1].t == 7.0);
}

test "Intersect a translated sphere with a ray" {
    const alctr = std.testing.allocator;
    const ray = Ray.init(Point.init(0, 0, -5), Vector.init(0, 0, 1));
    var s = Sphere.init();
    const vptr = VolPtr{ .sphere_idx = 7 };

    s.transform = trans.makeTranslation(5, 0, 0);

    const xs = intersect(vptr, s, ray, alctr);
    defer xs.deinit();

    try expect(xs.ixs.items.len == 0);
}

test "Intersect sets the object on the intersection" {
    const alctr = std.testing.allocator;
    const ray = Ray.init(Point.init(0, 0, -5), Vector.init(0, 0, 1));
    const sph = Sphere.init();
    const vptr = VolPtr{ .sphere_idx = 7 };

    const xs = intersect(vptr, sph, ray, alctr);
    defer xs.deinit();

    try expect(xs.ixs.items.len == 2);
    try expect(xs.ixs.items[0].vptr.sphere_idx == 7);
    try expect(xs.ixs.items[1].vptr.sphere_idx == 7);
}

test "The hit, where all intersections have positive t" {
    const alctr = std.testing.allocator;
    const int1 = Intersection.init(1.0, VolPtr{ .sphere_idx = 5 });
    const int2 = Intersection.init(2.0, VolPtr{ .sphere_idx = 5 });
    const xs = Intersections.init(alctr, .{ int1, int2 });
    defer xs.deinit();

    try expect(xs.hit().?.equals(int1));
}

test "The hit, where some intersections have negative t" {
    const alctr = std.testing.allocator;
    const int1 = Intersection.init(-1.0, VolPtr{ .sphere_idx = 5 });
    const int2 = Intersection.init(1.0, VolPtr{ .sphere_idx = 5 });
    const xs = Intersections.init(alctr, .{ int1, int2 });
    defer xs.deinit();

    try expect(xs.hit().?.equals(int2));
}

test "The hit, where all intersections have negative t" {
    const alctr = std.testing.allocator;
    const int1 = Intersection.init(-1.0, VolPtr{ .sphere_idx = 5 });
    const int2 = Intersection.init(-1.0, VolPtr{ .sphere_idx = 5 });
    const xs = Intersections.init(alctr, .{ int1, int2 });
    defer xs.deinit();

    try expect(xs.hit() == null);
}

test "The hit is always the lowest nonnegative intersection" {
    const alctr = std.testing.allocator;
    const int1 = Intersection.init(5.0, VolPtr{ .sphere_idx = 5 });
    const int2 = Intersection.init(7.0, VolPtr{ .sphere_idx = 5 });
    const int3 = Intersection.init(-3.0, VolPtr{ .sphere_idx = 5 });
    const int4 = Intersection.init(2.0, VolPtr{ .sphere_idx = 5 });
    const xs = Intersections.init(alctr, .{ int1, int2, int3, int4 });
    defer xs.deinit();

    try expect(xs.hit().?.equals(int4));
}
