//! intersection algorithms

/// get information that will be needed for shading based on an intersection
/// and the ray that generated it
pub const HitData = struct {
    /// the intersection that the HitData is based on
    intersection: Intersection,
    /// the point of the hit in world space
    point: Tuple,
    /// point slightly offset by normal_vector
    over_point: Tuple,
    /// point slightly offset by -normal_vector
    under_point: Tuple,
    /// the eye vector, going from point to the origin point of the
    /// ray that generated the intersection
    eye_vector: Tuple,
    /// the normal vector, going from the point outwards from the surface
    /// that was hit
    normal_vector: Tuple,
    /// the incoming ray bouncing off of over_point, reflected across normal_vector
    reflect_vector: Tuple,
    /// true if the normal vector points towards the inside of the shape
    inside: bool,
    /// index of refraction of the material that the hit is exiting n1: f64,
    n1: f64,
    /// index of refraction of the material that the hit is entering
    n2: f64,

    pub fn init(r: Ray, x: Intersection, normal: Tuple, n1: f64, n2: f64) This {
        const point = r.position(x.t);
        var normal_vector = normal;
        const eye_vector = r.direction.scaled(-1).normalized();
        const inside = eye_vector.dot(normal_vector) < 0;
        if (inside) normal_vector = normal_vector.scaled(-1);
        const over_point = point.plus(normal_vector.scaled(32 * mymath.floatTolerance));
        const under_point = point.plus(normal_vector.scaled(-32 * mymath.floatTolerance));
        const reflect_vector = r.direction.normalized().reflected(normal_vector);

        std.debug.assert(std.math.approxEqRel(
            f64,
            normal_vector.magnitude(),
            1.0,
            mymath.floatTolerance,
        ));
        std.debug.assert(std.math.approxEqRel(
            f64,
            eye_vector.magnitude(),
            1.0,
            mymath.floatTolerance,
        ));
        std.debug.assert(std.math.approxEqRel(
            f64,
            reflect_vector.magnitude(),
            1.0,
            mymath.floatTolerance,
        ));

        return .{
            .intersection = x,
            .point = point,
            .over_point = over_point,
            .under_point = under_point,
            .eye_vector = eye_vector,
            .normal_vector = normal_vector,
            .reflect_vector = reflect_vector,
            .inside = inside,
            .n1 = n1,
            .n2 = n2,
        };
    }

    pub fn schlick(self: This) f64 {
        // schlick approximation to determine total internal reflection

        // find the cosine of the angle btween the eye and normal
        var cos = self.eye_vector.dot(self.normal_vector);

        // total internal reflection can only occur if n1 > n2
        if (self.n1 > self.n2) {
            // TODO is infinity ok here?
            const ratio = self.n1 / self.n2;
            const sin2_t = ratio * ratio * (1.0 - cos * cos);
            if (sin2_t > 1.0) return 1;

            // compute cosine of theta_t using trig identity.
            // when n1 > n2, use cos(theta_t) instead of cos(eye dot normal)
            cos = @sqrt(1.0 - sin2_t);
        }

        const inner = (self.n1 - self.n2) / (self.n1 + self.n2);
        const r0 = inner * inner;
        return r0 + (1 - r0) * std.math.pow(f64, 1 - cos, 5);
    }

    const This = @This();
};

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
    containers: std.ArrayList(Intersection),

    pub fn init(alctr: std.mem.Allocator, inits: anytype) This {
        const InitsT = @TypeOf(inits);
        const inits_type_info = @typeInfo(InitsT);
        if (inits_type_info != .Struct) {
            @compileError("Init with .{} syntax");
        }
        const fields_info = inits_type_info.Struct.fields;

        var xs = This{
            .ixs = std.ArrayList(Intersection).init(alctr),
            .containers = std.ArrayList(Intersection).init(alctr),
        };

        inline for (fields_info) |field| {
            xs.ixs.append(@field(inits, field.name)) catch unreachable;
        }

        return xs;
    }

    pub fn deinit(self: This) void {
        self.ixs.deinit();
        self.containers.deinit();
    }

    pub fn intersect(self: *This, volu: anytype, vptr: VolPtr, ray: Ray) void {
        switch (@TypeOf(volu)) {
            vol.Sphere => {
                std.debug.assert(std.meta.activeTag(vptr) == .sphere_idx);
                self.intersectSphere(vptr, volu, ray);
            },
            vol.Plane => {
                std.debug.assert(std.meta.activeTag(vptr) == .plane_idx);
                self.intersectPlane(vptr, volu, ray);
            },
            else => unreachable,
        }
    }

    fn intersectSphere(self: *This, vptr: VolPtr, sphere: vol.Sphere, ray: Ray) void {
        const td_ray = ray.transformed(sphere.transform.inverse);
        const sphere_to_ray = td_ray.origin.minus(Point.init(0, 0, 0));

        const a = td_ray.direction.dot(td_ray.direction);
        const b = 2 * td_ray.direction.dot(sphere_to_ray);
        const c = sphere_to_ray.dot(sphere_to_ray) - 1;

        const discriminant = b * b - 4 * a * c;

        if (discriminant < 0) return;

        self.ixs.append(Intersection.init((-b - @sqrt(discriminant)) / (2 * a), vptr)) catch @panic("OOM");
        self.ixs.append(Intersection.init((-b + @sqrt(discriminant)) / (2 * a), vptr)) catch @panic("OOM");
    }

    fn intersectPlane(self: *This, vptr: VolPtr, plane: vol.Plane, ray: Ray) void {
        const td_ray = ray.transformed(plane.transform.inverse);

        if (@fabs(td_ray.direction.y()) < mymath.floatTolerance) {
            // ray is parallel (enough) to plane, no hit
            return;
        }

        self.ixs.append(Intersection.init(-td_ray.origin.y() / td_ray.direction.y(), vptr)) catch @panic("OOM");
    }

    pub fn clear(self: *This) void {
        self.ixs.clearRetainingCapacity();
    }

    /// sort intersections by increasing t value
    pub fn order(self: *This) void {
        const lt = struct {
            fn lt(_: void, lhs: Intersection, rhs: Intersection) bool {
                return lhs.t < rhs.t;
            }
        }.lt;

        std.sort.insertionSort(Intersection, self.ixs.items, {}, lt);
    }

    pub fn hit(self: *This) ?Intersection {
        // find smallest positive t
        var result: ?Intersection = null;

        for (self.ixs.items) |x| {
            if (x.t > 0 and (result == null or x.t < result.?.t)) {
                result = x;
            }
        }

        return result;
    }

    pub const Boundary = struct {
        // lesser has a smaller t value than greater
        lesser: ?VolPtr,
        greater: ?VolPtr,
    };

    /// Given an intersection, search internal list of intersections and return `exiting`
    /// and `entering`, where `exiting` is a VolPtr to the object that the hit is escaping,
    /// and `entering` is a VolPtr to the object that the hit is entering. Both fields of
    /// the return struct are optional, since a hit might be exiting into nothing or vice versa.
    pub fn findBoundaryObjects(self: *This, x: Intersection) Boundary {
        self.order();
        self.containers.clearRetainingCapacity();
        var result = Boundary{
            .lesser = null,
            .greater = null,
        };

        for (self.ixs.items) |my_x| {
            if (std.meta.eql(x, my_x)) {
                if (self.containers.items.len == 0) {
                    result.lesser = null;
                } else {
                    result.lesser = self.containers.items[self.containers.items.len - 1].vptr;
                }
            }

            const my_x_idx: ?usize = blk: {
                for (self.containers.items) |c, i| {
                    if (std.meta.eql(my_x.vptr, c.vptr)) {
                        break :blk i;
                    }
                }
                break :blk null;
            };

            if (my_x_idx) |idx| {
                _ = self.containers.orderedRemove(idx);
            } else {
                self.containers.append(my_x) catch unreachable;
            }

            if (std.meta.eql(x, my_x)) {
                if (self.containers.items.len == 0) {
                    result.greater = null;
                } else {
                    result.greater = self.containers.items[self.containers.items.len - 1].vptr;
                }
                break;
            }
        }

        return result;
    }

    const This = @This();
};

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
    const sphere = vol.Sphere.init();
    const vptr = VolPtr{ .sphere_idx = 0 };

    var xs = Intersections.init(alctr, .{});
    defer xs.deinit();
    xs.intersect(sphere, vptr, ray);

    try expect(xs.ixs.items.len == 2);
    try expect(xs.ixs.items[0].t == 4.0);
    try expect(xs.ixs.items[1].t == 6.0);
}

test "A ray intersects a sphere at a tangent" {
    const alctr = std.testing.allocator;
    const ray = Ray.init(Point.init(0, 1, -5), Vector.init(0, 0, 1));
    const sphere = vol.Sphere.init();
    const vptr = VolPtr{ .sphere_idx = 0 };

    var xs = Intersections.init(alctr, .{});
    defer xs.deinit();
    xs.intersect(sphere, vptr, ray);

    try expect(xs.ixs.items.len == 2);
    try expect(xs.ixs.items[0].t == 5.0);
    try expect(xs.ixs.items[1].t == 5.0);
}

test "A ray misses a sphere" {
    const alctr = std.testing.allocator;
    const ray = Ray.init(Point.init(0, 2, -5), Vector.init(0, 0, 1));
    const sphere = vol.Sphere.init();
    const vptr = VolPtr{ .sphere_idx = 0 };

    var xs = Intersections.init(alctr, .{});
    defer xs.deinit();
    xs.intersect(sphere, vptr, ray);

    try expect(xs.ixs.items.len == 0);
}

test "A ray originates inside a sphere" {
    const alctr = std.testing.allocator;
    const ray = Ray.init(Point.init(0, 0, 0), Vector.init(0, 0, 1));
    const sphere = vol.Sphere.init();
    const vptr = VolPtr{ .sphere_idx = 0 };

    var xs = Intersections.init(alctr, .{});
    defer xs.deinit();
    xs.intersect(sphere, vptr, ray);

    try expect(xs.ixs.items.len == 2);
    try expect(xs.ixs.items[0].t == -1.0);
    try expect(xs.ixs.items[1].t == 1.0);
}

test "A sphere behind a ray" {
    const alctr = std.testing.allocator;
    const ray = Ray.init(Point.init(0, 0, 5), Vector.init(0, 0, 1));
    const sphere = vol.Sphere.init();
    const vptr = VolPtr{ .sphere_idx = 0 };

    var xs = Intersections.init(alctr, .{});
    defer xs.deinit();
    xs.intersect(sphere, vptr, ray);

    try expect(xs.ixs.items.len == 2);
    try expect(xs.ixs.items[0].t == -6.0);
    try expect(xs.ixs.items[1].t == -4.0);
}

test "Intersect a scaled sphere with a ray" {
    const alctr = std.testing.allocator;
    const ray = Ray.init(Point.init(0, 0, -5), Vector.init(0, 0, 1));
    var s = vol.Sphere.init();
    const vptr = VolPtr{ .sphere_idx = 7 };

    s.transform = trans.makeScaling(2, 2, 2);

    var xs = Intersections.init(alctr, .{});
    defer xs.deinit();
    xs.intersect(s, vptr, ray);

    try expect(xs.ixs.items.len == 2);
    try expect(xs.ixs.items[0].t == 3.0);
    try expect(xs.ixs.items[1].t == 7.0);
}

test "Intersect a translated sphere with a ray" {
    const alctr = std.testing.allocator;
    const ray = Ray.init(Point.init(0, 0, -5), Vector.init(0, 0, 1));
    var s = vol.Sphere.init();
    const vptr = VolPtr{ .sphere_idx = 7 };

    s.transform = trans.makeTranslation(5, 0, 0);

    var xs = Intersections.init(alctr, .{});
    defer xs.deinit();
    xs.intersect(s, vptr, ray);

    try expect(xs.ixs.items.len == 0);
}

test "Intersect sets the object on the intersection" {
    const alctr = std.testing.allocator;
    const ray = Ray.init(Point.init(0, 0, -5), Vector.init(0, 0, 1));
    const sph = vol.Sphere.init();
    const vptr = VolPtr{ .sphere_idx = 7 };

    var xs = Intersections.init(alctr, .{});
    defer xs.deinit();
    xs.intersect(sph, vptr, ray);

    try expect(xs.ixs.items.len == 2);
    try expect(xs.ixs.items[0].vptr.sphere_idx == 7);
    try expect(xs.ixs.items[1].vptr.sphere_idx == 7);
}

test "The hit, where all intersections have positive t" {
    const alctr = std.testing.allocator;
    const int1 = Intersection.init(1.0, VolPtr{ .sphere_idx = 5 });
    const int2 = Intersection.init(2.0, VolPtr{ .sphere_idx = 5 });

    var xs = Intersections.init(alctr, .{ int1, int2 });
    defer xs.deinit();

    try expect(xs.hit().?.equals(int1));
}

test "The hit, where some intersections have negative t" {
    const alctr = std.testing.allocator;
    const int1 = Intersection.init(-1.0, VolPtr{ .sphere_idx = 5 });
    const int2 = Intersection.init(1.0, VolPtr{ .sphere_idx = 5 });
    var xs = Intersections.init(alctr, .{ int1, int2 });
    defer xs.deinit();

    try expect(xs.hit().?.equals(int2));
}

test "The hit, where all intersections have negative t" {
    const alctr = std.testing.allocator;
    const int1 = Intersection.init(-1.0, VolPtr{ .sphere_idx = 5 });
    const int2 = Intersection.init(-1.0, VolPtr{ .sphere_idx = 5 });
    var xs = Intersections.init(alctr, .{ int1, int2 });
    defer xs.deinit();

    try expect(xs.hit() == null);
}

test "The hit is always the lowest nonnegative intersection" {
    const alctr = std.testing.allocator;
    const int1 = Intersection.init(5.0, VolPtr{ .sphere_idx = 5 });
    const int2 = Intersection.init(7.0, VolPtr{ .sphere_idx = 5 });
    const int3 = Intersection.init(-3.0, VolPtr{ .sphere_idx = 5 });
    const int4 = Intersection.init(2.0, VolPtr{ .sphere_idx = 5 });
    var xs = Intersections.init(alctr, .{ int1, int2, int3, int4 });
    defer xs.deinit();

    try expect(xs.hit().?.equals(int4));
}

test "HitData" {
    const r = Ray.init(Point.init(0, 0, -5), Vector.init(0, 0, 1));
    const s = vol.Sphere.init();
    const x = Intersection.init(4.0, VolPtr{ .sphere_idx = 0 });
    const data = HitData.init(r, x, s.normalAt(r.position(x.t)), 1, 1);

    try expect(data.intersection.t == x.t);
    try expect(std.meta.eql(data.intersection.vptr, VolPtr{ .sphere_idx = 0 }));
    try expect(data.point.equals(Point.init(0, 0, -1)));
    try expect(data.eye_vector.equals(Vector.init(0, 0, -1)));
    try expect(data.normal_vector.equals(Vector.init(0, 0, -1)));
}

test "HitData: eye vector outside of the hit shape" {
    const r = Ray.init(Point.init(0, 0, -5), Vector.init(0, 0, 1));
    const s = vol.Sphere.init();
    const x = Intersection.init(4.0, VolPtr{ .sphere_idx = 0 });
    const data = HitData.init(r, x, s.normalAt(r.position(x.t)), 1, 1);

    try expect(data.inside == false);
}

test "HitData: eye vector inside of the hit shape" {
    const r = Ray.init(Point.init(0, 0, 0), Vector.init(0, 0, 1));
    const s = vol.Sphere.init();
    const x = Intersection.init(1.0, VolPtr{ .sphere_idx = 0 });
    const data = HitData.init(r, x, s.normalAt(r.position(x.t)), 1, 1);

    try expect(data.point.equals(Point.init(0, 0, 1)));
    try expect(data.inside == true);
    try expect(data.eye_vector.equals(Vector.init(0, 0, -1)));
    // normal would have been (0, 0, 1), but it is inverted since inside == true
    try expect(data.normal_vector.equals(Vector.init(0, 0, -1)));
}

test "HitData: the hit should offset the point" {
    const r = Ray.init(Point.init(0, 0, -5), Vector.init(0, 0, 1));
    var s = vol.Sphere.init();
    s.transform = trans.makeTranslation(0, 0, 1);
    const x = Intersection.init(5.0, VolPtr{ .sphere_idx = 0 });
    const data = HitData.init(r, x, s.normalAt(r.position(x.t)), 1, 1);

    try expect(data.over_point.z() < -mymath.floatTolerance * 16);
    try expect(data.point.z() > data.over_point.z());
}

test "HitData: The under point is offset below the surface" {
    const r = Ray.initC(0, 0, -5, 0, 0, 1);
    var s = vol.Sphere.init();
    s.transform = trans.makeTranslation(0, 0, 1);
    const sptr = VolPtr{ .sphere_idx = 0 };
    const x = Intersection.init(5, sptr);

    const data = HitData.init(r, x, s.normalAt(r.position(x.t)), 1, 1);

    try expect(data.under_point.z() > mymath.floatTolerance * 16);
    try expect(data.point.z() < data.under_point.z());
}

test "Intersect with a ray parallel to the plane should yield nothing" {
    const p = vol.Plane.init();
    const pptr = VolPtr{ .plane_idx = 0 };
    const r = Ray.init(Point.init(0, 1, 0), Vector.init(1, 0, 0));

    var xs = Intersections.init(std.testing.allocator, .{});
    defer xs.deinit();

    xs.intersect(p, pptr, r);

    try expect(xs.ixs.items.len == 0);
}

test "Intersect with a coplanar ray should yield nothing" {
    const p = vol.Plane.init();
    const pptr = VolPtr{ .plane_idx = 0 };
    const r = Ray.init(Point.init(0, 0, 0), Vector.init(1, 0, 0));

    var xs = Intersections.init(std.testing.allocator, .{});
    defer xs.deinit();

    xs.intersect(p, pptr, r);

    try expect(xs.ixs.items.len == 0);
}

test "A ray intersecting a plane from above" {
    const p = vol.Plane.init();
    const pptr = VolPtr{ .plane_idx = 0 };
    const r = Ray.init(Point.init(0, 10, 0), Vector.init(0, -1, 0));

    var xs = Intersections.init(std.testing.allocator, .{});
    defer xs.deinit();

    xs.intersect(p, pptr, r);

    try expect(xs.ixs.items.len == 1);
    try expect(xs.ixs.items[0].t == 10);
    try expect(std.meta.eql(xs.ixs.items[0].vptr, pptr));
}

test "A ray intersecting a plane from below" {
    const p = vol.Plane.init();
    const pptr = VolPtr{ .plane_idx = 0 };
    const r = Ray.init(Point.init(0, -3, 0), Vector.init(0, 1, 0));

    var xs = Intersections.init(std.testing.allocator, .{});
    defer xs.deinit();

    xs.intersect(p, pptr, r);

    try expect(xs.ixs.items.len == 1);
    try expect(xs.ixs.items[0].t == 3);
    try expect(std.meta.eql(xs.ixs.items[0].vptr, pptr));
}

test "A ray intersecting a plane from an angle" {
    const p = vol.Plane.init();
    const pptr = VolPtr{ .plane_idx = 0 };
    const r = Ray.init(Point.init(-3, 4, 0), Vector.init(3, -4, 0).normalized());

    var xs = Intersections.init(std.testing.allocator, .{});
    defer xs.deinit();

    xs.intersect(p, pptr, r);

    errdefer print(xs.ixs.items);

    try expect(xs.ixs.items.len == 1);
    try expect(xs.ixs.items[0].t == 5);
    try expect(std.meta.eql(xs.ixs.items[0].vptr, pptr));
}

test "Precomputing the reflection vector" {
    const pl = vol.Plane.init();
    const r = Ray.init(Point.init(0, 1, -1), Vector.init(0, -@sqrt(2.0) / 2.0, @sqrt(2.0) / 2.0));
    const i = Intersection.init(@sqrt(2.0), VolPtr{ .plane_idx = 0 });

    const data = HitData.init(r, i, pl.normalAt(r.position(i.t)), 1, 1);

    try expect(data.reflect_vector.equals(Vector.init(0, @sqrt(2.0) / 2.0, @sqrt(2.0) / 2.0)));
}

test "The Schlick approx. under total internal reflection" {
    const alctr = std.testing.allocator;

    var sh = vol.Sphere.init();
    prefab.toGlass(&sh.material);
    const r = Ray.initC(0, 0, @sqrt(2.0) / 2.0, 0, 1, 0);
    const xs = Intersections.init(alctr, .{
        Intersection.init(-@sqrt(2.0) / 2.0, VolPtr{ .sphere_idx = 0 }),
        Intersection.init(@sqrt(2.0) / 2.0, VolPtr{ .sphere_idx = 0 }),
    });
    defer xs.deinit();
    const data = HitData.init(r, xs.ixs.items[1], sh.normalAt(r.position(xs.ixs.items[1].t)), 1.5, 1);

    try expect(data.schlick() == 1.0);
}

test "The Schlick approx. with a perpendicular viewing angle" {
    const alctr = std.testing.allocator;

    var sh = vol.Sphere.init();
    prefab.toGlass(&sh.material);
    const r = Ray.initC(0, 0, 0, 0, 1, 0);
    const xs = Intersections.init(alctr, .{
        Intersection.init(-1, VolPtr{ .sphere_idx = 0 }),
        Intersection.init(1, VolPtr{ .sphere_idx = 0 }),
    });
    defer xs.deinit();
    const data = HitData.init(r, xs.ixs.items[1], sh.normalAt(r.position(xs.ixs.items[1].t)), 1.5, 1);

    try expect(std.math.approxEqAbs(f64, data.schlick(), 0.04, mymath.floatTolerance));
}

test "The Schlick approx. with a small angle and n2 > n1" {
    const alctr = std.testing.allocator;

    var sh = vol.Sphere.init();
    prefab.toGlass(&sh.material);
    const r = Ray.initC(0, 0.99, -2, 0, 0, 1);
    const xs = Intersections.init(alctr, .{
        Intersection.init(1.8589, VolPtr{ .sphere_idx = 0 }),
    });
    defer xs.deinit();
    const data = HitData.init(r, xs.ixs.items[0], sh.normalAt(r.position(xs.ixs.items[0].t)), 1, 1.5);

    errdefer std.debug.print("{}\n", .{data.schlick()});
    try expect(std.math.approxEqRel(f64, data.schlick(), 0.48873081012212, mymath.floatTolerance));
}

const std = @import("std");

const Tuple = @import("tuple.zig").Tuple;
const Vector = @import("tuple.zig").Vector;
const Point = @import("tuple.zig").Point;
const Ray = @import("ray.zig").Ray;

const VolPtr = @import("world.zig").VolumePtr;

const trans = @import("transform.zig");
const mtx = @import("matrix.zig");
const mymath = @import("mymath.zig");
const vol = @import("volume.zig");
const prefab = @import("prefab.zig");

const expect = std.testing.expect;
const print = @import("u.zig").print;
