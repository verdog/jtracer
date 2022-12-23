const std = @import("std");

const Tuple = @import("tuple.zig").Tuple;
const Vector = @import("tuple.zig").Vector;
const Point = @import("tuple.zig").Point;

const mtx = @import("matrix.zig");
const trans = @import("transform.zig");
const mymath = @import("mymath.zig");

pub const Ray = struct {
    origin: Tuple,
    direction: Tuple,

    pub fn init(origin: Tuple, direction: Tuple) This {
        std.debug.assert(origin.isPoint());
        std.debug.assert(direction.isVector());
        std.debug.assert(std.math.approxEqRel(f64, direction.magnitude(), 1.0, mymath.floatTolerance));

        return .{
            .origin = origin,
            .direction = direction,
        };
    }

    pub fn position(self: This, t: f64) Tuple {
        return self.origin.plus(self.direction.scaled(t));
    }

    pub fn transformed(self: This, t: mtx.Matrix(4, 4)) This {
        return .{
            .origin = t.mult(self.origin),
            .direction = t.mult(self.direction),
        };
    }

    const This = @This();
};

const expect = std.testing.expect;

test "Creating and querying a ray" {
    const origin = Point.init(1, 2, 3);
    const direction = Vector.init(4, 5, 6).normalized();

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

test "Translating a ray" {
    const ray = Ray.init(Point.init(1, 2, 3), Vector.init(0, 1, 0));
    const m = trans.makeTranslation(3, 4, 5);
    const ray2 = ray.transformed(m.t);

    try expect(ray2.origin.equals(Point.init(4, 6, 8)));
    try expect(ray2.direction.equals(Vector.init(0, 1, 0)));
}

test "Scaling a ray" {
    const ray = Ray.init(Point.init(1, 2, 3), Vector.init(0, 1, 0));
    const m = trans.makeScaling(2, 3, 4);
    const ray2 = ray.transformed(m.t);

    try expect(ray2.origin.equals(Point.init(2, 6, 12)));
    try expect(ray2.direction.equals(Vector.init(0, 3, 0)));
}
