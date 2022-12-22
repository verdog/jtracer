//! translation, rotation, scaling operations on matrices

const std = @import("std");

const Matrix = @import("matrix.zig").Matrix;
const Mat4 = Matrix(4, 4);
const Mat3 = Matrix(3, 3);

const Tuple = @import("tuple.zig").Tuple;
const Point = @import("tuple.zig").Point;
const Vector = @import("tuple.zig").Vector;

const pi = std.math.pi;

pub const Transform = struct {
    t: Matrix(4, 4) = Matrix(4, 4).identity(),
    inverse: Matrix(4, 4) = Matrix(4, 4).identity(),

    pub fn init(inits: std.meta.Tuple(&[_]type{f64} ** 16)) This {
        const t = Mat4.init(inits);
        const i = t.inverted() catch unreachable;
        return This{ .t = t, .inverse = i };
    }

    pub fn mult(self: This, other: This) This {
        const t = self.t.mult(other.t);
        const i = t.inverted() catch unreachable;
        return This{ .t = t, .inverse = i };
    }

    pub fn chain(self: This, others: anytype) This {
        var result = self.t;

        inline for (others) |tr| {
            result = result.mult(tr.t);
        }

        const iresult = result.inverted() catch unreachable;

        return This{ .t = result, .inverse = iresult };
    }

    const This = @This();
};

pub fn makeTranslation(x: f64, y: f64, z: f64) Transform {
    return Transform.init(.{
        1, 0, 0, x,
        0, 1, 0, y,
        0, 0, 1, z,
        0, 0, 0, 1,
    });
}

pub fn makeScaling(x: f64, y: f64, z: f64) Transform {
    return Transform.init(.{
        x, 0, 0, 0,
        0, y, 0, 0,
        0, 0, z, 0,
        0, 0, 0, 1,
    });
}

pub fn makeRotationX(rads: f64) Transform {
    return Transform.init(.{
        1, 0,          0,           0,
        0, @cos(rads), -@sin(rads), 0,
        0, @sin(rads), @cos(rads),  0,
        0, 0,          0,           1,
    });
}

pub fn makeRotationY(rads: f64) Transform {
    return Transform.init(.{
        @cos(rads),  0, @sin(rads), 0,
        0,           1, 0,          0,
        -@sin(rads), 0, @cos(rads), 0,
        0,           0, 0,          1,
    });
}

pub fn makeRotationZ(rads: f64) Transform {
    return Transform.init(.{
        @cos(rads), -@sin(rads), 0, 0,
        @sin(rads), @cos(rads),  0, 0,
        0,          0,           1, 0,
        0,          0,           0, 1,
    });
}

pub fn makeShearing(xy: f64, xz: f64, yx: f64, yz: f64, zx: f64, zy: f64) Transform {
    return Transform.init(.{
        1,  xy, xz, 0,
        yx, 1,  yz, 0,
        zx, zy, 1,  0,
        0,  0,  0,  1,
    });
}

pub fn makeView(from: Tuple, to: Tuple, up: Tuple) Transform {
    // page 99-100
    std.debug.assert(from.isPoint());
    std.debug.assert(to.isPoint());
    std.debug.assert(up.isVector());

    const forward = to.minus(from).normalized();
    const left = forward.cross(up.normalized());
    const true_up = left.cross(forward);

    const orientation = Mat4.init(.{
        left.x(),     left.y(),     left.z(),     0,
        true_up.x(),  true_up.y(),  true_up.z(),  0,
        -forward.x(), -forward.y(), -forward.z(), 0,
        0,            0,            0,            1,
    });

    const t = orientation.mult(makeTranslation(-from.x(), -from.y(), -from.z()).t);

    return Transform{ .t = t, .inverse = t.inverted() catch unreachable };
}

const expect = std.testing.expect;
const print = @import("u.zig").print;

test "Multiplying by a translation matrix" {
    const t = makeTranslation(5, -3, 2);
    const p = Point.init(-3, 4, 5);
    const td = t.t.mult(p);

    try expect(td.equals(Point.init(2, 1, 7)));
}

test "Multiplying by the inverse of a translation matrix" {
    const t = makeTranslation(5, -3, 2);
    const p = Point.init(-3, 4, 5);
    const td = t.inverse.mult(p);

    try expect(td.equals(Point.init(-8, 7, 3)));
}

test "Translation does not affect vectors" {
    const t = makeTranslation(5, -3, 2);
    const v = Vector.init(-3, 4, 5);
    const td = t.t.mult(v);

    try expect(td.equals(v));
}

test "A scaling matrix applied to a point" {
    const t = makeScaling(2, 3, 4);
    const p = Point.init(-4, 6, 8);
    const td = t.t.mult(p);

    try expect(td.equals(Point.init(-8, 18, 32)));
}

test "A scaling matrix applied to a vector" {
    const t = makeScaling(2, 3, 4);
    const v = Vector.init(-4, 6, 8);
    const td = t.t.mult(v);

    try expect(td.equals(Vector.init(-8, 18, 32)));
}

test "Multiplying by the inverse of a scaling matrix" {
    const t = makeScaling(2, 3, 4);
    const v = Vector.init(-4, 6, 8);
    const td = t.inverse.mult(v);

    try expect(td.equals(Vector.init(-2, 2, 2)));
}

test "Reflection is scaling by a negative value" {
    const t = makeScaling(-1, 1, 1);
    const v = Point.init(-4, 6, 8);
    const td = t.t.mult(v);

    try expect(td.equals(Point.init(4, 6, 8)));
}

test "Rotating a point around the x axis" {
    const p = Point.init(0, 1, 0);
    const half_quarter = makeRotationX(pi / 4.0);
    const full_quarter = makeRotationX(pi / 2.0);

    try expect(half_quarter.t.mult(p).equals(Point.init(0, @sqrt(2.0) / 2.0, @sqrt(2.0) / 2.0)));
    try expect(full_quarter.t.mult(p).equals(Point.init(0, 0, 1)));
}

test "Rotation a point around the x asix; inverted" {
    const p = Point.init(0, 1, 0);
    const half_quarter = makeRotationX(pi / 4.0);

    const hqid = half_quarter.inverse.mult(p);
    const eq = Point.init(0, @sqrt(2.0) / 2.0, -@sqrt(2.0) / 2.0);

    errdefer std.debug.print("{} != {}\n", .{ hqid.vec, eq.vec });

    try expect(hqid.equals(eq));
}

test "Rotating a point around the y axis" {
    const p = Point.init(0, 0, 1);
    const half_quarter = makeRotationY(pi / 4.0);
    const full_quarter = makeRotationY(pi / 2.0);

    try expect(half_quarter.t.mult(p).equals(Point.init(@sqrt(2.0) / 2.0, 0, @sqrt(2.0) / 2.0)));
    try expect(full_quarter.t.mult(p).equals(Point.init(1, 0, 0)));
}

test "Rotation a point around the y asix; inverted" {
    const p = Point.init(0, 0, 1);
    const half_quarter = makeRotationY(pi / 4.0);

    const hqid = half_quarter.inverse.mult(p);
    const eq = Point.init(-@sqrt(2.0) / 2.0, 0, @sqrt(2.0) / 2.0);

    errdefer std.debug.print("{} != {}\n", .{ hqid.vec, eq.vec });

    try expect(hqid.equals(eq));
}

test "Rotating a point around the z axis" {
    const p = Point.init(0, 1, 0);
    const half_quarter = makeRotationZ(pi / 4.0);
    const full_quarter = makeRotationZ(pi / 2.0);

    try expect(half_quarter.t.mult(p).equals(Point.init(-@sqrt(2.0) / 2.0, @sqrt(2.0) / 2.0, 0)));
    try expect(full_quarter.t.mult(p).equals(Point.init(-1, 0, 0)));
}

test "Rotation a point around the z asix; inverted" {
    const p = Point.init(0, 1, 0);
    const half_quarter = makeRotationZ(pi / 4.0);

    const hqid = half_quarter.inverse.mult(p);
    const eq = Point.init(@sqrt(2.0) / 2.0, @sqrt(2.0) / 2.0, 0);

    errdefer std.debug.print("{} != {}\n", .{ hqid.vec, eq.vec });

    try expect(hqid.equals(eq));
}

test "A shearing transformation moves x in proportion to y" {
    const p = Point.init(2, 3, 4);
    const t = makeShearing(1, 0, 0, 0, 0, 0);
    const td = t.t.mult(p);

    try expect(td.equals(Point.init(5, 3, 4)));
}

test "A shearing transformation moves x in proportion to z" {
    const p = Point.init(2, 3, 4);
    const t = makeShearing(0, 1, 0, 0, 0, 0);
    const td = t.t.mult(p);

    try expect(td.equals(Point.init(6, 3, 4)));
}

test "A shearing transformation moves y in proportion to x" {
    const p = Point.init(2, 3, 4);
    const t = makeShearing(0, 0, 1, 0, 0, 0);
    const td = t.t.mult(p);

    try expect(td.equals(Point.init(2, 5, 4)));
}

test "A shearing transformation moves y in proportion to z" {
    const p = Point.init(2, 3, 4);
    const t = makeShearing(0, 0, 0, 1, 0, 0);
    const td = t.t.mult(p);

    try expect(td.equals(Point.init(2, 7, 4)));
}

test "A shearing transformation moves z in proportion to x" {
    const p = Point.init(2, 3, 4);
    const t = makeShearing(0, 0, 0, 0, 1, 0);
    const td = t.t.mult(p);

    try expect(td.equals(Point.init(2, 3, 6)));
}

test "A shearing transformation moves z in proportion to y" {
    const p = Point.init(2, 3, 4);
    const t = makeShearing(0, 0, 0, 0, 0, 1);
    const td = t.t.mult(p);

    try expect(td.equals(Point.init(2, 3, 7)));
}

test "Individual transformations are applied in sequence" {
    var p = Point.init(1, 0, 1);
    const rot = makeRotationX(pi / 2.0);
    const scl = makeScaling(5, 5, 5);
    const trn = makeTranslation(10, 5, 7);

    p = rot.t.mult(p);
    try expect(p.equals(Point.init(1, -1, 0)));

    p = scl.t.mult(p);
    try expect(p.equals(Point.init(5, -5, 0)));

    p = trn.t.mult(p);
    try expect(p.equals(Point.init(15, 0, 7)));
}

test "Chained transformations are applied in reverse order" {
    var p = Point.init(1, 0, 1);

    const rot = makeRotationX(pi / 2.0);
    const scl = makeScaling(5, 5, 5);
    const trn = makeTranslation(10, 5, 7);

    const final = trn.t.mult(scl.t).mult(rot.t);

    p = final.mult(p);

    try expect(p.equals(Point.init(15, 0, 7)));
}

test "Default view transform" {
    const from = Point.init(0, 0, 0);
    const to = Point.init(0, 0, -1);
    const up = Vector.init(0, 1, 0);

    const t = makeView(from, to, up);

    try expect(t.t.equals(Mat4.identity()));
}

test "View transform looking in positive Z" {
    const from = Point.init(0, 0, 0);
    const to = Point.init(0, 0, 1);
    const up = Vector.init(0, 1, 0);

    const t = makeView(from, to, up);

    try expect(t.t.equals(makeScaling(-1, 1, -1).t));
}

test "The view transform moves the world, not the eye" {
    const from = Point.init(0, 0, 8);
    const to = Point.init(0, 0, 0);
    const up = Vector.init(0, 1, 0);

    const t = makeView(from, to, up);

    errdefer print(t);
    try expect(t.t.equals(makeTranslation(0, 0, -8).t));
}

test "An arbitrary view transform" {
    const from = Point.init(1, 3, 2);
    const to = Point.init(4, -2, 8);
    const up = Vector.init(1, 1, 0);

    const t = makeView(from, to, up);

    // TODO book tests are imprecise
    try expect(t.t.equalsTolerance(Mat4.init(.{
        -0.50709, 0.50709, 0.67612,  -2.36643,
        0.76772,  0.60609, 0.12122,  -2.82843,
        -0.35857, 0.59761, -0.71714, 0.0,
        0.0,      0.0,     0.0,      1.0,
    }), 100_000_000_000));
}

test "Transform: default initialization" {
    const t = Transform{};

    try expect(t.t.equals(Matrix(4, 4).identity()));
    try expect(t.inverse.equals(Matrix(4, 4).identity()));
}

test "Transform: mult" {
    const p = Point.init(0, 0, 0);
    const t = Transform{};
    const t2 = t.mult(makeTranslation(2, 4, 8));

    try expect(std.meta.eql(t2, makeTranslation(2, 4, 8)));
    try expect(t2.t.mult(p).equals(Point.init(2, 4, 8)));
}

test "Transform: chained mults are applied in reverse" {
    const p = Point.init(1, 0, 0);
    const t = Transform{};
    const t2 = t.mult(makeTranslation(2, 4, 0)).mult(makeScaling(2, 2, 2));

    try expect(t2.t.mult(p).equals(Point.init(4, 4, 0)));
}

test "Transform: chain() applied in reverse order" {
    const p = Point.init(1, 0, 0);
    const t = Transform{};
    const t2 = t.chain(.{
        makeTranslation(2, 4, 0),
        makeScaling(2, 2, 2),
    });

    try expect(t2.t.mult(p).equals(Point.init(4, 4, 0)));
}
