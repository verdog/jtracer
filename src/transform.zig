//! translation, rotation, scaling operations on matrices

const std = @import("std");

const Matrix = @import("matrix.zig").Matrix;
const Mat4 = Matrix(4, 4);
const Mat3 = Matrix(3, 3);

const Point = @import("tuple.zig").Point;
const Vector = @import("tuple.zig").Vector;

const pi = std.math.pi;

pub fn makeTranslation(x: f64, y: f64, z: f64) Mat4 {
    return Mat4.init(.{
        1, 0, 0, x,
        0, 1, 0, y,
        0, 0, 1, z,
        0, 0, 0, 1,
    });
}

pub fn makeScaling(x: f64, y: f64, z: f64) Mat4 {
    return Mat4.init(.{
        x, 0, 0, 0,
        0, y, 0, 0,
        0, 0, z, 0,
        0, 0, 0, 1,
    });
}

pub fn makeRotationX(rads: f64) Mat4 {
    return Mat4.init(.{
        1, 0,          0,           0,
        0, @cos(rads), -@sin(rads), 0,
        0, @sin(rads), @cos(rads),  0,
        0, 0,          0,           1,
    });
}

pub fn makeRotationY(rads: f64) Mat4 {
    return Mat4.init(.{
        @cos(rads),  0, @sin(rads), 0,
        0,           1, 0,          0,
        -@sin(rads), 0, @cos(rads), 0,
        0,           0, 0,          1,
    });
}

pub fn makeRotationZ(rads: f64) Mat4 {
    return Mat4.init(.{
        @cos(rads), -@sin(rads), 0, 0,
        @sin(rads), @cos(rads),  0, 0,
        0,          0,           1, 0,
        0,          0,           0, 1,
    });
}

pub fn makeShearing(xy: f64, xz: f64, yx: f64, yz: f64, zx: f64, zy: f64) Mat4 {
    return Mat4.init(.{
        1,  xy, xz, 0,
        yx, 1,  yz, 0,
        zx, zy, 1,  0,
        0,  0,  0,  1,
    });
}

const expect = std.testing.expect;

test "Multiplying by a translation matrix" {
    const t = makeTranslation(5, -3, 2);
    const p = Point.init(-3, 4, 5);
    const td = t.mult(p);

    try expect(td.equals(Point.init(2, 1, 7)));
}

test "Multiplying by the inverse of a translation matrix" {
    const t = makeTranslation(5, -3, 2);
    const i = try t.inverted();
    const p = Point.init(-3, 4, 5);
    const td = i.mult(p);

    try expect(td.equals(Point.init(-8, 7, 3)));
}

test "Translation does not affect vectors" {
    const t = makeTranslation(5, -3, 2);
    const v = Vector.init(-3, 4, 5);
    const td = t.mult(v);

    try expect(td.equals(v));
}

test "A scaling matrix applied to a point" {
    const t = makeScaling(2, 3, 4);
    const p = Point.init(-4, 6, 8);
    const td = t.mult(p);

    try expect(td.equals(Point.init(-8, 18, 32)));
}

test "A scaling matrix applied to a vector" {
    const t = makeScaling(2, 3, 4);
    const v = Vector.init(-4, 6, 8);
    const td = t.mult(v);

    try expect(td.equals(Vector.init(-8, 18, 32)));
}

test "Multiplying by the inverse of a scaling matrix" {
    const t = makeScaling(2, 3, 4);
    const i = try t.inverted();
    const v = Vector.init(-4, 6, 8);
    const td = i.mult(v);

    try expect(td.equals(Vector.init(-2, 2, 2)));
}

test "Reflection is scaling by a negative value" {
    const t = makeScaling(-1, 1, 1);
    const v = Point.init(-4, 6, 8);
    const td = t.mult(v);

    try expect(td.equals(Point.init(4, 6, 8)));
}

test "Rotating a point around the x axis" {
    const p = Point.init(0, 1, 0);
    const half_quarter = makeRotationX(pi / 4.0);
    const full_quarter = makeRotationX(pi / 2.0);

    try expect(half_quarter.mult(p).equals(Point.init(0, @sqrt(2.0) / 2.0, @sqrt(2.0) / 2.0)));
    try expect(full_quarter.mult(p).equals(Point.init(0, 0, 1)));
}

test "Rotation a point around the x asix; inverted" {
    const p = Point.init(0, 1, 0);
    const half_quarter = makeRotationX(pi / 4.0);
    const hqi = try half_quarter.inverted();

    const hqid = hqi.mult(p);
    const eq = Point.init(0, @sqrt(2.0) / 2.0, -@sqrt(2.0) / 2.0);

    errdefer std.debug.print("{} != {}\n", .{ hqid.vec, eq.vec });

    try expect(hqid.equals(eq));
}

test "Rotating a point around the y axis" {
    const p = Point.init(0, 0, 1);
    const half_quarter = makeRotationY(pi / 4.0);
    const full_quarter = makeRotationY(pi / 2.0);

    try expect(half_quarter.mult(p).equals(Point.init(@sqrt(2.0) / 2.0, 0, @sqrt(2.0) / 2.0)));
    try expect(full_quarter.mult(p).equals(Point.init(1, 0, 0)));
}

test "Rotation a point around the y asix; inverted" {
    const p = Point.init(0, 0, 1);
    const half_quarter = makeRotationY(pi / 4.0);
    const hqi = try half_quarter.inverted();

    const hqid = hqi.mult(p);
    const eq = Point.init(-@sqrt(2.0) / 2.0, 0, @sqrt(2.0) / 2.0);

    errdefer std.debug.print("{} != {}\n", .{ hqid.vec, eq.vec });

    try expect(hqid.equals(eq));
}

test "Rotating a point around the z axis" {
    const p = Point.init(0, 1, 0);
    const half_quarter = makeRotationZ(pi / 4.0);
    const full_quarter = makeRotationZ(pi / 2.0);

    try expect(half_quarter.mult(p).equals(Point.init(-@sqrt(2.0) / 2.0, @sqrt(2.0) / 2.0, 0)));
    try expect(full_quarter.mult(p).equals(Point.init(-1, 0, 0)));
}

test "Rotation a point around the z asix; inverted" {
    const p = Point.init(0, 1, 0);
    const half_quarter = makeRotationZ(pi / 4.0);
    const hqi = try half_quarter.inverted();

    const hqid = hqi.mult(p);
    const eq = Point.init(@sqrt(2.0) / 2.0, @sqrt(2.0) / 2.0, 0);

    errdefer std.debug.print("{} != {}\n", .{ hqid.vec, eq.vec });

    try expect(hqid.equals(eq));
}

test "A shearing transformation moves x in proportion to y" {
    const p = Point.init(2, 3, 4);
    const t = makeShearing(1, 0, 0, 0, 0, 0);
    const td = t.mult(p);

    try expect(td.equals(Point.init(5, 3, 4)));
}

test "A shearing transformation moves x in proportion to z" {
    const p = Point.init(2, 3, 4);
    const t = makeShearing(0, 1, 0, 0, 0, 0);
    const td = t.mult(p);

    try expect(td.equals(Point.init(6, 3, 4)));
}

test "A shearing transformation moves y in proportion to x" {
    const p = Point.init(2, 3, 4);
    const t = makeShearing(0, 0, 1, 0, 0, 0);
    const td = t.mult(p);

    try expect(td.equals(Point.init(2, 5, 4)));
}

test "A shearing transformation moves y in proportion to z" {
    const p = Point.init(2, 3, 4);
    const t = makeShearing(0, 0, 0, 1, 0, 0);
    const td = t.mult(p);

    try expect(td.equals(Point.init(2, 7, 4)));
}

test "A shearing transformation moves z in proportion to x" {
    const p = Point.init(2, 3, 4);
    const t = makeShearing(0, 0, 0, 0, 1, 0);
    const td = t.mult(p);

    try expect(td.equals(Point.init(2, 3, 6)));
}

test "A shearing transformation moves z in proportion to y" {
    const p = Point.init(2, 3, 4);
    const t = makeShearing(0, 0, 0, 0, 0, 1);
    const td = t.mult(p);

    try expect(td.equals(Point.init(2, 3, 7)));
}

test "Individual transformations are applied in sequence" {
    var p = Point.init(1, 0, 1);
    const rot = makeRotationX(pi / 2.0);
    const scl = makeScaling(5, 5, 5);
    const trn = makeTranslation(10, 5, 7);

    p = rot.mult(p);
    try expect(p.equals(Point.init(1, -1, 0)));

    p = scl.mult(p);
    try expect(p.equals(Point.init(5, -5, 0)));

    p = trn.mult(p);
    try expect(p.equals(Point.init(15, 0, 7)));
}

test "Chained transformations are applied in reverse order" {
    var p = Point.init(1, 0, 1);

    const rot = makeRotationX(pi / 2.0);
    const scl = makeScaling(5, 5, 5);
    const trn = makeTranslation(10, 5, 7);

    const final = trn.mult(scl).mult(rot);

    p = final.mult(p);

    try expect(p.equals(Point.init(15, 0, 7)));
}
