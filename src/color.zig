//! color type. a color is a 3-tuple of floats in the range [0, 1]. the color values can
//! leave the [0, 1] range during math operations, but when displayed, will be clamped back
//! down. the order of the tuple items is red, green, blue.

const std = @import("std");

const mymath = @import("mymath.zig");

pub const Color = struct {
    pub fn init(r: f64, g: f64, b: f64) This {
        return .{ .vec = .{ r, g, b } };
    }

    pub fn equals(self: This, other: This) bool {
        const diff = @fabs(self.vec - other.vec);
        const max_diff = @reduce(.Max, diff);
        return max_diff <= mymath.floatTolerance;
    }

    pub fn red(self: This) f64 {
        return self.vec[0];
    }

    pub fn green(self: This) f64 {
        return self.vec[1];
    }

    pub fn blue(self: This) f64 {
        return self.vec[2];
    }

    pub fn plus(self: This, other: This) This {
        return This{ .vec = self.vec + other.vec };
    }

    pub fn minus(self: This, other: This) This {
        return This{ .vec = self.vec - other.vec };
    }

    pub fn scaled(self: This, f: f64) This {
        return This{ .vec = self.vec * @splat(3, f) };
    }

    pub fn multiplied(self: This, other: This) This {
        return This{ .vec = self.vec * other.vec };
    }

    vec: @Vector(3, f64),

    const This = @This();
};

const expect = std.testing.expect;

test "A color is a r, g, b tuple" {
    const c = Color.init(-0.5, 0.4, 1.7);

    try expect(c.red() == -0.5);
    try expect(c.green() == 0.4);
    try expect(c.blue() == 1.7);
}

test "Testing Equality" {
    try expect(Color.init(1.1, 2.2, 3.3).equals(Color.init(1.1, 2.2, 3.3)));
    try expect(!Color.init(1.0, 2.2, 3.3).equals(Color.init(1.1, 2.2, 3.3)));
    try expect(!Color.init(1.00000000000001, 2.2000000000002, 3.30000000000000003).equals(Color.init(1.1, 2.2, 3.3)));

    const c = Color.init(1.1, 4.4, 6.0);
    const e = std.math.floatEps(f64);
    var c2 = c;

    c2 = c2.plus(Color.init(e, 0, 0));
    try expect(c.equals(c2));
    try expect(c2.equals(c));

    c2 = c2.plus(Color.init(0, e, 0));
    try expect(c.equals(c2));
    try expect(c2.equals(c));

    c2 = c2.plus(Color.init(0, 0, e));
    try expect(c.equals(c2));
    try expect(c2.equals(c));
}

test "Adding colors" {
    const c = Color.init(-0.5, 0.4, 1.7);
    const c2 = Color.init(-0.7, 0.2, 0.3);

    const sum = c.plus(c2);
    errdefer std.debug.print("sum: {}\n", .{sum});

    try expect(sum.equals(Color.init(-1.2, 0.6, 2.0)));
}

test "Subtracting colors" {
    const c = Color.init(-0.5, 0.4, 1.7);
    const c2 = Color.init(-0.7, 0.2, 0.3);

    try expect(c.minus(c2).equals(Color.init(0.2, 0.2, 1.4)));
}

test "Scalar multiplication of colors" {
    const c = Color.init(-0.5, 0.4, 1.7);
    const c2 = Color.init(-1.0, 0.8, 3.4);

    try expect(c.scaled(2).equals(c2));
}

test "hadamard product (color multiplication)" {
    const c = Color.init(1, 0.2, 0.4);
    const c2 = Color.init(-0.9, 1.0, 0.1);

    try expect(c.multiplied(c2).equals(Color.init(-0.9, 0.2, 0.04)));
}
