//! materials - phong model

const std = @import("std");

const Tuple = @import("tuple.zig").Tuple;
const Point = @import("tuple.zig").Point;
const Color = @import("color.zig").Color;

pub const FlatColor = struct {
    color: Color,

    pub fn init(color: Color) ColorMap {
        return ColorMap{ .flat = FlatColor{ .color = color } };
    }

    pub fn at(self: This, point: Tuple) Color {
        _ = point; // flat color is always the same
        return self.color;
    }

    const This = @This();
};

pub const StripedColor = struct {
    a: Color,
    b: Color,

    pub fn init(a: Color, b: Color) ColorMap {
        return ColorMap{ .striped = StripedColor{
            .a = a,
            .b = b,
        } };
    }

    pub fn at(self: This, point: Tuple) Color {
        return if (std.math.floor(@mod(point.x(), 2)) == 0) self.a else self.b;
    }

    const This = @This();
};

pub const ColorMap = union(enum) {
    flat: FlatColor,
    striped: StripedColor,

    pub fn at(self: ColorMap, point: Tuple) Color {
        switch (self) {
            inline else => |inner| return inner.at(point),
        }
    }

    pub fn atC(self: ColorMap, x: f64, y: f64, z: f64) Color {
        switch (self) {
            inline else => |_| return self.at(Point.init(x, y, z)),
        }
    }
};

pub const Material = struct {
    color_map: ColorMap,
    ambient: f64,
    diffuse: f64,
    specular: f64,
    shininess: f64,

    pub fn init() This {
        return .{
            .color_map = FlatColor.init(Color.init(1, 1, 1)),
            .ambient = 0.1,
            .diffuse = 0.9,
            .specular = 0.9,
            .shininess = 200.0,
        };
    }

    const This = @This();
};

const expect = std.testing.expect;

test "The default material" {
    const m = Material.init();

    try expect(m.color_map.atC(0, 0, 0).equals(Color.init(1, 1, 1)));
    try expect(m.ambient == 0.1);
    try expect(m.diffuse == 0.9);
    try expect(m.specular == 0.9);
    try expect(m.shininess == 200.0);
}

test "Flat color is the same everywhere" {
    const c = FlatColor.init(Color.init(1, 1, 1));

    try expect(c.at(Point.init(0, 0, 0)).equals(Color.init(1, 1, 1)));
    try expect(c.at(Point.init(1, 0, 0)).equals(Color.init(1, 1, 1)));
    try expect(c.at(Point.init(1, 1, 0)).equals(Color.init(1, 1, 1)));
    try expect(c.at(Point.init(1, 1, 1)).equals(Color.init(1, 1, 1)));
    try expect(c.at(Point.init(0.9, 20, 123)).equals(Color.init(1, 1, 1)));
}

test "Stripe pattern" {
    const sc = StripedColor.init(
        Color.init(1, 1, 1),
        Color.init(0, 0, 0),
    );

    try expect(sc.striped.a.equals(Color.init(1, 1, 1)));
    try expect(sc.striped.b.equals(Color.init(0, 0, 0)));
}

test "Stripe pattern is constant in y" {
    const sc = StripedColor.init(
        Color.init(1, 1, 1),
        Color.init(0, 0, 0),
    );

    try expect(sc.atC(0, 0, 0).equals(Color.init(1, 1, 1)));
    try expect(sc.atC(0, 1, 0).equals(Color.init(1, 1, 1)));
    try expect(sc.atC(0, 2, 0).equals(Color.init(1, 1, 1)));
}

test "Stripe pattern is constant in z" {
    const sc = StripedColor.init(
        Color.init(1, 1, 1),
        Color.init(0, 0, 0),
    );

    try expect(sc.atC(0, 0, 0).equals(Color.init(1, 1, 1)));
    try expect(sc.atC(0, 0, 1).equals(Color.init(1, 1, 1)));
    try expect(sc.atC(0, 0, 2).equals(Color.init(1, 1, 1)));
}

test "Stripe pattern alternates in x" {
    const sc = StripedColor.init(
        Color.init(1, 1, 1),
        Color.init(0, 0, 0),
    );

    try expect(sc.atC(0, 0, 0).equals(Color.init(1, 1, 1)));
    try expect(sc.atC(0.9, 0, 0).equals(Color.init(1, 1, 1)));
    try expect(sc.atC(1, 0, 0).equals(Color.init(0, 0, 0)));
    try expect(sc.atC(-0.1, 0, 0).equals(Color.init(0, 0, 0)));
    try expect(sc.atC(-1, 0, 0).equals(Color.init(0, 0, 0)));
    try expect(sc.atC(-1.1, 0, 0).equals(Color.init(1, 1, 1)));
}
