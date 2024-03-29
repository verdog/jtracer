//! materials - phong model

const std = @import("std");

const trans = @import("transform.zig");
const mymath = @import("mymath.zig");

const Tuple = @import("tuple.zig").Tuple;
const Point = @import("tuple.zig").Point;
const Color = @import("color.zig").Color;

pub const TestColor = struct {
    pub fn init() ColorMap {
        return ColorMap{ .@"test" = TestColor{} };
    }

    pub fn at(self: This, point: Tuple) Color {
        _ = self;
        return Color.init(point.x(), point.y(), point.z());
    }

    const This = @This();
};

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

pub const StripeColor = struct {
    a: Color,
    b: Color,

    pub fn init(a: Color, b: Color) ColorMap {
        return ColorMap{ .stripe = StripeColor{
            .a = a,
            .b = b,
        } };
    }

    pub fn initSingle(a: Color) ColorMap {
        return ColorMap{ .stripe = StripeColor{
            .a = a,
            .b = a.scaled(0.5),
        } };
    }

    pub fn at(self: This, point: Tuple) Color {
        return if (std.math.floor(@mod(point.x(), 2)) == 0) self.a else self.b;
    }

    const This = @This();
};

pub const RingColor = struct {
    a: Color,
    b: Color,

    pub fn init(a: Color, b: Color) ColorMap {
        return ColorMap{ .ring = RingColor{
            .a = a,
            .b = b,
        } };
    }

    pub fn initSingle(a: Color) ColorMap {
        return ColorMap{ .ring = RingColor{
            .a = a,
            .b = a.scaled(0.5),
        } };
    }

    pub fn at(self: This, point: Tuple) Color {
        const x = point.x();
        const z = point.z();
        return if (std.math.floor(@mod(@sqrt(x * x + z * z), 2)) == 0) self.a else self.b;
    }

    const This = @This();
};

pub const GradientColor = struct {
    a: Color,
    b: Color,

    pub fn init(a: Color, b: Color) ColorMap {
        return ColorMap{ .gradient = GradientColor{
            .a = a,
            .b = b,
        } };
    }

    pub fn initSingle(a: Color) ColorMap {
        return ColorMap{ .gradient = GradientColor{
            .a = a,
            .b = a.scaled(0.5),
        } };
    }

    pub fn at(self: This, point: Tuple) Color {
        const d = self.b.minus(self.a);
        const fraction = point.x() - std.math.floor(point.x());
        return self.a.plus(d.scaled(fraction));
    }

    const This = @This();
};

pub const TwoDCheckedColor = struct {
    a: Color,
    b: Color,

    pub fn init(a: Color, b: Color) ColorMap {
        return ColorMap{ .twodchecked = TwoDCheckedColor{
            .a = a,
            .b = b,
        } };
    }

    pub fn initSingle(a: Color) ColorMap {
        return ColorMap{ .twodchecked = TwoDCheckedColor{
            .a = a,
            .b = a.scaled(0.5),
        } };
    }

    pub fn at(self: This, point: Tuple) Color {
        const x = std.math.floor(@mod(point.x(), 2));
        const z = std.math.floor(@mod(point.z(), 2));
        return if (x != z) self.a else self.b;
    }

    const This = @This();
};

pub const ThreeDCheckedColor = struct {
    a: Color,
    b: Color,

    pub fn init(a: Color, b: Color) ColorMap {
        return ColorMap{ .threedchecked = ThreeDCheckedColor{
            .a = a,
            .b = b,
        } };
    }

    pub fn initSingle(a: Color) ColorMap {
        return ColorMap{ .threedchecked = ThreeDCheckedColor{
            .a = a,
            .b = a.scaled(0.5),
        } };
    }

    pub fn at(self: This, point: Tuple) Color {
        // we offset the coords a little to fight acne in the common case of a flat object
        // like a cube or plane using this material. in those cases the seams are right on
        // the boundaries.
        const sign = struct {
            fn f(fl: f64) f64 {
                if (std.math.approxEqAbs(f64, fl, 0, mymath.floatTolerance)) return 1;
                return std.math.sign(fl);
            }
        }.f;
        const nudge = 2;
        const x = std.math.floor(point.x() + sign(point.x()) * nudge * mymath.floatTolerance);
        const y = std.math.floor(point.y() + sign(point.y()) * nudge * mymath.floatTolerance);
        const z = std.math.floor(point.z() + sign(point.z()) * nudge * mymath.floatTolerance);
        return if (@mod(x + y + z, 2) == 0) self.a else self.b;
    }

    const This = @This();
};

pub const ColorMap = union(enum) {
    @"test": TestColor,
    flat: FlatColor,
    stripe: StripeColor,
    ring: RingColor,
    gradient: GradientColor,
    twodchecked: TwoDCheckedColor,
    threedchecked: ThreeDCheckedColor,

    pub fn at(self: ColorMap, point: Tuple) Color {
        switch (self) {
            inline else => |inner| return inner.at(point),
        }
    }

    pub fn atT(
        self: ColorMap,
        point: Tuple,
        obj_tfm: trans.Transform,
        mat_tfm: trans.Transform,
    ) Color {
        const obj_space_pt = obj_tfm.inverse.mult(point);
        const mat_space_pt = mat_tfm.inverse.mult(obj_space_pt);
        return self.at(mat_space_pt);
    }

    pub fn atC(self: ColorMap, x: f64, y: f64, z: f64) Color {
        switch (self) {
            inline else => |_| return self.at(Point.init(x, y, z)),
        }
    }
};

pub const Material = struct {
    color_map: ColorMap,
    transform: trans.Transform = trans.Transform{},
    ambient: f64,
    diffuse: f64,
    specular: f64,
    shininess: f64,
    reflective: f64,
    transparency: f64,
    // refractive index examples for various materials:
    // - vacuum: 1
    // - air: 1.00029
    // - water: 1.333
    // - glass: 1.52
    // - diamond: 2.417
    refractive_index: f64,

    pub fn init() This {
        return .{
            .color_map = FlatColor.init(Color.init(1, 1, 1)),
            .ambient = 0.1,
            .diffuse = 0.9,
            .specular = 0.9,
            .shininess = 200.0,
            .reflective = 0.0,
            .transparency = 0.0,
            .refractive_index = 1.0,
        };
    }

    const This = @This();
};

const expect = std.testing.expect;
const print = @import("u.zig").print;

test "The default material" {
    const m = Material.init();

    try expect(m.color_map.atC(0, 0, 0).equals(Color.init(1, 1, 1)));
    try expect(m.ambient == 0.1);
    try expect(m.diffuse == 0.9);
    try expect(m.specular == 0.9);
    try expect(m.shininess == 200.0);
}

test "Test pattern with an object transformation" {
    const c = TestColor.init();
    const obj_tfm = trans.makeScaling(2, 2, 2);
    const mat_tfm = trans.Transform{};

    try expect(c.atT(Point.init(2, 3, 4), obj_tfm, mat_tfm).equals(Color.init(1, 1.5, 2)));
}

test "Test pattern with a pattern transformation" {
    const c = TestColor.init();
    const obj_tfm = trans.Transform{};
    const mat_tfm = trans.makeScaling(2, 2, 2);

    try expect(c.atT(Point.init(2, 3, 4), obj_tfm, mat_tfm).equals(Color.init(1, 1.5, 2)));
}

test "Test pattern with both a pattern and object transformation" {
    const c = TestColor.init();
    const obj_tfm = trans.makeScaling(2, 2, 2);
    const mat_tfm = trans.makeTranslation(0.5, 1, 1.5);

    try expect(c.atT(Point.init(2.5, 3, 3.5), obj_tfm, mat_tfm).equals(Color.init(0.75, 0.5, 0.25)));
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
    const sc = StripeColor.init(
        Color.init(1, 1, 1),
        Color.init(0, 0, 0),
    );

    try expect(sc.stripe.a.equals(Color.init(1, 1, 1)));
    try expect(sc.stripe.b.equals(Color.init(0, 0, 0)));
}

test "Stripe pattern is constant in y" {
    const sc = StripeColor.init(
        Color.init(1, 1, 1),
        Color.init(0, 0, 0),
    );

    try expect(sc.atC(0, 0, 0).equals(Color.init(1, 1, 1)));
    try expect(sc.atC(0, 1, 0).equals(Color.init(1, 1, 1)));
    try expect(sc.atC(0, 2, 0).equals(Color.init(1, 1, 1)));
}

test "Stripe pattern is constant in z" {
    const sc = StripeColor.init(
        Color.init(1, 1, 1),
        Color.init(0, 0, 0),
    );

    try expect(sc.atC(0, 0, 0).equals(Color.init(1, 1, 1)));
    try expect(sc.atC(0, 0, 1).equals(Color.init(1, 1, 1)));
    try expect(sc.atC(0, 0, 2).equals(Color.init(1, 1, 1)));
}

test "Stripe pattern alternates in x" {
    const sc = StripeColor.init(
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

test "Stripes with an object transformation" {
    const obj_tfm = trans.makeScaling(2, 2, 2);
    const mat_tfm = trans.Transform{};
    const sp = StripeColor.init(
        Color.init(1, 1, 1),
        Color.init(0, 0, 0),
    );

    try expect(sp.atT(Point.init(1.5, 0, 0), obj_tfm, mat_tfm).equals(Color.init(1, 1, 1)));
}

test "Stripes with an pattern transformation" {
    const obj_tfm = trans.Transform{};
    const mat_tfm = trans.makeScaling(2, 2, 2);
    const sp = StripeColor.init(
        Color.init(1, 1, 1),
        Color.init(0, 0, 0),
    );

    try expect(sp.atT(Point.init(1.5, 0, 0), obj_tfm, mat_tfm).equals(Color.init(1, 1, 1)));
}

test "Stripes with an object and pattern transformation" {
    const obj_tfm = trans.makeScaling(2, 2, 2);
    const mat_tfm = trans.makeTranslation(0.5, 0, 0);
    const sp = StripeColor.init(
        Color.init(1, 1, 1),
        Color.init(0, 0, 0),
    );

    try expect(sp.atT(Point.init(2.5, 0, 0), obj_tfm, mat_tfm).equals(Color.init(1, 1, 1)));
}

test "Ring pattern should extend in both x and z" {
    const rp = RingColor.init(
        Color.init(1, 1, 1),
        Color.init(0, 0, 0),
    );

    try expect(rp.atC(0, 0, 0).equals(Color.init(1, 1, 1)));
    try expect(rp.atC(1, 0, 0).equals(Color.init(0, 0, 0)));
    try expect(rp.atC(0, 0, 1).equals(Color.init(0, 0, 0)));
    // slightly more than @sqrt(2)
    try expect(rp.atC(0.708, 0, 0.708).equals(Color.init(0, 0, 0)));
}

test "A gradient linearly interpolates between colors" {
    const gp = GradientColor.init(
        Color.init(1, 1, 1),
        Color.init(0, 0, 0),
    );

    {
        const c = gp.atC(0, 0, 0);
        errdefer print(c);
        try expect(c.equals(Color.init(1, 1, 1)));
    }
    {
        const c = gp.atC(0.25, 0, 0);
        errdefer print(c);
        try expect(c.equals(Color.init(0.75, 0.75, 0.75)));
    }
    {
        const c = gp.atC(0.5, 0, 0);
        errdefer print(c);
        try expect(c.equals(Color.init(0.5, 0.5, 0.5)));
    }
    {
        const c = gp.atC(0.75, 0, 0);
        errdefer print(c);
        try expect(c.equals(Color.init(0.25, 0.25, 0.25)));
    }
    {
        const c = gp.atC(1, 0, 0);
        errdefer print(c);
        try expect(c.equals(Color.init(1, 1, 1)));
    }
}

test "A TwoDCheckedColor is checkered in xz but constant in y" {
    const cp = TwoDCheckedColor.init(
        Color.init(1, 1, 1),
        Color.init(0, 0, 0),
    );

    {
        const c = cp.atC(0.25, 0, 0.25);
        errdefer print(c);
        try expect(c.equals(Color.init(0, 0, 0)));
    }
    {
        const c = cp.atC(-0.25, 0, 0.25);
        errdefer print(c);
        try expect(c.equals(Color.init(1, 1, 1)));
    }
    {
        const c = cp.atC(-0.25, 0, -0.25);
        errdefer print(c);
        try expect(c.equals(Color.init(0, 0, 0)));
    }
    {
        const c = cp.atC(0.25, 0, -0.25);
        errdefer print(c);
        try expect(c.equals(Color.init(1, 1, 1)));
    }
    {
        const c = cp.atC(0.25, 1.25, 0.25);
        errdefer print(c);
        try expect(c.equals(Color.init(0, 0, 0)));
    }
}

test "A ThreeDCheckedColor is checkered in all dimensions" {
    const cp = ThreeDCheckedColor.init(
        Color.init(1, 1, 1),
        Color.init(0, 0, 0),
    );

    {
        const c = cp.atC(0.25, 0, 0.25);
        errdefer print(c);
        try expect(c.equals(Color.init(1, 1, 1)));
    }
    {
        const c = cp.atC(-0.25, 0, 0.25);
        errdefer print(c);
        try expect(c.equals(Color.init(0, 0, 0)));
    }
    {
        const c = cp.atC(-0.25, 0, -0.25);
        errdefer print(c);
        try expect(c.equals(Color.init(1, 1, 1)));
    }
    {
        const c = cp.atC(0.25, 0, -0.25);
        errdefer print(c);
        try expect(c.equals(Color.init(0, 0, 0)));
    }
    {
        const c = cp.atC(0.25, 1.25, 0.25);
        errdefer print(c);
        try expect(c.equals(Color.init(0, 0, 0)));
    }
}

test "A ThreeDCheckedColor should repeat in x" {
    const cp = ThreeDCheckedColor.init(
        Color.init(1, 1, 1),
        Color.init(0, 0, 0),
    );

    {
        const c = cp.atC(0, 0, 0);
        errdefer print(c);
        try expect(c.equals(Color.init(1, 1, 1)));
    }
    {
        const c = cp.atC(0.99, 0, 0);
        errdefer print(c);
        try expect(c.equals(Color.init(1, 1, 1)));
    }
    {
        const c = cp.atC(1.01, 0, 0);
        errdefer print(c);
        try expect(c.equals(Color.init(0, 0, 0)));
    }
}

test "A ThreeDCheckedColor should repeat in y" {
    const cp = ThreeDCheckedColor.init(
        Color.init(1, 1, 1),
        Color.init(0, 0, 0),
    );

    {
        const c = cp.atC(0, 0, 0);
        errdefer print(c);
        try expect(c.equals(Color.init(1, 1, 1)));
    }
    {
        const c = cp.atC(0, 0.99, 0);
        errdefer print(c);
        try expect(c.equals(Color.init(1, 1, 1)));
    }
    {
        const c = cp.atC(0, 1.01, 0);
        errdefer print(c);
        try expect(c.equals(Color.init(0, 0, 0)));
    }
}

test "A ThreeDCheckedColor should repeat in z" {
    const cp = ThreeDCheckedColor.init(
        Color.init(1, 1, 1),
        Color.init(0, 0, 0),
    );

    {
        const c = cp.atC(0, 0, 0);
        errdefer print(c);
        try expect(c.equals(Color.init(1, 1, 1)));
    }
    {
        const c = cp.atC(0, 0, 0.99);
        errdefer print(c);
        try expect(c.equals(Color.init(1, 1, 1)));
    }
    {
        const c = cp.atC(0, 0, 1.01);
        errdefer print(c);
        try expect(c.equals(Color.init(0, 0, 0)));
    }
}

test "Reflectivity for the default material" {
    const m = Material.init();
    try expect(m.reflective == 0.0);
}

test "Transparency and refractive index for the default material" {
    const m = Material.init();
    try expect(m.transparency == 0.0);
    try expect(m.refractive_index == 1.0);
}
