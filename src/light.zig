//! light sources

const std = @import("std");

const Color = @import("color.zig").Color;
const Point = @import("tuple.zig").Point;
const Vector = @import("tuple.zig").Vector;
const Tuple = @import("tuple.zig").Tuple;

const mat = @import("material.zig");
const Material = @import("material.zig").Material;
const trans = @import("transform.zig");

pub const PointLight = struct {
    position: Tuple,
    intensity: Color,

    pub fn init(pos: Tuple, int: Color) This {
        return .{ .position = pos, .intensity = int };
    }

    const This = @This();
};

/// input vectors are:
/// - p : The point on the volume
/// - e : The eye vector (point -> eye)
/// - n : The normal vector at p
pub fn lighting(
    m: Material,
    obj_tfm: trans.Transform,
    l: PointLight,
    p: Tuple,
    e: Tuple,
    n: Tuple,
    in_shadow: bool,
) Color {
    const effective_color = m.color_map.atT(p, obj_tfm, m.transform).multiplied(l.intensity);

    var lightv = l.position.minus(p);
    if (lightv.magnitude() == 0) {
        // light is sitting on the surface. pretend the light direction
        // is the same as the normal vector.
        lightv = n;
    }
    lightv = lightv.normalized();

    const ambient = effective_color.scaled(m.ambient);

    var diffuse = Color.init(0, 0, 0);
    var specular = Color.init(0, 0, 0);

    const ldn = lightv.dot(n);
    if (ldn < 0 or in_shadow) {
        // light is behind the surface
        diffuse = diffuse.scaled(0);
        specular = specular.scaled(0);
    } else {
        // light is in front of surface
        diffuse = effective_color.scaled(m.diffuse).scaled(ldn);

        const reflectv = lightv.scaled(-1).reflected(n);
        const rde = reflectv.dot(e);
        if (rde <= 0) {
            // light is reflected away from eye
            specular = specular.scaled(0);
        } else {
            const factor = std.math.pow(f64, rde, m.shininess);
            specular = l.intensity.scaled(m.specular).scaled(factor);
        }
    }

    return ambient.plus(diffuse).plus(specular);
}

const expect = std.testing.expect;

test "A PointLight has position and intensity" {
    const int = Color.init(1, 1, 1);
    const pos = Point.init(0, 0, 0);
    const light = PointLight.init(pos, int);

    try expect(light.position.equals(pos));
    try expect(light.intensity.equals(int));
}

test "Lighting with the eye between the light and the surface" {
    const m = Material.init();
    const p = Point.init(0, 0, 0);

    const eye = Vector.init(0, 0, -1);
    const n = Vector.init(0, 0, -1);
    const l = PointLight.init(Point.init(0, 0, -10), Color.init(1, 1, 1));

    try expect(lighting(m, trans.Transform{}, l, p, eye, n, false).equals(Color.init(1.9, 1.9, 1.9)));
}

test "Lighting with the surface in the shadow" {
    const m = Material.init();
    const p = Point.init(0, 0, 0);

    const eye = Vector.init(0, 0, -1);
    const n = Vector.init(0, 0, -1);
    const l = PointLight.init(Point.init(0, 0, -10), Color.init(1, 1, 1));
    const in_shadow = true;

    try expect(lighting(m, trans.Transform{}, l, p, eye, n, in_shadow).equals(Color.init(0.1, 0.1, 0.1)));
}

test "Lighting with the eye between light and surface, eye offset 45d" {
    const m = Material.init();
    const p = Point.init(0, 0, 0);

    const eye = Vector.init(0, @sqrt(2.0) / 2.0, -@sqrt(2.0) / 2.0);
    const n = Vector.init(0, 0, -1);
    const l = PointLight.init(Point.init(0, 0, -10), Color.init(1, 1, 1));

    try expect(lighting(m, trans.Transform{}, l, p, eye, n, false).equals(Color.init(1.0, 1.0, 1.0)));
}

test "Lighting with eye opposite surface, light offset 45d" {
    const m = Material.init();
    const p = Point.init(0, 0, 0);

    const eye = Vector.init(0, 0, -1);
    const n = Vector.init(0, 0, -1);
    const l = PointLight.init(Point.init(0, 10, -10), Color.init(1, 1, 1));

    // book examples are much less precise than f64
    try expect(lighting(m, trans.Transform{}, l, p, eye, n, false).equalsTolerance(Color.init(0.7364, 0.7364, 0.7364), 100_000_000_000));
}

test "Lighting with eye in the path of the reflection vector" {
    const m = Material.init();
    const p = Point.init(0, 0, 0);

    const eye = Vector.init(0, -@sqrt(2.0) / 2.0, -@sqrt(2.0) / 2.0);
    const n = Vector.init(0, 0, -1);
    const l = PointLight.init(Point.init(0, 10, -10), Color.init(1, 1, 1));

    // book examples are much less precise than f64
    try expect(lighting(m, trans.Transform{}, l, p, eye, n, false).equalsTolerance(Color.init(1.6364, 1.6364, 1.6364), 100_000_000_000));
}

test "Lighting with the light behind the surface" {
    const m = Material.init();
    const p = Point.init(0, 0, 0);

    const eye = Vector.init(0, 0, -1);
    const n = Vector.init(0, 0, -1);
    const l = PointLight.init(Point.init(0, 0, 10), Color.init(1, 1, 1));

    try expect(lighting(m, trans.Transform{}, l, p, eye, n, false).equals(Color.init(0.1, 0.1, 0.1)));
}

test "Lighting with a pattern applied" {
    var m = Material.init();
    m.color_map = mat.StripeColor.init(Color.init(1, 1, 1), Color.init(0, 0, 0));
    m.ambient = 1;
    m.diffuse = 0;
    m.specular = 0;

    const eye = Vector.init(0, 0, -1);
    const n = Vector.init(0, 0, -1);
    const l = PointLight.init(Point.init(0, 0, -10), Color.init(1, 1, 1));

    {
        const c = lighting(m, trans.Transform{}, l, Point.init(0.9, 0, 0), eye, n, false);
        try expect(c.equals(Color.init(1, 1, 1)));
    }
    {
        const c = lighting(m, trans.Transform{}, l, Point.init(1.1, 0, 0), eye, n, false);
        try expect(c.equals(Color.init(0, 0, 0)));
    }
}
