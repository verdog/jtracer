//! materials - phong model

const std = @import("std");

const Tuple = @import("tuple.zig").Tuple;
const Color = @import("color.zig").Color;

pub const Material = struct {
    color: Color,
    ambient: f64,
    diffuse: f64,
    specular: f64,
    shininess: f64,

    pub fn init() This {
        return .{
            .color = Color.init(1, 1, 1),
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

    try expect(m.color.equals(Color.init(1, 1, 1)));
    try expect(m.ambient == 0.1);
    try expect(m.diffuse == 0.9);
    try expect(m.specular == 0.9);
    try expect(m.shininess == 200.0);
}
