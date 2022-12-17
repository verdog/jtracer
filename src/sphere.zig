const std = @import("std");

const mtx = @import("matrix.zig");
const trans = @import("transform.zig");

pub const Sphere = struct {
    var next_id: usize = 1;

    id: usize,
    transform: mtx.Matrix(4, 4),

    pub fn init() This {
        const id = This.next_id;
        This.next_id += 1;
        return .{
            .id = id,
            .transform = mtx.Matrix(4, 4).identity(),
        };
    }

    const This = @This();
};

const expect = std.testing.expect;

test "Spheres have unique ids" {
    const s1 = Sphere.init();
    const s2 = Sphere.init();

    try expect(s1.id != s2.id);
}

test "Sphere's have a default transform" {
    const s = Sphere.init();

    try expect(s.transform.equals(mtx.Matrix(4, 4).identity()));
}

test "Changing a sphere's transformation" {
    var s = Sphere.init();
    var t = trans.makeTranslation(2, 3, 4);
    s.transform = t;

    try expect(s.transform.equals(t));
}
