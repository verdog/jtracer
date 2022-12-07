const std = @import("std");

pub const Sphere = struct {
    var next_id: usize = 1;

    id: usize,

    pub fn init() This {
        const id = This.next_id;
        This.next_id += 1;
        return .{ .id = id };
    }

    const This = @This();
};

const expect = std.testing.expect;

test "Spheres have unique ids" {
    const s1 = Sphere.init();
    const s2 = Sphere.init();

    try expect(s1.id != s2.id);
}
