//! holds all volumes

const std = @import("std");

pub const World = struct {
    pub fn init(alctr: std.mem.Allocator) This {
        return .{
            .spheres_buf = std.ArrayList(Sphere).init(alctr),
            .lights_buf = std.ArrayList(PointLight).init(alctr),
        };
    }

    pub fn deinit(self: This) void {
        self.spheres_buf.deinit();
        self.lights_buf.deinit();
    }

    pub fn add(self: *This, comptime T: type) VolumePtr {
        switch (T) {
            Sphere => {
                self.spheres_buf.append(Sphere.init()) catch unreachable;
                return VolumePtr{ .sphere_idx = self.spheres_buf.items.len - 1 };
            },
            PointLight => {
                const point = @import("tuple.zig").Point.init(0, 0, 0);
                const color = @import("color.zig").Color.init(1, 1, 1);
                self.lights_buf.append(PointLight.init(point, color)) catch unreachable;
                return VolumePtr{ .light_idx = self.lights_buf.items.len - 1 };
            },
            else => unreachable,
        }
    }

    pub fn get(self: This, comptime T: type, vptr: VolumePtr) *T {
        return switch (T) {
            Sphere => return &self.spheres_buf.items[vptr.sphere_idx],
            PointLight => return &self.lights_buf.items[vptr.light_idx],
            else => unreachable,
        };
    }

    // TODO change these to a list that won't move data
    // around upon extension (e.g. a chunky linked list)?
    spheres_buf: std.ArrayList(Sphere),
    lights_buf: std.ArrayList(PointLight),

    const This = @This();
};

const Sphere = @import("sphere.zig").Sphere;

pub const VolumePtr = union(enum) {
    none: void,
    sphere_idx: usize,
    light_idx: usize,
};

const PointLight = @import("light.zig").PointLight;

const expect = std.testing.expect;

test "World creation/destruction" {
    const alctr = std.testing.allocator;
    var world = World.init(alctr);
    defer world.deinit();
}

test "World: add/get sphere" {
    const alctr = std.testing.allocator;
    var world = World.init(alctr);
    defer world.deinit();

    const sphereidx = world.add(Sphere);
    var sphereptr = world.get(Sphere, sphereidx);

    try expect(sphereptr.id != 0);
}
