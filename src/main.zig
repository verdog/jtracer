const std = @import("std");
const sdl2 = @import("sdl2");

const imgio = @import("imgio.zig");
const trans = @import("transform.zig");
const rndr = @import("render.zig");

const World = @import("world.zig").World;
const Ray = @import("ray.zig").Ray;
const Sphere = @import("sphere.zig").Sphere;
const PointLight = @import("light.zig").PointLight;

const Point = @import("tuple.zig").Point;
const Vector = @import("tuple.zig").Vector;
const Qanvas = @import("qanvas.zig").Qanvas;
const Color = @import("color.zig").Color;

var gpa_impl = std.heap.GeneralPurposeAllocator(.{}){};
const gpa = gpa_impl.allocator();

pub fn main() !void {
    defer if (!gpa_impl.detectLeaks()) std.debug.print("(No leaks)\n", .{});

    var qan = try Qanvas.init(gpa, 1024, 512);
    defer qan.deinit();

    try sdl2.init(.{
        .video = true,
        .events = true,
    });
    defer sdl2.quit();
    try sdl2.image.init(.{ .png = true });
    defer sdl2.image.quit();

    // create window
    var sdl_window = try sdl2.createWindow(
        "jtracer",
        .{ .centered = {} },
        .{ .centered = {} },
        qan.width,
        qan.height,
        .{ .vis = .shown },
    );
    defer sdl_window.destroy();

    var sdl_renderer = try sdl2.createRenderer(sdl_window, null, .{ .accelerated = true });
    defer sdl_renderer.destroy();

    // pipeline is
    // qanvas -> qixels -> sdl texture

    var qix_tex = try imgio.encodeQanvasAsQixel(qan, gpa);
    defer gpa.free(qix_tex);

    var sdl_tex = try imgio.encodeQixelsAsSdl(qix_tex, sdl_renderer, @intCast(u32, qan.width), @intCast(u32, qan.height));
    defer sdl_tex.destroy();

    var frame: f64 = 0;

    {
        var world = World.init(gpa);
        defer world.deinit();

        var lgt = world.add(PointLight);
        var lgtptr = world.get(PointLight, lgt);
        lgtptr.position = Point.init(-0.5, 1.2, -1.1);
        lgtptr.intensity = Color.init(1, 1, 1);

        var sph = world.add(Sphere);
        var sphptr = world.get(Sphere, sph);
        sphptr.material.color = Color.init(1, 0, 0);
        sphptr.transform = sphptr.transform.mult(trans.makeTranslation(-3, 0, 1));
        sphptr.transform = sphptr.transform.mult(trans.makeScaling(1, 2, 1));

        sph = world.add(Sphere);
        sphptr = world.get(Sphere, sph);
        sphptr.material.color = Color.init(0, 1, 0);
        sphptr.transform = sphptr.transform.mult(trans.makeTranslation(0, 0, 1));
        sphptr.transform = sphptr.transform.mult(trans.makeScaling(1, 1, 1));

        sph = world.add(Sphere);
        sphptr = world.get(Sphere, sph);
        sphptr.material.color = Color.init(0, 0, 1);
        sphptr.transform = sphptr.transform.mult(trans.makeTranslation(3, 0, 1));
        sphptr.transform = sphptr.transform.mult(trans.makeRotationZ(std.math.pi / 4.0));
        sphptr.transform = sphptr.transform.mult(trans.makeScaling(1, 0.4, 1));

        const window_center = Point.init(0, 0, -3);
        const heading = Vector.init(0, 0, 1);

        // zig fmt: off
        const render_thread = std.Thread.spawn(.{}, rndr.startRenderEngine,
            .{world, window_center, heading, &qan, gpa}
        ) catch unreachable;
        render_thread.detach();
        // zig fmt: on

        main_loop: while (true) {
            defer {
                frame += 1;
            }

            while (sdl2.pollEvent()) |ev| {
                switch (ev) {
                    .quit => break :main_loop,
                    .key_up => |key_up| {
                        if (key_up.keycode == .escape) break :main_loop;
                    },
                    else => {
                        // std.debug.print("Unhandled event: {}\n\n", .{ev});
                    },
                }
            }

            try imgio.encodeQanvasAsQixelUpdate(qan, qix_tex);
            try imgio.encodeQixelsAsSdlUpdate(qix_tex, &sdl_tex, @intCast(u32, qan.width));

            try sdl_renderer.copy(sdl_tex, null, null);
            sdl_renderer.present();

            sdl2.delay(17);
        }
    }
}

test {
    _ = @import("tuple.zig");
    _ = @import("color.zig");
    _ = @import("qanvas.zig");
    _ = @import("qoi.zig");
    _ = @import("imgio.zig");
    _ = @import("matrix.zig");
    _ = @import("transform.zig");
    _ = @import("ray.zig");
    _ = @import("sphere.zig");
    _ = @import("world.zig");
    _ = @import("intersect.zig");
    _ = @import("light.zig");
    _ = @import("material.zig");
    _ = @import("end2end.zig");
}
