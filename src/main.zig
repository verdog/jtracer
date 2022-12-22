const std = @import("std");
const sdl2 = @import("sdl2");

const imgio = @import("imgio.zig");
const trans = @import("transform.zig");
const rndr = @import("render.zig");

const World = @import("world.zig").World;
const Ray = @import("ray.zig").Ray;
const Sphere = @import("sphere.zig").Sphere;
const PointLight = @import("light.zig").PointLight;
const Camera = @import("world.zig").Camera;

const Point = @import("tuple.zig").Point;
const Vector = @import("tuple.zig").Vector;
const Qanvas = @import("qanvas.zig").Qanvas;
const Color = @import("color.zig").Color;

var gpa_impl = std.heap.GeneralPurposeAllocator(.{}){};
const gpa = gpa_impl.allocator();

pub fn main() !void {
    defer if (!gpa_impl.detectLeaks()) std.debug.print("(No leaks)\n", .{});

    var qan = try Qanvas.init(gpa, 1200, 800);
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

        var lgt = world.addLight(PointLight);
        lgt.ptr.position = Point.init(-3, 1, 1);
        lgt.ptr.intensity = Color.init(0.2, 0.2, 0.4);

        lgt = world.addLight(PointLight);
        lgt.ptr.position = Point.init(3, 1, 3);
        lgt.ptr.intensity = Color.init(0.2, 0.4, 0.2);

        lgt = world.addLight(PointLight);
        lgt.ptr.position = Point.init(0, 3, 2);
        lgt.ptr.intensity = Color.init(0.6, 0.2, 0.2);

        {
            var x: f64 = -2.5;
            while (x <= 2.5) : (x += 0.4) {
                var z: f64 = -0.5;
                while (z <= 4.5) : (z += 0.4) {
                    var sph = world.addVolume(Sphere);
                    sph.ptr.material.color = Color.init(0.75, 0.75, 1);
                    sph.ptr.transform = sph.ptr.transform.chain(.{
                        trans.makeTranslation(x, @sin(x + z) / 4 + z / 2, z),
                        trans.makeRotationX(-std.math.pi / 3.0),
                        trans.makeScaling(0.1, 0.05, 0.1),
                    });
                }
            }
        }

        var cam = Camera.init(@intCast(i64, qan.width), @intCast(i64, qan.height), std.math.pi / 2.0);
        const from = Point.init(0, 1.5, -1);
        const to = Point.init(0, 0.5, 1);
        const up = Vector.init(0, 1, 0);
        cam.transform = trans.makeView(from, to, up);

        var thr = std.Thread.spawn(.{}, World.render, .{ world, cam, &qan, gpa }) catch unreachable;
        thr.detach();

        main_loop: while (true) {
            defer {
                frame += 1;
            }

            while (sdl2.pollEvent()) |ev| {
                switch (ev) {
                    .quit => break :main_loop,
                    .key_up => |key_up| {
                        if (key_up.keycode == .escape) break :main_loop;
                        if (key_up.keycode == .s) {
                            std.debug.print("Saving... ", .{});
                            const buf = try imgio.encodeQanvasAsQoi(qan, gpa);
                            defer gpa.free(buf);
                            try imgio.saveBufAsFile(buf, "out.qoi");
                            std.debug.print(" {s}\n", .{"out.qoi"});
                        }
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
