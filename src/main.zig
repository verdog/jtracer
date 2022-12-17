const std = @import("std");
const sdl2 = @import("sdl2");

const imgio = @import("imgio.zig");
const trans = @import("transform.zig");
const intx = @import("intersect.zig");

const World = @import("world.zig").World;
const Ray = @import("ray.zig").Ray;
const Sphere = @import("sphere.zig").Sphere;

const Point = @import("tuple.zig").Point;
const Vector = @import("tuple.zig").Vector;
const Qanvas = @import("qanvas.zig").Qanvas;
const Color = @import("color.zig").Color;

var gpa_impl = std.heap.GeneralPurposeAllocator(.{}){};
const gpa = gpa_impl.allocator();

pub fn main() !void {
    defer if (!gpa_impl.detectLeaks()) std.debug.print("(No leaks)\n", .{});

    var qan = try Qanvas.init(gpa, 600, 400);
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

    {
        std.debug.print("rendering...\n", .{});

        var prog_ctx = std.Progress{};
        var prog = prog_ctx.start("Pixels", qan.width * qan.height);

        var world = World.init(gpa);
        defer world.deinit();

        var sph = world.add(Sphere);
        var sphptr = world.get(Sphere, sph);

        const camera = Point.init(0, 0, -3);
        const heading = Vector.init(0, 0, 1).normalized();
        const window_center = camera.plus(heading);

        var pix_y: i64 = 0;
        while (pix_y < qan.height) : (pix_y += 1) {
            var pix_x: i64 = 0;
            while (pix_x < qan.width) : (pix_x += 1) {
                const pixels_per_unit = 512;
                // get point on window
                const pix_x_w = @intToFloat(f64, pix_x - @divTrunc(@intCast(i64, qan.width), 2)) / pixels_per_unit;
                const pix_y_w = @intToFloat(f64, pix_y - @divTrunc(@intCast(i64, qan.height), 2)) / pixels_per_unit;
                // note the -pix_y_w. pixel space y increases downwards
                // but object space y increases upwards
                const aim_point = window_center.plus(Vector.init(pix_x_w, -pix_y_w, 0));
                const aim_ray = Ray.init(camera, aim_point.minus(camera));

                const xs = intx.intersect(sph, sphptr.*, aim_ray, gpa);
                defer xs.deinit();

                // std.debug.print("{d: <4.3} {d: <4.3}", .{ aim_x, aim_y });
                if (xs.hit()) |_| {
                    qan.write(Color.init(1, 0, 0), pix_x, pix_y);
                }

                prog.completeOne();
            }
        }

        prog.end();

        std.debug.print("done.\n", .{});
    }

    // pipeline is
    // qanvas -> qixels -> sdl texture

    var qix_tex = try imgio.encodeQanvasAsQixel(qan, gpa);
    defer gpa.free(qix_tex);

    var sdl_tex = try imgio.encodeQixelsAsSdl(qix_tex, sdl_renderer, @intCast(u32, qan.width), @intCast(u32, qan.height));
    defer sdl_tex.destroy();

    var frame: f64 = 0;

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
    _ = @import("end2end.zig");
}
