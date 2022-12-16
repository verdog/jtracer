const std = @import("std");
const sdl2 = @import("sdl2");

const imgio = @import("imgio.zig");
const tlate = @import("transform.zig");

const Point = @import("tuple.zig").Point;
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

        {
            qan.clear();

            const tick = 256;
            const rot = tlate.makeRotationZ(2 * std.math.pi / @intToFloat(f64, tick));
            const shrF = @sin(frame / 20.0);
            const shr = tlate.makeShearing(0, 0, shrF, 0, 0, 0);
            const anchor = tlate.makeTranslation(512 + 128 * @cos(frame / 4), 256 + 128 * @sin(frame / 4), 0);
            var clock_hand = Point.init(0, -128, 0);

            var i: usize = 0;
            while (i < tick) : (i += 1) {
                const p = anchor.mult(shr).mult(clock_hand);
                qan.write(Color.init(0, 1, 1), @floatToInt(i64, p.x()), @floatToInt(i64, p.y()));
                clock_hand = rot.mult(clock_hand);
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
    _ = @import("intersect.zig");
}
