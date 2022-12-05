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

    {
        var i: usize = 0;
        while (i < qan.width and i < qan.height) : (i += 1) {
            qan.write(Color.init(0, 1, 0), i, i);
            qan.write(Color.init(1, 1, 0), i, qan.height - i - 1);
        }
    }

    {
        const tick = 512;
        const rot = tlate.makeRotationZ(2 * std.math.pi / @intToFloat(f64, tick));
        const shr = tlate.makeShearing(0, 0, 1, 0, 0, 0);
        const anchor = tlate.makeTranslation(512, 256, 0);
        var clock_hand = Point.init(0, -128, 0);

        var i: usize = 0;
        while (i < tick) : (i += 1) {
            const p = anchor.mult(shr).mult(clock_hand);
            qan.write(Color.init(0, 1, 1), @floatToInt(u32, p.x()), @floatToInt(u32, p.y()));
            clock_hand = rot.mult(clock_hand);
        }
    }

    var encoded = try imgio.encodeQanvasAsQoi(qan, gpa);
    defer gpa.free(encoded);

    try imgio.saveBufAsFile(encoded, "out.qoi");

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

    var sdl_tex = try imgio.encodeQanvasAsSdl(qan, gpa, sdl_renderer);
    defer sdl_tex.destroy();

    main_loop: while (true) {
        while (sdl2.pollEvent()) |ev| {
            switch (ev) {
                .quit => break :main_loop,
                else => {},
            }
        }

        try sdl_renderer.copy(sdl_tex, null, null);
        sdl_renderer.present();
        sdl2.delay(7); // ~144hz
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
}
