const std = @import("std");
const sdl2 = @import("sdl2");

const imgio = @import("imgio.zig");
const trans = @import("transform.zig");
const rndr = @import("render.zig");
const vol = @import("volume.zig");
const mat = @import("material.zig");
const file = @import("world_file.zig");

const World = @import("world.zig").World;
const Ray = @import("ray.zig").Ray;
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

    const filename = blk: {
        if (std.os.argv.len > 1) {
            break :blk std.mem.span(std.os.argv[1]);
        }

        const default_scene_file = "scenes/demo_scene.txt";
        std.debug.print("No scene file specified, using {s}.\n", .{default_scene_file});
        break :blk default_scene_file;
    };
    const scene_file = file.parseWorldFile(filename, gpa) catch |e| switch (e) {
        error.FileNotFound => {
            std.debug.print("File not found: \"{s}\".\n", .{filename});
            return;
        },
        else => return e,
    };

    defer gpa.free(scene_file.text);
    defer gpa.free(scene_file.sections);
    defer scene_file.world.deinit();

    var qan = try Qanvas.init(
        gpa,
        @intCast(usize, scene_file.camera.width),
        @intCast(usize, scene_file.camera.height),
    );
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

    {
        var thr = std.Thread.spawn(
            .{},
            World.render,
            .{ scene_file.world, scene_file.camera, &qan, gpa },
        ) catch unreachable;
        thr.detach();

        main_loop: while (true) {
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
    _ = @import("volume.zig");
    _ = @import("world.zig");
    _ = @import("world_file.zig");
    _ = @import("intersect.zig");
    _ = @import("light.zig");
    _ = @import("material.zig");
    _ = @import("prefab.zig");
}
