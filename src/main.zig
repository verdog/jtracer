const std = @import("std");
const sdl2 = @import("sdl2");

const imgio = @import("imgio.zig");
const trans = @import("transform.zig");
const rndr = @import("render.zig");
const vol = @import("volume.zig");
const mat = @import("material.zig");

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

    const scale = 1;
    var qan = try Qanvas.init(gpa, 1200 * scale, 1200 * scale);
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
        lgt.ptr.position = Point.init(-4, 10, -5);
        lgt.ptr.intensity = Color.init(1, 1, 1);

        {
            // checkered cube
            var cube = world.addVolume(vol.Cube);
            cube.ptr.material.color_map = mat.ThreeDCheckedColor.init(
                Color.init(0.6, 0.2, 0.2),
                Color.init(0.2, 0.2, 0.6),
            );

            cube.ptr.transform = cube.ptr.transform.chain(.{
                trans.makeTranslation(3, 2, 0),
                trans.makeScaling(0.7, 2, 0.7),
                trans.makeRotationY(std.math.pi / -12.0),
            });
        }

        {
            // glass rectangle
            var cube = world.addVolume(vol.Cube);
            cube.ptr.material.color_map = mat.FlatColor.init(
                Color.init(0.2, 0.2, 0.2),
            );
            cube.ptr.material.diffuse = 0.2;
            cube.ptr.material.ambient = 0.1;
            cube.ptr.material.reflective = 0.9;
            cube.ptr.material.transparency = 0.9;
            cube.ptr.material.refractive_index = 1.52;

            cube.ptr.transform = cube.ptr.transform.chain(.{
                trans.makeTranslation(0, 0.7001, -3),
                trans.makeRotationY(0.3),
                trans.makeRotationX(std.math.pi / 2.0),
                trans.makeScaling(3, 0.1, 0.7),
            });
        }

        {
            // checkered cube
            var cube = world.addVolume(vol.Cube);
            cube.ptr.material.color_map = mat.ThreeDCheckedColor.init(
                Color.init(0.6, 0.2, 0.2),
                Color.init(0.2, 0.2, 0.6),
            );
            cube.ptr.material.diffuse = 0.1;
            cube.ptr.material.ambient = 0;
            cube.ptr.material.reflective = 1;

            cube.ptr.transform = cube.ptr.transform.chain(.{
                trans.makeTranslation(0, 1, 0),
                trans.makeRotationY(std.math.pi / 4.0 + 0.2),
            });
        }

        {
            var sph = world.addVolume(vol.Sphere);
            sph.ptr.material.color_map = mat.FlatColor.init(Color.init(0.2, 0.2, 0.7));

            sph.ptr.transform = sph.ptr.transform.chain(.{
                trans.makeTranslation(0, 3, 0),
            });
        }

        {
            // right sphere
            var sph = world.addVolume(vol.Sphere);
            sph.ptr.material.color_map = mat.FlatColor.init(Color.init(0.2, 0.2, 0.7));

            sph.ptr.transform = sph.ptr.transform.chain(.{
                trans.makeTranslation(3, 1, -2.5),
            });
        }

        {
            // checkered cube
            var cube = world.addVolume(vol.Cube);
            cube.ptr.material.color_map = mat.ThreeDCheckedColor.init(
                Color.init(0.6, 0.2, 0.2),
                Color.init(0.2, 0.2, 0.6),
            );

            cube.ptr.transform = cube.ptr.transform.chain(.{
                trans.makeTranslation(-3, 0.7, 0),
                trans.makeScaling(1.2, 0.7, 1.2),
                // trans.makeRotationY(std.math.pi / -12.0),
            });
        }

        {
            // small cube on left cube
            var cube = world.addVolume(vol.Cube);
            cube.ptr.material.color_map = mat.RingColor.init(
                Color.init(0.6, 0.2, 0.2),
                Color.init(0.2, 0.2, 0.6),
            );

            cube.ptr.transform = cube.ptr.transform.chain(.{
                trans.makeTranslation(-3.3, 1.7, 0.3),
                trans.makeRotationY(42),
                trans.makeScaling(0.3, 0.3, 0.3),
            });
        }

        {
            // small cylinder on left cube
            var cyl = world.addVolume(vol.Cylinder);
            cyl.ptr.material.color_map = mat.ThreeDCheckedColor.init(
                Color.init(0.6, 0.2, 0.2),
                Color.init(0.2, 0.2, 0.6),
            );

            cyl.ptr.closed = true;
            cyl.ptr.length = 1.5;

            cyl.ptr.transform = cyl.ptr.transform.chain(.{
                trans.makeTranslation(-2.2, 2.15, -0.8),
                trans.makeRotationY(32),
                trans.makeScaling(0.3, 1, 0.3),
            });
        }

        {
            // tall cylinder on left cube
            var cyl = world.addVolume(vol.Cylinder);
            cyl.ptr.material.color_map = mat.FlatColor.init(
                Color.init(0.1, 0.1, 0.1),
            );

            cyl.ptr.material.diffuse = 0.2;
            cyl.ptr.material.ambient = 0.1;
            cyl.ptr.material.reflective = 0.9;
            cyl.ptr.material.transparency = 0.9;
            cyl.ptr.material.refractive_index = 1.52;

            cyl.ptr.closed = true;
            cyl.ptr.length = 2;

            cyl.ptr.transform = cyl.ptr.transform.chain(.{
                trans.makeTranslation(-3, 1.7, -0.5),
                trans.makeRotationY(std.math.pi / 4.0),
                trans.makeRotationX(std.math.pi / 2.0),
                trans.makeScaling(0.25, 2, 0.25),
            });
        }

        {
            // upright cylinder
            var cyl = world.addVolume(vol.Cylinder);
            cyl.ptr.material.color_map = mat.ThreeDCheckedColor.init(
                Color.init(0.6, 0.2, 0.2),
                Color.init(0.2, 0.2, 0.6),
            );

            cyl.ptr.length = 1;

            cyl.ptr.transform = cyl.ptr.transform.chain(.{
                trans.makeTranslation(-3.5, 0.8, -3.2),
                trans.makeRotationY(-1.1),
                trans.makeRotationX(std.math.pi / 2.0),
                trans.makeScaling(0.8, 0.4, 0.8),
            });
        }

        {
            var cone = world.addVolume(vol.Cone);
            cone.ptr.material.color_map = mat.FlatColor.init(
                Color.init(0.2, 0.2, 0.7),
            );

            cone.ptr.max = 1;
            cone.ptr.min = 0;
            cone.ptr.closed = false;

            cone.ptr.transform = cone.ptr.transform.chain(.{
                trans.makeTranslation(3, 4, 0),
            });
        }

        {
            // checkered floor
            var pln = world.addVolume(vol.Plane);
            pln.ptr.material.color_map = mat.ThreeDCheckedColor.initSingle(Color.init(
                0.5,
                0.5,
                0.5,
            ));
            pln.ptr.material.specular = 0.1;
            pln.ptr.material.diffuse = 1;
            pln.ptr.material.reflective = 0.12;
            pln.ptr.material.transform = pln.ptr.material.transform.chain(.{
                trans.makeScaling(8, 8, 8),
                trans.makeTranslation(0.5, 0, 0.5),
            });
        }

        var cam = Camera.init(
            @intCast(i64, qan.width),
            @intCast(i64, qan.height),
            std.math.pi / 2.5,
        );

        // slightly down
        const from = Point.init(-0.2, 3, -10);
        const to = Point.init(-0.2, 2, 0);
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
    _ = @import("volume.zig");
    _ = @import("world.zig");
    _ = @import("intersect.zig");
    _ = @import("light.zig");
    _ = @import("material.zig");
    _ = @import("prefab.zig");
}
