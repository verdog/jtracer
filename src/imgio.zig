//! functions to save qanvas objects to disk as qoi images

const std = @import("std");
const sdl2 = @import("sdl2");

const Qanvas = @import("qanvas.zig").Qanvas;
const qoi = @import("qoi.zig");

pub fn encodeQanvasAsQoi(qan: Qanvas, alctr: std.mem.Allocator) ![]u8 {
    var qixels = try alctr.alloc(qoi.Qixel, qan.pixels.len);
    defer alctr.free(qixels);

    for (qixels) |*q, i| {
        q.* = qoi.Qixel.fromColor(qan.pixels[i]);
    }

    return try qoi.encode(qixels, alctr, @intCast(u32, qan.width), @intCast(u32, qan.height), .rgb, .all_linear);
}

pub fn encodeQanvasAsSdl(qan: Qanvas, alctr: std.mem.Allocator, renderer: sdl2.Renderer) !sdl2.Texture {
    var qixels = try alctr.alloc(qoi.Qixel, qan.pixels.len);
    defer alctr.free(qixels);

    for (qixels) |*q, i| {
        q.* = qoi.Qixel.fromColor(qan.pixels[i]);
    }

    var tex = try sdl2.createTexture(renderer, .abgr8888, .static, qan.width, qan.height);
    errdefer tex.destroy();

    try tex.update(@ptrCast(*[]u8, &qixels).*, qan.width * @sizeOf(qoi.Qixel), null);

    return tex;
}

pub fn saveBufAsFile(buf: []const u8, filename: []const u8) !void {
    var cwd = std.fs.cwd();

    var file = try cwd.createFile(filename, .{ .truncate = true });
    defer file.close();

    _ = try file.write(buf);
}
