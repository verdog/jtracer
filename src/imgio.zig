//! functions to change one image format to another

const std = @import("std");
const sdl2 = @import("sdl2");

const Qanvas = @import("qanvas.zig").Qanvas;
const qoi = @import("qoi.zig");

pub fn encodeQanvasAsQoi(qan: Qanvas, alctr: std.mem.Allocator) ![]u8 {
    var qixels = try alctr.alloc(qoi.Qixel, qan.pixels.len);
    defer alctr.free(qixels);

    for (qixels, 0..) |*q, i| {
        q.* = qoi.Qixel.fromColor(qan.pixels[i]);
    }

    return try qoi.encode(qixels, alctr, @intCast(qan.width), @intCast(qan.height), .rgb, .all_linear);
}

pub fn encodeQanvasAsQixel(qan: Qanvas, alctr: std.mem.Allocator) ![]qoi.Qixel {
    var qixels = try alctr.alloc(qoi.Qixel, qan.pixels.len);
    errdefer alctr.free(qixels);

    try encodeQanvasAsQixelUpdate(qan, qixels);

    return qixels;
}

pub fn encodeQanvasAsQixelUpdate(qan: Qanvas, to_update: []qoi.Qixel) !void {
    for (to_update, 0..) |*q, i| {
        q.* = qoi.Qixel.fromColor(qan.pixels[i]);
    }
}

pub fn encodeQixelsAsSdl(qix: []qoi.Qixel, renderer: sdl2.Renderer, width: u32, height: u32) !sdl2.Texture {
    var tex = try sdl2.createTexture(renderer, .abgr8888, .static, width, height);
    errdefer tex.destroy();

    try encodeQixelsAsSdlUpdate(qix, &tex, width);

    return tex;
}

pub fn encodeQixelsAsSdlUpdate(qix: []qoi.Qixel, tex: *sdl2.Texture, width: u32) !void {
    try tex.update(@as(*const []u8, @ptrCast(&qix)).*, width * @sizeOf(qoi.Qixel), null);
}

pub fn saveBufAsFile(buf: []const u8, filename: []const u8) !void {
    var cwd = std.fs.cwd();

    var file = try cwd.createFile(filename, .{ .truncate = true });
    defer file.close();

    _ = try file.write(buf);
}
