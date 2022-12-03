//! functions to save qanvas objects to disk as qoi images

const std = @import("std");

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

pub fn saveBufAsFile(buf: []const u8, filename: []const u8) !void {
    var cwd = std.fs.cwd();

    var file = try cwd.createFile(filename, .{ .truncate = true });
    defer file.close();

    _ = try file.write(buf);
}
