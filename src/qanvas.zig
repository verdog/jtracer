//! qanvas is a collection of data in memory representing an image.
//! it supports some basic drawing operations.
//! when written to disk, it is written in the form of a qoi image.

const std = @import("std");
const Color = @import("color.zig").Color;

pub const Qanvas = struct {
    width: usize,
    height: usize,
    pixels: []Color,
    alctr: std.mem.Allocator,

    pub fn init(alctr: std.mem.Allocator, width: usize, height: usize) !This {
        var q = Qanvas{
            .width = width,
            .height = height,
            .pixels = try alctr.alloc(Color, width * height),
            .alctr = alctr,
        };

        for (q.pixels) |*qix, i| {
            const xy = q.ithPixelCoords(i);
            if (@divTrunc(xy.x, 64) & 1 == @divTrunc(xy.y, 64) & 1) {
                // magenta
                qix.* = Color.init(0.1, 0.1, 0.1);
            } else {
                // black
                qix.* = Color.init(0, 0, 0);
            }
        }

        return q;
    }

    pub fn deinit(self: *This) void {
        self.alctr.free(self.pixels);
    }

    pub fn ithPixelCoords(self: This, i: usize) struct { x: usize, y: usize } {
        return .{
            .x = i % self.width,
            .y = @divTrunc(i, self.width),
        };
    }

    pub fn at(self: This, x: usize, y: usize) Color {
        return self.pixels[y * self.width + x];
    }

    pub fn write(self: This, color: Color, x: usize, y: usize) void {
        self.pixels[y * self.width + x] = color;
    }

    const This = @This();
};

const expect = std.testing.expect;

test "Construct/Destruct Qanvas" {
    var alloc = std.testing.allocator;
    var q = try Qanvas.init(alloc, 1024, 512);
    defer q.deinit();

    try expect(q.width == 1024);
    try expect(q.height == 512);
}

test "Write pixel" {
    var alloc = std.testing.allocator;
    var q = try Qanvas.init(alloc, 64, 64);
    defer q.deinit();

    const red = Color.init(1, 0, 0);

    q.write(red, 24, 24);
    try expect(q.at(24, 24).equals(red));
}

test "ithPixelCoords" {
    var alloc = std.testing.allocator;
    var q = try Qanvas.init(alloc, 32, 8);
    defer q.deinit();

    try expect(q.ithPixelCoords(0).x == 0);
    try expect(q.ithPixelCoords(0).y == 0);

    try expect(q.ithPixelCoords(31).x == 31);
    try expect(q.ithPixelCoords(31).y == 0);

    try expect(q.ithPixelCoords(63).x == 31);
    try expect(q.ithPixelCoords(63).y == 1);
}
