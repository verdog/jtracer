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
    background_color: Color,

    pub fn init(alctr: std.mem.Allocator, width: usize, height: usize) !This {
        var q = Qanvas{
            .width = width,
            .height = height,
            .pixels = try alctr.alloc(Color, width * height),
            .background_color = Color.init(0.5, 0.5, 0.5),
            .alctr = alctr,
        };

        q.clear();

        return q;
    }

    pub fn clear(self: *This) void {
        for (self.pixels) |*qix| {
            qix.* = self.background_color;
        }
    }

    pub fn deinit(self: *This) void {
        self.alctr.free(self.pixels);
    }

    pub fn ithPixelCoords(self: This, i: usize) struct { x: i64, y: i64 } {
        return .{
            .x = @intCast(i64, i % self.width),
            .y = @intCast(i64, @divTrunc(i, self.width)),
        };
    }

    pub fn at(self: This, x: i64, y: i64) Color {
        const idx = y * @intCast(i64, self.width) + x;
        if (idx < 0 or idx >= @intCast(i64, self.pixels.len)) return Color.init(0, 0, 0);

        return self.pixels[@intCast(usize, idx)];
    }

    pub fn write(self: This, color: Color, x: i64, y: i64) void {
        const idx = y * @intCast(i64, self.width) + x;
        if (idx < 0 or idx >= @intCast(i64, self.pixels.len)) return;

        self.pixels[@intCast(usize, idx)] = color;
    }

    pub fn fill(self: This, x: i64, y: i64, width: i64, height: i64, color: Color) void {
        if (width <= 0 or height <= 0) return;

        const start_x: u64 = @intCast(u64, @max(0, x));
        const start_y: u64 = @intCast(u64, @max(0, y));
        const end_x: u64 = @intCast(u64, @min(@intCast(i64, self.width) - 1, x + width - 1));
        const end_y: u64 = @intCast(u64, @min(@intCast(i64, self.height) - 1, y + height - 1));

        for (start_y..end_y + 1) |y_pix| {
            for (start_x..end_x + 1) |x_pix| {
                self.write(color, @intCast(i64, x_pix), @intCast(i64, y_pix));
            }
        }
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

test "Write pixel out of bounds" {
    var alloc = std.testing.allocator;
    var q = try Qanvas.init(alloc, 64, 64);
    defer q.deinit();

    const red = Color.init(1, 0, 0);

    q.write(red, 128, 128);
    q.write(red, -128, -128);
}

test "Read out of bounds returns black" {
    var alloc = std.testing.allocator;
    var q = try Qanvas.init(alloc, 64, 64);
    defer q.deinit();

    const black = Color.init(0, 0, 0);

    try expect(q.at(-1, -1).equals(black));
    try expect(q.at(64, 64).equals(black));
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

test "fill" {
    var alloc = std.testing.allocator;
    var q = try Qanvas.init(alloc, 64, 64);
    defer q.deinit();

    const blue = Color.init(0, 0, 1);

    q.fill(0, 0, 32, 32, blue);

    try expect(q.at(0, 0).equals(blue));
    try expect(q.at(31, 31).equals(blue));
    try expect(q.at(63, 63).equals(q.background_color));
}

test "fill oob" {
    var alloc = std.testing.allocator;
    var q = try Qanvas.init(alloc, 64, 64);
    defer q.deinit();

    const blue = Color.init(0, 0, 1);

    q.fill(-16, -16, 32, 32, blue);

    try expect(q.at(0, 0).equals(blue));
    try expect(q.at(15, 15).equals(blue));
    try expect(q.at(31, 31).equals(q.background_color));
    try expect(q.at(63, 63).equals(q.background_color));
}
