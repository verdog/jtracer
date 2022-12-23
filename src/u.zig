//! utilities

const std = @import("std");

pub fn print(it: anytype) void {
    switch (@TypeOf(it)) {
        @import("tuple.zig").Tuple => {
            std.debug.print("({d:.5}, {d:.5}, {d:.5}, {d:.5})\n", .{ it.vec[0], it.vec[1], it.vec[2], it.vec[3] });
        },
        @import("matrix.zig").Matrix(4, 4) => {
            const fmt =
                \\[
                \\  {d: <8.5} {d: <8.5} {d: <8.5} {d: <8.5}
                \\  {d: <8.5} {d: <8.5} {d: <8.5} {d: <8.5}
                \\  {d: <8.5} {d: <8.5} {d: <8.5} {d: <8.5}
                \\  {d: <8.5} {d: <8.5} {d: <8.5} {d: <8.5}
                \\]
                \\
            ;
            std.debug.print(fmt, .{
                it.vec[0],  it.vec[1],  it.vec[2],  it.vec[3],
                it.vec[4],  it.vec[5],  it.vec[6],  it.vec[7],
                it.vec[8],  it.vec[9],  it.vec[10], it.vec[11],
                it.vec[12], it.vec[13], it.vec[14], it.vec[15],
            });
        },
        @import("color.zig").Color => {
            std.debug.print("({d:.5}, {d:.5}, {d:.5})\n", .{
                it.vec[0],
                it.vec[1],
                it.vec[2],
            });
        },
        f64 => {
            std.debug.print("{d: <8.5}\n", .{it});
        },
        else => {
            std.debug.print("{}\n", .{it});
        },
    }
}
