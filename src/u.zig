//! utilities

const std = @import("std");

pub fn print(it: anytype) void {
    switch (@TypeOf(it)) {
        @import("tuple.zig").Tuple => {
            std.debug.print("({d:.5}, {d:.5}, {d:.5}, {d:.5})\n", .{ it.vec[0], it.vec[1], it.vec[2], it.vec[3] });
        },
        else => {
            std.debug.print("{}\n", .{it});
        },
    }
}
