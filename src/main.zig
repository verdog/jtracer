const std = @import("std");

const Qanvas = @import("qanvas.zig").Qanvas;
const Color = @import("color.zig").Color;
const imgio = @import("imgio.zig");

var gpa_impl = std.heap.GeneralPurposeAllocator(.{}){};
const gpa = gpa_impl.allocator();

pub fn main() !void {
    defer if (!gpa_impl.detectLeaks()) std.debug.print("(No leaks)\n", .{});

    var qan = try Qanvas.init(gpa, 256, 256);
    defer qan.deinit();

    var i: usize = 0;
    while (i < qan.width and i < qan.height) : (i += 1) {
        qan.write(Color.init(0, 1, 0), i, i);
        qan.write(Color.init(1, 1, 0), i, qan.height - i - 1);
    }

    var encoded = try imgio.encodeQanvasAsQoi(qan, gpa);
    defer gpa.free(encoded);

    try imgio.saveBufAsFile(encoded, "out.qoi");
}

test {
    _ = @import("tuple.zig");
    _ = @import("color.zig");
    _ = @import("qanvas.zig");
    _ = @import("qoi.zig");
}
