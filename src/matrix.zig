//! simple matrix implementation

const std = @import("std");

const Tuple = @import("tuple.zig").Tuple;

pub fn Matrix(comptime _rows: usize, comptime _cols: usize) type {
    if (_rows != _cols)
        @compileError("Only square matrices are supported");

    return struct {
        const compt_columns = _cols;
        const compt_rows = _rows;

        /// memory representation of the matrix
        vec: @Vector(compt_rows * compt_columns, f64),

        pub fn init(inits: std.meta.Tuple(&[_]type{f64} ** (_rows * _cols))) This {
            const zero: f64 = 0;
            var m = This{ .vec = @splat(_rows * _cols, zero) };

            inline for (inits) |val, i| {
                m.vec[i] = val;
            }

            return m;
        }

        pub fn identity() This {
            // asumming a square matrix
            const inits = comptime blk: {
                var inits: std.meta.Tuple(&[_]type{f64} ** (_rows * _cols)) = .{0} ** _rows ** _cols;
                var i: usize = 0;
                while (i < _rows) : (i += 1) {
                    @field(inits, std.fmt.comptimePrint("{}", .{i * _cols + i})) = 1;
                }
                break :blk inits;
            };

            return This.init(inits);
        }

        pub fn rows(_: This) usize {
            return This.compt_rows;
        }

        pub fn columns(_: This) usize {
            return This.compt_columns;
        }

        pub fn at(self: This, row: usize, col: usize) f64 {
            return self.vec[row * self.columns() + col];
        }

        pub fn equals(self: This, other: This) bool {
            const diff = @fabs(self.vec - other.vec);
            const max_diff = @reduce(.Max, diff);
            return max_diff <= std.math.floatEps(f64);
        }

        fn sliceRows(self: This) [This.compt_rows]@Vector(This.compt_columns, f64) {
            const zero: f64 = 0;
            var result: [This.compt_rows]@Vector(This.compt_columns, f64) = .{@splat(This.compt_columns, zero)} ** This.compt_rows;

            comptime var i: usize = 0;
            inline while (i < This.compt_columns * This.compt_rows) : (i += 1) {
                result[i / This.compt_rows][i % This.compt_columns] = self.vec[i];
            }

            return result;
        }

        pub fn mult(self: This, other: This) This {
            const othert = other.transposed();

            const selfRows = self.sliceRows();
            const otherRows = othert.sliceRows();

            // asumming a square matrix
            var inits: std.meta.Tuple(&[_]type{f64} ** (_rows * _cols)) = .{0} ** _rows ** _cols;

            {
                comptime var r: usize = 0;
                inline while (r < _rows) : (r += 1) {
                    comptime var c: usize = 0;
                    inline while (c < _cols) : (c += 1) {
                        // zig fmt: off
                        @field(inits, std.fmt.comptimePrint("{}", .{r * _cols + c})) =
                            @reduce(.Add, selfRows[r] * otherRows[c]);
                        // zig fmt: on
                    }
                }
            }

            return This.init(inits);
        }

        pub fn multTuple(self: This, tuple: Tuple) Tuple {
            if (This.compt_columns != 4)
                @compileError("multTuple matrix must have 4 columns");

            const myRows = self.sliceRows();

            return Tuple.init(
                @reduce(.Add, myRows[0] * tuple.vec),
                @reduce(.Add, myRows[1] * tuple.vec),
                @reduce(.Add, myRows[2] * tuple.vec),
                @reduce(.Add, myRows[3] * tuple.vec),
            );
        }

        pub fn transposed(self: This) This {
            // asumming a square matrix
            var inits: std.meta.Tuple(&[_]type{f64} ** (_rows * _cols)) = .{0} ** _rows ** _cols;

            {
                comptime var y: usize = 0;
                inline while (y < _rows) : (y += 1) {
                    comptime var x: usize = 0;
                    inline while (x < _cols) : (x += 1) {
                        @field(inits, std.fmt.comptimePrint("{}", .{x * _cols + y})) = self.at(y, x);
                    }
                }
            }

            return This.init(inits);
        }

        const This = @This();
    };
}

const expect = std.testing.expect;
const expectEq = std.testing.expectEqual;

test "Represent a 4x4 matrix" {
    // zig fmt: off
    const m = Matrix(4, 4).init(.{
        1,    2,    3,    4,
        5.5,  6.5,  7.5,  8.5,
        9,    10,   11,   12,
        13.5, 14.5, 15.5, 16.5,
    });
    // zig fmt: on

    try expectEq(m.at(0, 0), 1);
    try expectEq(m.at(0, 1), 2);
    try expectEq(m.at(0, 2), 3);
    try expectEq(m.at(0, 3), 4);

    try expectEq(m.at(1, 0), 5.5);
    try expectEq(m.at(1, 1), 6.5);
    try expectEq(m.at(1, 2), 7.5);
    try expectEq(m.at(1, 3), 8.5);

    try expectEq(m.at(2, 0), 9);
    try expectEq(m.at(2, 1), 10);
    try expectEq(m.at(2, 2), 11);
    try expectEq(m.at(2, 3), 12);

    try expectEq(m.at(3, 0), 13.5);
    try expectEq(m.at(3, 1), 14.5);
    try expectEq(m.at(3, 2), 15.5);
    try expectEq(m.at(3, 3), 16.5);
}

test "Represent a 3x3 matrix" {
    // zig fmt: off
    const m = Matrix(3, 3).init(.{
        1,   2,   3,
        5.5, 6.5, 7.5,
        9,   10,  11,
    });
    // zig fmt: on

    try expectEq(m.at(0, 0), 1);
    try expectEq(m.at(0, 1), 2);
    try expectEq(m.at(0, 2), 3);

    try expectEq(m.at(1, 0), 5.5);
    try expectEq(m.at(1, 1), 6.5);
    try expectEq(m.at(1, 2), 7.5);

    try expectEq(m.at(2, 0), 9);
    try expectEq(m.at(2, 1), 10);
    try expectEq(m.at(2, 2), 11);
}

test "Represent a 2x2 matrix" {
    // zig fmt: off
    const m = Matrix(2, 2).init(.{
        1,   2,
        5.5, 6.5,
    });
    // zig fmt: on

    try expectEq(m.at(0, 0), 1);
    try expectEq(m.at(0, 1), 2);

    try expectEq(m.at(1, 0), 5.5);
    try expectEq(m.at(1, 1), 6.5);
}

test "Matrix equals" {
    // zig fmt: off
    const m = Matrix(4, 4).init(.{
        1,    2,    3,    4,
        5.5,  6.5,  7.5,  8.5,
        9,    10,   11,   12,
        13.5, 14.5, 15.5, 16.5,
    });
    // zig fmt: on

    try expect(m.equals(m));
}

test "Matrix !equals" {
    // zig fmt: off
    const m = Matrix(4, 4).init(.{
        1,    2,    3,    4,
        5.5,  6.5,  7.5,  8.5,
        9,    10,   11,   12,
        13.5, 14.5, 15.5, 16.5,
    });
    const m2 = Matrix(4, 4).init(.{
        9,    10,   11,   12,
        13.5, 14.5, 15.5, 16.5,
        1,    2,    3,    4,
        5.5,  6.5,  7.5,  8.5,
    });
    // zig fmt: on

    try expect(!m.equals(m2));
    try expect(!m2.equals(m));
}

test "matrix multiplication" {
    // zig fmt: off
    const m = Matrix(4, 4).init(.{
        1,    2,    3,    4,
        5,    6,    7,    8,
        9,    8,    7,    6,
        5,    4,    3,    2,
    });
    const m2 = Matrix(4, 4).init(.{
       -2,    1,    2,    3,
        3,    2,    1,   -1,
        4,    3,    6,    5,
        1,    2,    7,    8,
    });
    const m3 = Matrix(4, 4).init(.{
        20,    22,    50,    48,
        44,    54,    114,   108,
        40,    58,    110,   102,
        16,    26,    46,    42,
    });
    // zig fmt: on

    try expect(m.mult(m2).equals(m3));
}

test "matrix multiply by tuple" {
    // zig fmt: off
    const m = Matrix(4, 4).init(.{
        1, 2, 3, 4,
        2, 4, 4, 2,
        8, 6, 4, 1,
        0, 0, 0, 1,
    });
    // zig fmt: on

    const t = Tuple.init(1, 2, 3, 1);
    const a = Tuple.init(18, 24, 33, 1);

    try expect(m.multTuple(t).equals(a));
}

test "matrix transpose" {
    // zig fmt: off
    const m = Matrix(4, 4).init(.{
        1,    2,    3,    4,
        5,    6,    7,    8,
        9,    8,    7,    6,
        5,    4,    3,    2,
    });
    const t = Matrix(4, 4).init(.{
        1,    5,    9,    5,
        2,    6,    8,    4,
        3,    7,    7,    3,
        4,    8,    6,    2,
    });
    // zig fmt: on

    try expect(m.transposed().equals(t));
}

test "matrix identity" {
    // zig fmt: off
    const m = Matrix(4, 4).init(.{
        1,    0,    0,    0,
        0,    1,    0,    0,
        0,    0,    1,    0,
        0,    0,    0,    1,
    });
    // zig fmt: on

    try expect(m.equals(Matrix(4, 4).identity()));
}

test "matrix transpose of identity" {
    const m = Matrix(4, 4).identity();

    try expect(m.transposed().equals(m));
}

test "matrix multiplication of identity" {
    // zig fmt: off
    const m = Matrix(4, 4).init(.{
        1,    2,    3,    4,
        5,    6,    7,    8,
        9,    8,    7,    6,
        5,    4,    3,    2,
    });
    // zig fmt: on

    try expect(m.mult(Matrix(4, 4).identity()).equals(m));
}
