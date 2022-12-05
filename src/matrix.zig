//! simple matrix implementation

const std = @import("std");
const mymath = @import("mymath.zig");

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

        pub fn fromVector(inits: @Vector(compt_rows * compt_columns, f64)) This {
            return .{ .vec = inits };
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
            return max_diff <= mymath.floatTolerance;
        }

        pub fn equalsTolerance(self: This, other: This, tolerance: usize) bool {
            const diff = @fabs(self.vec - other.vec);
            const max_diff = @reduce(.Max, diff);
            return max_diff <= std.math.floatEps(f64) * @intToFloat(f64, tolerance);
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

        pub fn mult(self: This, other: anytype) @TypeOf(other) {
            switch (@TypeOf(other)) {
                This => { // matrix x matrix
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
                },
                Tuple => { // matrix x tuple
                    if (This.compt_columns != 4)
                        @compileError("mult matrix must have 4 columns");

                    const myRows = self.sliceRows();

                    return Tuple.init(
                        @reduce(.Add, myRows[0] * other.vec),
                        @reduce(.Add, myRows[1] * other.vec),
                        @reduce(.Add, myRows[2] * other.vec),
                        @reduce(.Add, myRows[3] * other.vec),
                    );
                },
                else => @compileError("Matrix.mult can only be used on Matrices or Tuples of compatible sizes"),
            }
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

        pub fn determinant(self: This) f64 {
            if (This.compt_rows != This.compt_columns)
                @compileError("Only square matrices are supported");

            switch (This.compt_rows) {
                2 => return self.vec[0] * self.vec[3] - self.vec[1] * self.vec[2],
                3...4 => {
                    var elemVec = @splat(This.compt_columns, self.vec[0]);
                    var cofVec = @splat(This.compt_columns, self.cofactor(0, 0));

                    comptime var i: usize = 1;
                    inline while (i < This.compt_columns) : (i += 1) {
                        elemVec[i] = self.vec[i];
                        cofVec[i] = self.cofactor(0, i);
                    }

                    return @reduce(.Add, elemVec * cofVec);
                },
                else => @compileError("Unsupported size"),
            }
        }

        pub fn isInvertible(self: This) bool {
            // TODO use epsilon
            return self.determinant() != 0;
        }

        pub fn inverted(self: This) !This {
            if (!self.isInvertible()) return error.NotInvertible;

            var cofacts = blk: {
                var result = self.vec;

                comptime var i: usize = 0;
                inline while (i < This.compt_rows * This.compt_columns) : (i += 1) {
                    const this_row = i / This.compt_rows;
                    const this_col = i % This.compt_columns;
                    result[i] = self.cofactor(this_row, this_col);
                }

                break :blk result;
            };

            const size = This.compt_rows * This.compt_columns;
            cofacts /= @splat(size, self.determinant());

            return This.fromVector(cofacts).transposed();
        }

        pub fn submatrix(self: This, row: usize, column: usize) Matrix(This.compt_rows - 1, This.compt_columns - 1) {
            // asumming a square matrix
            var inits: std.meta.Tuple(&[_]type{f64} ** (_rows - 1) ** (_cols - 1)) = .{0} ** (_rows - 1) ** (_cols - 1);

            var self_i: usize = 0;
            inline for (inits) |*val| {
                while ((self_i / This.compt_rows) == row or (self_i % This.compt_columns) == column) self_i += 1;
                val.* = self.vec[self_i];
                self_i += 1;
            }

            return Matrix(This.compt_rows - 1, This.compt_columns - 1).init(inits);
        }

        pub fn minor(self: This, row: usize, column: usize) f64 {
            return self.submatrix(row, column).determinant();
        }

        pub fn cofactor(self: This, row: usize, column: usize) f64 {
            const min = self.minor(row, column);
            if ((row ^ column) & 1 == 1) {
                return min * -1.0;
            }
            return min;
        }

        const This = @This();
    };
}

const expect = std.testing.expect;
const expectEq = std.testing.expectEqual;

test "Represent a 4x4 matrix" {
    const m = Matrix(4, 4).init(.{
        1,    2,    3,    4,
        5.5,  6.5,  7.5,  8.5,
        9,    10,   11,   12,
        13.5, 14.5, 15.5, 16.5,
    });

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
    const m = Matrix(3, 3).init(.{
        1,   2,   3,
        5.5, 6.5, 7.5,
        9,   10,  11,
    });

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
    const m = Matrix(2, 2).init(.{
        1,   2,
        5.5, 6.5,
    });

    try expectEq(m.at(0, 0), 1);
    try expectEq(m.at(0, 1), 2);

    try expectEq(m.at(1, 0), 5.5);
    try expectEq(m.at(1, 1), 6.5);
}

test "Matrix equals" {
    const m = Matrix(4, 4).init(.{
        1,    2,    3,    4,
        5.5,  6.5,  7.5,  8.5,
        9,    10,   11,   12,
        13.5, 14.5, 15.5, 16.5,
    });

    try expect(m.equals(m));
}

test "Matrix !equals" {
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

    try expect(!m.equals(m2));
    try expect(!m2.equals(m));
}

test "matrix multiplication" {
    const m = Matrix(4, 4).init(.{
        1, 2, 3, 4,
        5, 6, 7, 8,
        9, 8, 7, 6,
        5, 4, 3, 2,
    });
    const m2 = Matrix(4, 4).init(.{
        -2, 1, 2, 3,
        3,  2, 1, -1,
        4,  3, 6, 5,
        1,  2, 7, 8,
    });
    const m3 = Matrix(4, 4).init(.{
        20, 22, 50,  48,
        44, 54, 114, 108,
        40, 58, 110, 102,
        16, 26, 46,  42,
    });

    try expect(m.mult(m2).equals(m3));
}

test "matrix multiply by tuple" {
    const m = Matrix(4, 4).init(.{
        1, 2, 3, 4,
        2, 4, 4, 2,
        8, 6, 4, 1,
        0, 0, 0, 1,
    });

    const t = Tuple.init(1, 2, 3, 1);
    const a = Tuple.init(18, 24, 33, 1);

    try expect(m.mult(t).equals(a));
}

test "matrix transpose" {
    const m = Matrix(4, 4).init(.{
        1, 2, 3, 4,
        5, 6, 7, 8,
        9, 8, 7, 6,
        5, 4, 3, 2,
    });
    const t = Matrix(4, 4).init(.{
        1, 5, 9, 5,
        2, 6, 8, 4,
        3, 7, 7, 3,
        4, 8, 6, 2,
    });

    try expect(m.transposed().equals(t));
}

test "matrix identity" {
    const m = Matrix(4, 4).init(.{
        1, 0, 0, 0,
        0, 1, 0, 0,
        0, 0, 1, 0,
        0, 0, 0, 1,
    });

    try expect(m.equals(Matrix(4, 4).identity()));
}

test "matrix transpose of identity" {
    const m = Matrix(4, 4).identity();

    try expect(m.transposed().equals(m));
}

test "matrix multiplication of identity" {
    const m = Matrix(4, 4).init(.{
        1, 2, 3, 4,
        5, 6, 7, 8,
        9, 8, 7, 6,
        5, 4, 3, 2,
    });

    try expect(m.mult(Matrix(4, 4).identity()).equals(m));
}

test "2x2 matrix determinant" {
    const m = Matrix(2, 2).init(.{ 1, 5, -3, 2 });

    try expectEq(m.determinant(), 17);
}

test "3x3 matrix determinant" {
    const m = Matrix(3, 3).init(.{
        1,  2, 6,
        -5, 8, -4,
        2,  6, 4,
    });

    try expectEq(m.cofactor(0, 0), 56);
    try expectEq(m.cofactor(0, 1), 12);
    try expectEq(m.cofactor(0, 2), -46);
    try expectEq(m.determinant(), -196);
}

test "4x4 matrix determinant" {
    const m = Matrix(4, 4).init(.{
        -2, -8, 3,  5,
        -3, 1,  7,  3,
        1,  2,  -9, 6,
        -6, 7,  7,  -9,
    });

    try expectEq(m.cofactor(0, 0), 690);
    try expectEq(m.cofactor(0, 1), 447);
    try expectEq(m.cofactor(0, 2), 210);
    try expectEq(m.cofactor(0, 3), 51);
    try expectEq(m.determinant(), -4071);
}

test "A submatrix of a 4x4 matrix is a 3x3 matrix" {
    const m = Matrix(4, 4).init(.{
        1, 2, 3, 4,
        5, 6, 7, 8,
        9, 8, 7, 6,
        5, 4, 3, 2,
    });

    const a = Matrix(3, 3).init(.{
        1, 3, 4,
        5, 7, 8,
        5, 3, 2,
    });

    const sub = m.submatrix(2, 1);

    try expect(sub.equals(a));
}

test "A submatrix of a 3x3 matrix is a 2x2 matrix" {
    const m = Matrix(3, 3).init(.{
        1, 3, 4,
        5, 7, 8,
        5, 3, 2,
    });

    const a = Matrix(2, 2).init(.{
        1, 4,
        5, 2,
    });

    const sub = m.submatrix(1, 1);

    try expect(sub.equals(a));
}

test "3x3 matrix minor" {
    const m = Matrix(3, 3).init(.{
        3, 5,  0,
        2, -1, -7,
        6, -1, 5,
    });

    const sub = m.submatrix(1, 0);

    try expect(sub.determinant() == 25);
    try expect(m.minor(1, 0) == 25);
}

test "3x3 matrix cofactor" {
    const m = Matrix(3, 3).init(.{
        3, 5,  0,
        2, -1, -7,
        6, -1, 5,
    });

    try expect(m.minor(0, 0) == -12);
    try expect(m.cofactor(0, 0) == -12);
    try expect(m.minor(1, 0) == 25);
    try expect(m.cofactor(1, 0) == -25);
}

test "check for invertibility: invertible" {
    const m = Matrix(4, 4).init(.{
        6, 4,  4, 4,
        5, 5,  7, 6,
        4, -9, 3, -7,
        9, 1,  7, -6,
    });

    try expectEq(m.determinant(), -2120);
    try expect(m.isInvertible());
}

test "check for invertibility: not invertible" {
    const m = Matrix(4, 4).init(.{
        -4, 2,  -2, -3,
        9,  6,  2,  6,
        0,  -5, 1,  -5,
        0,  0,  0,  0,
    });

    try expectEq(m.determinant(), 0);
    try expect(!m.isInvertible());
}

test "inverting a matrix" {
    const m = Matrix(4, 4).init(.{
        -5, 2,  6,  -8,
        1,  -5, 1,  8,
        7,  7,  -6, -7,
        1,  -3, 7,  4,
    });

    const ii = Matrix(4, 4).init(.{
        0.21805,  0.45113,  0.24060,  -0.04511,
        -0.80827, -1.45677, -0.44361, 0.52068,
        -0.07895, -0.22368, -0.05263, 0.19737,
        -0.52256, -0.81391, -0.30075, 0.30639,
    });

    const i = try m.inverted();

    try expectEq(m.determinant(), 532.0);
    try expectEq(m.cofactor(2, 3), -160.0);
    try expectEq(i.at(3, 2), -160.0 / 532.0);
    try expectEq(m.cofactor(3, 2), 105.0);
    try expectEq(i.at(2, 3), 105.0 / 532.0);

    // book examples are much less precise than f64
    try expect(i.equalsTolerance(ii, 100_000_000_000));
}

test "inverting a matrix 2" {
    const m = Matrix(4, 4).init(.{
        8,  -5, 9,  2,
        7,  5,  6,  1,
        -6, 0,  9,  6,
        -3, 0,  -9, -4,
    });

    const ii = Matrix(4, 4).init(.{
        -0.15385, -0.15385, -0.28205, -0.53846,
        -0.07692, 0.12308,  0.02562,  0.03077,
        0.35897,  0.35897,  0.43590,  0.92308,
        -0.69231, -0.69231, -0.76923, -1.92308,
    });

    const i = try m.inverted();

    // book examples are much less precise than f64
    try expect(i.equalsTolerance(ii, 100_000_000_000));
}

test "inverting a matrix 3" {
    const m = Matrix(4, 4).init(.{
        9,  3,  0,  9,
        -5, -2, -6, -3,
        -4, 9,  6,  4,
        -7, 6,  6,  2,
    });

    const ii = Matrix(4, 4).init(.{
        -0.04074, -0.07778, 0.14444,  -0.22222,
        -0.07778, 0.03333,  0.36667,  -0.33333,
        -0.02901, -0.14630, -0.10926, 0.12963,
        0.17778,  0.06667,  -0.26667, 0.33333,
    });

    const i = try m.inverted();

    // book examples are much less precise than f64
    try expect(i.equalsTolerance(ii, 100_000_000_000));
}

test "matrix multiplied by its inverse is itself" {
    const m = Matrix(4, 4).init(.{
        1, 2, 3, 4,
        5, 6, 7, 8,
        9, 8, 7, 6,
        5, 4, 3, 2,
    });
    const m2 = Matrix(4, 4).init(.{
        -2, 1, 2, 3,
        3,  2, 1, -1,
        4,  3, 6, 5,
        1,  2, 7, 8,
    });

    const m3 = m.mult(m2);
    const inverted = m3.mult(try m2.inverted());

    try expect(inverted.equalsTolerance(m, 512));
}
