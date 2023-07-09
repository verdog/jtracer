//! basic tuple type. a tuple is four elements long: x, y, z, w, where w is 1 if the tuple
//! represents a point and 0 if the tuple represents a vector.

const std = @import("std");

const mymath = @import("mymath.zig");

pub const Tuple = struct {
    pub fn init(xx: f64, yy: f64, zz: f64, ww: f64) This {
        return .{
            .vec = .{ xx, yy, zz, ww },
        };
    }

    pub fn equals(self: This, other: This) bool {
        const diff = @fabs(self.vec - other.vec);
        const max_diff = @reduce(.Max, diff);
        return max_diff <= mymath.floatTolerance;
    }

    pub fn equalsTolerance(self: This, other: This, tolerance: usize) bool {
        const diff = @fabs(self.vec - other.vec);
        const max_diff = @reduce(.Max, diff);
        return max_diff <= std.math.floatEps(f64) * @as(f64, @floatFromInt(tolerance));
    }

    pub fn plus(self: This, other: This) This {
        const t = This{ .vec = self.vec + other.vec };
        return t;
    }

    pub fn minus(self: This, other: This) This {
        const t = This{ .vec = self.vec - other.vec };
        return t;
    }

    pub fn negated(self: This) This {
        const zero: f64 = 0;
        return This{ .vec = @splat(4, zero) - self.vec };
    }

    pub fn scaled(self: This, f: f64) This {
        return This{ .vec = self.vec * @splat(4, f) };
    }

    pub fn divided(self: This, f: f64) This {
        return This{ .vec = self.vec / @splat(4, f) };
    }

    pub fn magnitude(self: This) f64 {
        return @sqrt(@reduce(.Add, self.vec * self.vec));
    }

    pub fn normalized(self: This) This {
        std.debug.assert(self.magnitude() != 0);
        var r = This{ .vec = self.vec / @splat(4, self.magnitude()) };
        r.vec[3] = self.vec[3];
        return r;
    }

    pub fn dot(self: This, other: This) f64 {
        std.debug.assert(self.isVector());
        std.debug.assert(other.isVector());
        return @reduce(.Add, self.vec * other.vec);
    }

    pub fn cross(self: This, other: This) This {
        std.debug.assert(self.isVector());
        std.debug.assert(other.isVector());

        const first = @Vector(4, f64){ self.y() * other.z(), self.z() * other.x(), self.x() * other.y(), 0 };
        const second = @Vector(4, f64){ self.z() * other.y(), self.x() * other.z(), self.y() * other.x(), 0 };
        return This{ .vec = first - second };
    }

    pub fn reflected(self: This, normal: This) This {
        const two: f64 = 2.0;
        const vtwo = @splat(4, two);
        return This{ .vec = self.vec - normal.vec * vtwo * @splat(4, @reduce(.Add, self.vec * normal.vec)) };
    }

    pub fn x(self: This) f64 {
        return self.vec[0];
    }

    pub fn y(self: This) f64 {
        return self.vec[1];
    }

    pub fn z(self: This) f64 {
        return self.vec[2];
    }

    pub fn isPoint(self: This) bool {
        return self.vec[3] == 1;
    }

    pub fn isVector(self: This) bool {
        return self.vec[3] == 0;
    }

    vec: @Vector(4, f64),

    const This = @This();
};

pub const Point = struct {
    pub fn init(x: f64, y: f64, z: f64) Tuple {
        return Tuple.init(x, y, z, 1);
    }
};

pub const Vector = struct {
    pub fn init(x: f64, y: f64, z: f64) Tuple {
        return Tuple.init(x, y, z, 0);
    }
};

const expect = std.testing.expect;

test "A tuple with w == 1 is a point" {
    const p = Tuple.init(4.3, -4.2, 3.1, 1.0);

    try expect(p.x() == 4.3);
    try expect(p.y() == -4.2);
    try expect(p.z() == 3.1);
    try expect(p.isPoint() == true);
    try expect(p.isVector() == false);
}

test "A tuple with w == 0 is a vector" {
    const v = Tuple.init(5.3, -5.2, 4.1, 0.0);

    try expect(v.x() == 5.3);
    try expect(v.y() == -5.2);
    try expect(v.z() == 4.1);
    try expect(v.isPoint() == false);
    try expect(v.isVector() == true);
}

test "Point.init creates a correct tuple" {
    const p = Point.init(3, 6.0, -11);

    try expect(p.x() == 3);
    try expect(p.y() == 6.0);
    try expect(p.z() == -11);
    try expect(p.isPoint() == true);
    try expect(p.isVector() == false);
}

test "Vector.init creates a correct tuple" {
    const v = Vector.init(9, 11.11, -3.23);

    try expect(v.x() == 9);
    try expect(v.y() == 11.11);
    try expect(v.z() == -3.23);
    try expect(v.isPoint() == false);
    try expect(v.isVector() == true);
}

test "Tuple.equals: point != vector" {
    const v = Vector.init(1, 2, 3);
    const p = Point.init(1, 2, 3);

    try expect(!v.equals(p));
    try expect(!p.equals(v));
}

test "Tuple.equals: point == point" {
    const p = Point.init(1, 2, 3);
    const p2 = Point.init(1, 2, 3);

    try expect(p.equals(p2));
    try expect(p2.equals(p));
}

test "Tuple.equals: vector == vector" {
    const v = Vector.init(1, 2, 3);
    const v2 = Vector.init(1, 2, 3);

    try expect(v.equals(v2));
    try expect(v2.equals(v));
}

test "Tuple.equals: exact" {
    const p = Point.init(1, 2, 3);
    const p2 = Point.init(1, 2, 3);

    try expect(p.equals(p));
    try expect(p.equals(p2));
    try expect(p2.equals(p));
}

test "Tuple.equals: epsilon" {
    const t = Tuple.init(1.1, 4.4, 6.0, -19);
    const e = std.math.floatEps(f64);
    var t2 = t;

    t2 = t2.plus(Tuple.init(e, 0, 0, 0));
    try expect(t.equals(t2));
    try expect(t2.equals(t));

    t2 = t2.plus(Tuple.init(0, e, 0, 0));
    try expect(t.equals(t2));
    try expect(t2.equals(t));

    t2 = t2.plus(Tuple.init(0, 0, e, 0));
    try expect(t.equals(t2));
    try expect(t2.equals(t));

    t2 = t2.plus(Tuple.init(0, 0, 0, e));
    try expect(t.equals(t2));
    try expect(t2.equals(t));
}

test "Tuple.plus: point plus vector equals point" {
    const p = Point.init(1, 2, 3);
    const v = Vector.init(4, 5, 6);
    const csum = Tuple.init(5, 7, 9, 1);

    const sum = p.plus(v);
    const rsum = v.plus(p);

    errdefer (std.debug.print("{} {} {}\n", .{ sum, rsum, csum }));

    try expect(sum.isPoint());
    try expect(rsum.isPoint());

    try expect(sum.equals(rsum));
    try expect(rsum.equals(sum));

    try expect(sum.equals(csum));
    try expect(rsum.equals(csum));
}

test "Tuple.plus: vector plus vector equals vector" {
    const v = Vector.init(1, 2, 3);
    const v2 = Vector.init(4, 5, 6);
    const csum = Tuple.init(5, 7, 9, 0);

    const sum = v.plus(v2);
    const rsum = v2.plus(v);

    errdefer (std.debug.print("{} {} {}\n", .{ sum, rsum, csum }));

    try expect(sum.isVector());
    try expect(rsum.isVector());

    try expect(sum.equals(rsum));
    try expect(rsum.equals(sum));

    try expect(sum.equals(csum));
    try expect(rsum.equals(csum));
}

test "Tuple.minus: point minus vector equals point" {
    const p = Point.init(1, 2, 3);
    const v = Vector.init(4, 5, 6);
    const cdiff = Tuple.init(-3, -3, -3, 1);

    const diff = p.minus(v);

    errdefer (std.debug.print("{} {}\n", .{ diff, cdiff }));

    try expect(diff.isPoint());
    try expect(diff.equals(cdiff));
}

test "Tuple.minus: vector minus vector equals vector" {
    const v = Vector.init(1, 2, 3);
    const v2 = Vector.init(4, 5, 6);
    const cdiff = Tuple.init(-3, -3, -3, 0);
    const rcdiff = Tuple.init(3, 3, 3, 0);

    const diff = v.minus(v2);
    const rdiff = v2.minus(v);

    errdefer (std.debug.print("{} {} {} {}\n", .{ diff, rdiff, cdiff, rcdiff }));

    try expect(diff.isVector());
    try expect(rdiff.isVector());

    try expect(diff.equals(cdiff));
    try expect(rdiff.equals(rcdiff));
}

test "Tuple.minus: subtracting from the zero vector" {
    const zero = Vector.init(0, 0, 0);
    const v = Vector.init(1, -2, 3);
    const neg = zero.minus(v);

    try expect(neg.equals(Vector.init(-1, 2, -3)));
}

test "Tuple.negated" {
    const t = Tuple.init(1, -2, 3, -4);
    const neg = t.negated();

    try expect(neg.equals(Tuple.init(-1, 2, -3, 4)));
}

test "Tuple.scaled: > 1" {
    const t = Tuple.init(1, -2, 3, -4);
    const sca = t.scaled(3.5);

    try expect(sca.equals(Tuple.init(3.5, -7, 10.5, -14)));
}

test "Tuple.scaled: < 1" {
    const t = Tuple.init(1, -2, 3, -4);
    const sca = t.scaled(0.5);

    try expect(sca.equals(Tuple.init(0.5, -1, 1.5, -2)));
}

test "Tuple.divided" {
    const t = Tuple.init(1, -2, 3, -4);
    const div = t.divided(2);

    try expect(div.equals(Tuple.init(0.5, -1, 1.5, -2)));
}

test "Vector.magnitude of (0, 0, 0)" {
    const v = Vector.init(0, 0, 0);
    try expect(v.magnitude() == 0);
}

test "Vector.magnitude of (1, 0, 0)" {
    const v = Vector.init(1, 0, 0);
    try expect(v.magnitude() == 1);
}

test "Vector.magnitude of (0, 1, 0)" {
    const v = Vector.init(0, 1, 0);
    try expect(v.magnitude() == 1);
}

test "Vector.magnitude of (0, 0, 1)" {
    const v = Vector.init(0, 0, 1);
    try expect(v.magnitude() == 1);
}

test "Vector.magnitude of (1, 2, 3)" {
    const v = Vector.init(1, 2, 3);
    try expect(v.magnitude() == @sqrt(14.0));
}

test "Vector.magnitude of (-1, -2, -3)" {
    const v = Vector.init(-1, -2, -3);
    try expect(v.magnitude() == @sqrt(14.0));
}

test "Vector.normalized (4, 0, 0)" {
    const v = Vector.init(4, 0, 0);
    try expect(v.normalized().equals(Vector.init(1, 0, 0)));
}

test "Vector.normalized (4, 0, 0) has magnitude 1" {
    const v = Vector.init(4, 0, 0);
    try expect(v.normalized().magnitude() == 1.0);
}

test "Vector.normalized (1, 2, 3)" {
    const v = Vector.init(1, 2, 3);
    const normalized = Vector.init(1.0 / @sqrt(14.0), 2.0 / @sqrt(14.0), 3.0 / @sqrt(14.0));
    try expect(v.normalized().equals(normalized));
}

test "Vector.dot" {
    const v = Vector.init(1, 2, 3);
    const u = Vector.init(2, 3, 4);
    try expect(v.dot(u) == 20);
}

test "Vector.cross" {
    const v = Vector.init(1, 2, 3);
    const u = Vector.init(2, 3, 4);
    try expect(v.cross(u).equals(Vector.init(-1, 2, -1)));
    try expect(u.cross(v).equals(Vector.init(1, -2, 1)));
}

test "Reflecting a vector at 45d" {
    const v = Vector.init(1, -1, 0);
    const n = Vector.init(0, 1, 0);
    const r = v.reflected(n);
    try expect(r.equals(Vector.init(1, 1, 0)));
}

test "Reflecting a vector off a slanted surfact" {
    const v = Vector.init(0, -1, 0);
    const n = Vector.init(@sqrt(2.0) / 2.0, @sqrt(2.0) / 2.0, 0);
    const r = v.reflected(n);
    try expect(r.equals(Vector.init(1, 0, 0)));
}
