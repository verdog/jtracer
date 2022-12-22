const std = @import("std");

const mtx = @import("matrix.zig");
const trans = @import("transform.zig");

const Tuple = @import("tuple.zig").Tuple;
const Point = @import("tuple.zig").Point;
const Vector = @import("tuple.zig").Vector;

const Material = @import("material.zig").Material;

pub const Sphere = struct {
    var next_id: u16 = 1;

    id: u16,
    transform: mtx.Matrix(4, 4),
    // TODO make this field more data-oriented
    material: Material,

    pub fn init() This {
        const id = This.next_id;
        This.next_id += 1;
        return .{
            .id = id,
            .transform = mtx.Matrix(4, 4).identity(),
            .material = Material.init(),
        };
    }

    pub fn normalAt(self: This, point: Tuple) Tuple {
        const trans_invs = self.transform.inverted() catch unreachable;
        const obj_space_point = trans_invs.mult(point);
        const obj_space_normal = obj_space_point.minus(Point.init(0, 0, 0));
        var wld_space_normal = trans_invs.transposed().mult(obj_space_normal);
        // hack: technically we should multiply by submatrix(self.transform, 3, 3), to avoid
        //       messing up the w component of the resulting vector. we can get away with
        //       the above as long as we set the w component back to 0 manually
        wld_space_normal.vec[3] = 0;
        return wld_space_normal.normalized();
    }

    const This = @This();
};

const expect = std.testing.expect;
const print = @import("u.zig").print;

test "Spheres have unique ids" {
    const s1 = Sphere.init();
    const s2 = Sphere.init();

    try expect(s1.id != s2.id);
}

test "Sphere's have a default transform" {
    const s = Sphere.init();

    try expect(s.transform.equals(mtx.Matrix(4, 4).identity()));
}

test "Changing a sphere's transformation" {
    var s = Sphere.init();
    var t = trans.makeTranslation(2, 3, 4);
    s.transform = t;

    try expect(s.transform.equals(t));
}

test "The normal on a sphere at a point on the x axis" {
    const s = Sphere.init();
    const n = s.normalAt(Point.init(1, 0, 0));

    try expect(n.equals(Vector.init(1, 0, 0)));
}

test "The normal on a sphere at a point on the y axis" {
    const s = Sphere.init();
    const n = s.normalAt(Point.init(0, 1, 0));

    try expect(n.equals(Vector.init(0, 1, 0)));
}

test "The normal on a sphere at a point on the z axis" {
    const s = Sphere.init();
    const n = s.normalAt(Point.init(0, 0, 1));

    try expect(n.equals(Vector.init(0, 0, 1)));
}
test "The normal on a sphere at a nonaxial point" {
    const s = Sphere.init();
    const n = s.normalAt(Point.init(@sqrt(3.0) / 3.0, @sqrt(3.0) / 3.0, @sqrt(3.0) / 3.0));

    try expect(n.equals(Vector.init(@sqrt(3.0) / 3.0, @sqrt(3.0) / 3.0, @sqrt(3.0) / 3.0)));
}

test "The normal on a sphere is a normalized vector" {
    const s = Sphere.init();
    const n = s.normalAt(Point.init(@sqrt(3.0) / 3.0, @sqrt(3.0) / 3.0, @sqrt(3.0) / 3.0));

    try expect(n.equals(n.normalized()));
}

test "The normal on a translated sphere" {
    var s = Sphere.init();
    s.transform = trans.makeTranslation(0, 1, 0);

    const n = s.normalAt(Point.init(0, 1.70711, -0.70711));

    errdefer print(n);

    // book examples are much less precise than f64
    try expect(n.equalsTolerance(Vector.init(0, 0.70711, -0.70711), 100_000_000_000));
}

test "The normal on a transformed sphere" {
    var s = Sphere.init();
    s.transform = trans.makeScaling(1, 0.5, 1).mult(trans.makeRotationZ(std.math.pi / 5.0));

    const n = s.normalAt(Point.init(0, @sqrt(2.0) / 2.0, -@sqrt(2.0) / 2.0));

    errdefer print(n);

    // book examples are much less precise than f64
    try expect(n.equalsTolerance(Vector.init(0, 0.97014, -0.24254), 100_000_000_000));
}

test "A sphere has a default material" {
    const s = Sphere.init();
    const m = s.material;

    try expect(std.meta.eql(s.material, m));
}

test "A sphere may be assigned a material" {
    var s = Sphere.init();
    var m = Material.init();
    m.ambient = 1;
    s.material = m;

    try expect(std.meta.eql(s.material, m));
}
