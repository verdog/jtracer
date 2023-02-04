//! definitions of volumes that can be placed in a scene

pub const Sphere = struct {
    transform: trans.Transform,

    // TODO make this field more data-oriented
    material: Material,

    pub fn init() This {
        return .{
            .transform = trans.Transform{},
            .material = Material.init(),
        };
    }

    pub fn normalAt(self: This, point: Tuple) Tuple {
        const trans_invs = self.transform.inverse;
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

pub const Plane = struct {
    transform: trans.Transform,
    material: Material,

    pub fn init() This {
        return .{
            .transform = trans.Transform{},
            .material = Material.init(),
        };
    }

    pub fn normalAt(self: This, point: Tuple) Tuple {
        _ = point; // the normal vector on a plane is the same everywhere
        return self.transform.t.mult(Vector.init(0, 1, 0)).normalized();
    }

    const This = @This();
};

pub const Cube = struct {
    transform: trans.Transform,
    material: Material,

    pub fn init() This {
        return .{
            .transform = trans.Transform{},
            .material = Material.init(),
        };
    }

    pub fn normalAt(self: This, world_space_point: Tuple) Tuple {
        const obj_space_point = self.transform.inverse.mult(world_space_point);

        const obj_space_normal = blk: {
            // cubes' normals are always the same for a given face. so finding the normal on a
            // cube reduces to finding the face that the point is on. there are two patterns
            // we can exploit:
            //   1. a point on the cube will always have a component with either 1.0 or -1.0
            //   2. the 1.0 (or -1.0) will always be the greatest (or smallest) of all the components.
            // e.g.
            //   (-1.0, 0.3, -0.5) => -x
            //   (0.3, -0.9, 1.0) => +z
            //   (-0.6, 1.0, 0.7) => +y
            //   (0.4, -1.0, 0.2) => -y
            // because of floating point rounding, we can't assume the component will exactly
            // equal one, so we use the second pattern in practice. for corners, we return the
            // normal from the +-x face.

            const max_component = std.math.max3(
                @fabs(obj_space_point.x()),
                @fabs(obj_space_point.y()),
                @fabs(obj_space_point.z()),
            );

            if (max_component == @fabs(obj_space_point.x()))
                break :blk Vector.init(obj_space_point.x(), 0, 0);
            if (max_component == @fabs(obj_space_point.y()))
                break :blk Vector.init(0, obj_space_point.y(), 0);
            // then assume max_component == @fabs(obj_space_point.z())
            break :blk Vector.init(0, 0, obj_space_point.z());
        };

        var world_space_normal = self.transform.inverse.transposed().mult(obj_space_normal);
        // hack: technically we should multiply by submatrix(self.transform, 3, 3), to avoid
        //       messing up the w component of the resulting vector. we can get away with
        //       the above as long as we set the w component back to 0 manually
        world_space_normal.vec[3] = 0;

        return world_space_normal.normalized();
    }

    const This = @This();
};

pub const Cylinder = struct {
    transform: trans.Transform,
    material: Material,
    length: f64,
    closed: bool,

    pub fn init() This {
        return .{
            .transform = trans.Transform{},
            .material = Material.init(),
            .length = std.math.inf(f64),
            .closed = false,
        };
    }

    pub fn normalAt(self: This, world_space_point: Tuple) Tuple {
        // into object space
        const obj_space_point = self.transform.inverse.mult(world_space_point);

        const obj_space_normal = blk: {
            // first check if we are on a cap or a wall
            const dist = std.math.pow(f64, obj_space_point.x(), 2) + std.math.pow(f64, obj_space_point.z(), 2);
            if (dist < 1 and obj_space_point.y() >= self.length / 2.0 - mymath.floatTolerance) {
                // on top cap
                break :blk Vector.init(0, 1, 0);
            }
            if (dist < 1 and obj_space_point.y() <= -self.length / 2.0 + mymath.floatTolerance) {
                // on bottom cap
                break :blk Vector.init(0, -1, 0);
            }

            // not on cap...

            // in normal calculation, even if not technically true, we treat cylinders as
            // having infinite height, and assume that any point passed into this function is
            // on the surface of the cylinder. to get the normal under those conditions, simply
            // remove the y component of the point.
            break :blk Vector.init(obj_space_point.x(), 0, obj_space_point.z());
        };

        // back to world space
        var world_space_normal = self.transform.inverse.transposed().mult(obj_space_normal);
        // hack: technically we should multiply by submatrix(self.transform, 3, 3), to avoid
        //       messing up the w component of the resulting vector. we can get away with
        //       the above as long as we set the w component back to 0 manually
        world_space_normal.vec[3] = 0;

        return world_space_normal.normalized();
    }

    const This = @This();
};

pub const Cone = struct {
    transform: trans.Transform,
    material: Material,
    min: f64,
    max: f64,
    closed: bool,

    pub fn init() This {
        return .{
            .transform = trans.Transform{},
            .material = Material.init(),
            .min = -std.math.inf(f64),
            .max = std.math.inf(f64),
            .closed = false,
        };
    }

    pub fn normalAt(self: This, world_space_point: Tuple) Tuple {
        // into object space
        const obj_space_point = self.transform.inverse.mult(world_space_point);

        const obj_space_normal = blk: {
            if (std.math.approxEqRel(f64, obj_space_point.x(), 0, mymath.floatTolerance) and
                std.math.approxEqRel(f64, obj_space_point.y(), 0, mymath.floatTolerance) and
                std.math.approxEqRel(f64, obj_space_point.z(), 0, mymath.floatTolerance))
            {
                // special case: (0, 0, 0). arbitrarily pick (1, 0, 0)
                break :blk Vector.init(1, 0, 0);
            }

            // first check if we are on a cap or a wall
            const dist = std.math.pow(f64, obj_space_point.x(), 2) + std.math.pow(f64, obj_space_point.z(), 2);
            if (dist < std.math.pow(f64, obj_space_point.y(), 2) and obj_space_point.y() >= self.max - mymath.floatTolerance) {
                // on top cap
                break :blk Vector.init(0, 1, 0);
            }
            if (dist < std.math.pow(f64, obj_space_point.y(), 2) and obj_space_point.y() <= self.min + mymath.floatTolerance) {
                // on bottom cap
                break :blk Vector.init(0, -1, 0);
            }

            // not on cap...

            const y = yblk: {
                var ry = @sqrt(std.math.pow(f64, obj_space_point.x(), 2) + std.math.pow(f64, obj_space_point.z(), 2));
                break :yblk if (obj_space_point.y() > 0) -ry else ry;
            };
            break :blk Vector.init(obj_space_point.x(), y, obj_space_point.z());
        };

        // back to world space
        var world_space_normal = self.transform.inverse.transposed().mult(obj_space_normal);
        // hack: technically we should multiply by submatrix(self.transform, 3, 3), to avoid
        //       messing up the w component of the resulting vector. we can get away with
        //       the above as long as we set the w component back to 0 manually
        world_space_normal.vec[3] = 0;

        return world_space_normal.normalized();
    }

    const This = @This();
};

pub const Triangle = struct {
    p1: Tuple,
    p2: Tuple,
    p3: Tuple,
    e1: Tuple,
    e2: Tuple,
    normal: Tuple,

    transform: trans.Transform,
    material: Material,

    pub fn init(p1: Tuple, p2: Tuple, p3: Tuple) This {
        std.debug.assert(p1.isPoint());
        std.debug.assert(p2.isPoint());
        std.debug.assert(p3.isPoint());

        const e1 = p2.minus(p1);
        const e2 = p3.minus(p1);

        std.debug.assert(e1.isVector());
        std.debug.assert(e2.isVector());

        return .{
            .p1 = p1,
            .p2 = p2,
            .p3 = p3,
            .e1 = e1,
            .e2 = e2,
            .normal = e2.cross(e1).normalized(),
            .transform = trans.Transform{},
            .material = Material.init(),
        };
    }

    pub fn normalAt(self: This, world_space_point: Tuple) Tuple {
        // normal is the same everywhere
        _ = world_space_point;
        return self.transform.t.mult(self.normal).normalized();
    }

    pub fn equals(self: This, other: This) bool {
        return self.p1.equals(other.p1) and
            self.p2.equals(other.p2) and
            self.p3.equals(other.p3);
    }

    const This = @This();
};

test "The normal of a plane is constant everywhere" {
    const p = Plane.init();

    const n1 = p.normalAt(Point.init(0, 0, 0));
    const n2 = p.normalAt(Point.init(10, 0, -10));
    const n3 = p.normalAt(Point.init(100000, 0, 3));

    const n = Vector.init(0, 1, 0);

    errdefer print(n);
    errdefer print(n3);
    errdefer print(n2);
    errdefer print(n1);

    try expect(n.equals(n1));
    try expect(n.equals(n2));
    try expect(n.equals(n3));
}

test "The normal of a scaled plane is constant everywhere" {
    var p = Plane.init();
    p.transform = p.transform.mult(trans.makeScaling(1, 10, 1));

    const n1 = p.normalAt(Point.init(0, 0, 0));
    const n2 = p.normalAt(Point.init(10, 0, -10));
    const n3 = p.normalAt(Point.init(100000, 0, 3));

    const n = Vector.init(0, 1, 0);

    errdefer print(n3);
    errdefer print(n2);
    errdefer print(n1);
    errdefer print(n);

    try expect(n.equals(n1));
    try expect(n.equals(n2));
    try expect(n.equals(n3));
}

test "The normal of a rotated plane is constant everywhere" {
    var p = Plane.init();
    p.transform = p.transform.mult(trans.makeRotationX(std.math.pi / 2.0));

    const n1 = p.normalAt(Point.init(0, 0, 0));
    const n2 = p.normalAt(Point.init(10, 0, -10));
    const n3 = p.normalAt(Point.init(100000, 0, 3));

    const n = Vector.init(0, 0, 1);

    errdefer print(n3);
    errdefer print(n2);
    errdefer print(n1);
    errdefer print(n);

    try expect(n.equals(n1));
    try expect(n.equals(n2));
    try expect(n.equals(n3));
}

test "The normal on a cube" {
    const tst = struct {
        fn t(p: Tuple, expected_n: Tuple) !void {
            const c = Cube.init();
            const n = c.normalAt(p);
            errdefer {
                print(n);
                print(expected_n);
            }
            try expect(n.equals(expected_n));
        }
    }.t;

    try tst(Point.init(1, 0.5, -0.8), Vector.init(1, 0, 0));
    try tst(Point.init(-1, -0.2, 0.9), Vector.init(-1, 0, 0));
    try tst(Point.init(-0.4, 1, -0.1), Vector.init(0, 1, 0));
    try tst(Point.init(0.3, -1, -0.7), Vector.init(0, -1, 0));
    try tst(Point.init(-0.6, 0.3, 1), Vector.init(0, 0, 1));
    try tst(Point.init(0.4, 0.4, -1), Vector.init(0, 0, -1));
    try tst(Point.init(1, 1, 1), Vector.init(1, 0, 0)); // corners assume x face
    try tst(Point.init(-1, -1, -1), Vector.init(-1, 0, 0)); // corners assume x face
}

test "The normal on a scaled cube" {
    const tst = struct {
        fn t(p: Tuple, expected_n: Tuple) !void {
            var c = Cube.init();
            c.transform = c.transform.mult(trans.makeScaling(0.5, 1, 0.5));
            const n = c.normalAt(p);
            errdefer {
                print(n);
                print(expected_n);
            }
            try expect(n.equals(expected_n));
        }
    }.t;

    try tst(Point.init(0.5, 0.3, -0.2), Vector.init(1, 0, 0));
    try tst(Point.init(0.24, -1, 0.09), Vector.init(0, -1, 0));
    try tst(Point.init(-0.4, 0.3, 0.5), Vector.init(0, 0, 1));
}

test "Sphere's have a default transform" {
    const s = Sphere.init();

    try expect(s.transform.t.equals(mtx.Matrix(4, 4).identity()));
}

test "Changing a sphere's transformation" {
    var s = Sphere.init();
    var t = trans.makeTranslation(2, 3, 4);
    s.transform = t;

    try expect(s.transform.t.equals(t.t));
    try expect(s.transform.inverse.equals(t.inverse));
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

    const scale = trans.makeScaling(1, 0.5, 1);
    const rot = trans.makeRotationZ(std.math.pi / 5.0);

    s.transform = s.transform.chain(.{ scale, rot });

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

test "The normal on a cylinder" {
    const tst = struct {
        fn t(p: Tuple, expected_n: Tuple) !void {
            const c = Cylinder.init();
            const n = c.normalAt(p);
            errdefer {
                print(n);
                print(expected_n);
            }
            try expect(n.equals(expected_n));
        }
    }.t;

    try tst(Point.init(1, 0, 0), Vector.init(1, 0, 0));
    try tst(Point.init(0, 5, -1), Vector.init(0, 0, -1));
    try tst(Point.init(0, -2, 1), Vector.init(0, 0, 1));
    try tst(Point.init(-1, 1, 0), Vector.init(-1, 0, 0));
    try tst(
        Point.init(@sqrt(2.0) / 2.0, 10, @sqrt(2.0) / 2.0),
        Vector.init(@sqrt(2.0) / 2.0, 0, @sqrt(2.0) / 2.0),
    );
}

test "The normal on a scaled cylinder" {
    const tst = struct {
        fn t(p: Tuple, expected_n: Tuple) !void {
            var c = Cylinder.init();
            c.transform = c.transform.mult(trans.makeScaling(2, 1, 2));
            const n = c.normalAt(p);
            errdefer {
                print(n);
                print(expected_n);
            }
            try expect(n.equals(expected_n));
        }
    }.t;

    try tst(Point.init(2, 0, 0), Vector.init(1, 0, 0));
    try tst(Point.init(0, 5, -2), Vector.init(0, 0, -1));
    try tst(Point.init(0, -2, 2), Vector.init(0, 0, 1));
    try tst(Point.init(-2, 1, 0), Vector.init(-1, 0, 0));
}

test "The normal on a rotated cylinder" {
    const tst = struct {
        fn t(p: Tuple, expected_n: Tuple) !void {
            var c = Cylinder.init();
            c.transform = c.transform.mult(trans.makeRotationZ(std.math.pi / 2.0));
            const n = c.normalAt(p);
            errdefer {
                print(n);
                print(expected_n);
            }
            try expect(n.equals(expected_n));
        }
    }.t;

    try tst(Point.init(0, 1, 0), Vector.init(0, 1, 0));
    try tst(Point.init(5, 1, 0), Vector.init(0, 1, 0));
    try tst(Point.init(5, 0, -1), Vector.init(0, 0, -1));
    try tst(Point.init(5, 0, 1), Vector.init(0, 0, 1));
}

test "The normal on the caps of a cylinder" {
    const tst = struct {
        fn t(p: Tuple, expected_n: Tuple) !void {
            var c = Cylinder.init();
            c.length = 1;
            c.closed = true;
            // c.transform = c.transform.mult(trans.makeRotationZ(std.math.pi / 2.0));
            const n = c.normalAt(p);
            errdefer {
                print(n);
                print(expected_n);
            }
            try expect(n.equals(expected_n));
        }
    }.t;

    try tst(Point.init(0, 1, 0), Vector.init(0, 1, 0));
    try tst(Point.init(0.5, 1, 0), Vector.init(0, 1, 0));
    try tst(Point.init(0, 1, 0.5), Vector.init(0, 1, 0));
    try tst(Point.init(0, -1, 0), Vector.init(0, -1, 0));
    try tst(Point.init(-0.5, -1, 0), Vector.init(0, -1, 0));
    try tst(Point.init(0, -1, -0.5), Vector.init(0, -1, 0));
}

test "A cylinder has a length" {
    const cyl = Cylinder.init();
    try expect(cyl.length == std.math.inf(f64));
}

test "A cylinder can be closed or open" {
    const cyl = Cylinder.init();
    try expect(cyl.closed == false);
}

test "The normal of a cone" {
    const tst = struct {
        fn t(p: Tuple, expected_n: Tuple) !void {
            var c = Cone.init();
            c.closed = true;
            c.min = -2;
            c.max = 2;
            const n = c.normalAt(p);
            errdefer {
                print(n);
                print(expected_n);
            }
            try expect(n.equals(expected_n));
        }
    }.t;

    // normal at origin: always (1, 0, 0)
    try tst(Point.init(0, 0, 0), Vector.init(1, 0, 0));
    // on sides
    try tst(Point.init(1, 1, 1), Vector.init(1, -@sqrt(2.0), 1).normalized());
    try tst(Point.init(-1, -1, 0), Vector.init(-1, 1, 0).normalized());
    // on caps
    try tst(Point.init(0.5, 2, 0.5), Vector.init(0, 1, 0));
    try tst(Point.init(-0.5, -2, -0.5), Vector.init(0, -1, 0));
    try tst(Point.init(0, -2, 1.99), Vector.init(0, -1, 0));
}

test "Constructing a triangle" {
    const p1 = Point.init(0, 1, 0);
    const p2 = Point.init(-1, 0, 0);
    const p3 = Point.init(1, 0, 0);

    const t = Triangle.init(p1, p2, p3);

    try expect(t.p1.equals(p1));
    try expect(t.p2.equals(p2));
    try expect(t.p3.equals(p3));
    try expect(t.e1.equals(Vector.init(-1, -1, 0)));
    try expect(t.e2.equals(Vector.init(1, -1, 0)));
    try expect(t.normal.equals(Vector.init(0, 0, -1)));
}

test "The normal vector on a triangle is the same everywhere" {
    const p1 = Point.init(0, 1, 0);
    const p2 = Point.init(-1, 0, 0);
    const p3 = Point.init(1, 0, 0);

    const t = Triangle.init(p1, p2, p3);

    try expect(t.normalAt(Point.init(0, 0.5, 0)).equals(t.normal));
    try expect(t.normalAt(Point.init(-0.5, 0.75, 0)).equals(t.normal));
    try expect(t.normalAt(Point.init(0.5, 0.25, 0)).equals(t.normal));
}

const std = @import("std");

const mtx = @import("matrix.zig");
const trans = @import("transform.zig");
const mymath = @import("mymath.zig");

const Tuple = @import("tuple.zig").Tuple;
const Point = @import("tuple.zig").Point;
const Vector = @import("tuple.zig").Vector;

const Material = @import("material.zig").Material;

const expect = std.testing.expect;
const print = @import("u.zig").print;
