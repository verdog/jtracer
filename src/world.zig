//! holds all volumes

const std = @import("std");

pub const World = struct {
    pub fn init(alctr: std.mem.Allocator) This {
        return .{
            .spheres_buf = std.ArrayList(Sphere).init(alctr),
            .lights_buf = std.ArrayList(PointLight).init(alctr),
        };
    }

    pub fn deinit(self: This) void {
        self.spheres_buf.deinit();
        self.lights_buf.deinit();
    }
    pub fn addVolume(self: *This, comptime T: type) struct { handle: VolumePtr, ptr: *T } {
        switch (T) {
            Sphere => {
                self.spheres_buf.append(Sphere.init()) catch unreachable;
                const last = self.spheres_buf.items.len - 1;
                return .{ .handle = VolumePtr{ .sphere_idx = last }, .ptr = &self.spheres_buf.items[last] };
            },
            else => unreachable,
        }
    }

    pub fn addLight(self: *This, comptime T: type) struct { handle: LightPtr, ptr: *T } {
        switch (T) {
            PointLight => {
                const point = @import("tuple.zig").Point.init(0, 0, 0);
                const color = @import("color.zig").Color.init(1, 1, 1);
                self.lights_buf.append(PointLight.init(point, color)) catch unreachable;
                const last = self.lights_buf.items.len - 1;
                return .{
                    .handle = LightPtr{ .light_idx = last },
                    .ptr = &self.lights_buf.items[last],
                };
            },
            else => unreachable,
        }
    }

    pub fn getVolume(self: This, comptime T: type, vptr: VolumePtr) *T {
        return switch (T) {
            Sphere => return &self.spheres_buf.items[vptr.sphere_idx],
            else => unreachable,
        };
    }

    pub fn getLight(self: This, comptime T: type, lptr: LightPtr) *T {
        return switch (T) {
            PointLight => return &self.lights_buf.items[lptr.light_idx],
            else => unreachable,
        };
    }

    pub fn isShadowed(self: This, point: Tuple, lht: PointLight, alctr: std.mem.Allocator) bool {
        const p2l = lht.position.minus(point);
        const distance = p2l.magnitude();

        var ixs = Intersections.init(alctr, .{});
        defer ixs.deinit();
        self.intersect(Ray.init(point, p2l.normalized()), &ixs);

        if (ixs.hit()) |hit| {
            return hit.t < distance;
        }

        return false;
    }

    pub fn intersect(self: This, ray: Ray, ixs: *Intersections) void {
        ixs.clear();

        for (self.spheres_buf.items) |*sphptr, i| {
            const vptr = VolumePtr{ .sphere_idx = i };
            ixs.intersect(sphptr.*, vptr, ray);
        }

        ixs.order();
    }

    pub fn shadeHit(self: This, data: HitData, alctr: std.mem.Allocator) Color {
        // TODO smart shape lookup
        const s = self.getVolume(Sphere, data.vptr);

        var color = Color.init(0, 0, 0);

        for (self.lights_buf.items) |l| {
            const in_shadow = self.isShadowed(data.over_point, l, alctr);
            color = color.plus(light.lighting(s.material, l, data.point, data.eye_vector, data.normal_vector, in_shadow));
        }

        return color;
    }

    pub fn colorAt(self: This, ray: Ray, alctr: std.mem.Allocator) Color {
        var ixs = Intersections.init(alctr, .{});
        defer ixs.deinit();

        self.intersect(ray, &ixs);

        if (ixs.hit()) |hit| {
            const s = self.getVolume(Sphere, hit.vptr);
            const data = HitData.init(ray, hit, s);
            return self.shadeHit(data, alctr);
        } else {
            return Color.init(0, 0, 0);
        }
    }

    pub fn render(self: This, cam: Camera, qan: *Qanvas, alctr: std.mem.Allocator) void {
        rdr.startRenderEngine(self, cam, qan, alctr);
    }

    pub fn numVolumes(self: This) usize {
        return self.spheres_buf.items.len;
    }

    pub fn numLights(self: This) usize {
        return self.lights_buf.items.len;
    }

    // TODO change these to a list that won't move data
    // around upon extension (e.g. a chunky linked list)?
    spheres_buf: std.ArrayList(Sphere),
    lights_buf: std.ArrayList(PointLight),

    const This = @This();
};

pub const Camera = struct {
    width: i64,
    height: i64,
    half_width: f64,
    half_height: f64,
    fov: f64,
    transform: trans.Transform,
    pixel_size: f64,

    pub fn init(width: i64, height: i64, fov: f64) This {
        const aspect = @intToFloat(f64, width) / @intToFloat(f64, height);

        var half_width: f64 = undefined;
        var half_height: f64 = undefined;

        const pixel_size = blk: {
            const half_view = @tan(fov / 2.0);

            if (aspect >= 1.0) {
                half_width = half_view;
                half_height = half_view / aspect;
            } else {
                half_width = half_view * aspect;
                half_height = half_view;
            }

            break :blk (half_width * 2.0) / @intToFloat(f64, width);
        };

        return .{
            .width = width,
            .height = height,
            .half_width = half_width,
            .half_height = half_height,
            .fov = fov,
            .transform = trans.Transform{},
            .pixel_size = pixel_size,
        };
    }

    pub fn rayForPixel(self: This, x: i64, y: i64) Ray {
        // page 104
        // the offset from the edge of the canvas to the pixel's center
        const x_offset = (@intToFloat(f64, x) + 0.5) * self.pixel_size;
        const y_offset = (@intToFloat(f64, y) + 0.5) * self.pixel_size;

        // the untransformed coordinates of the pixel in world space.
        // (camera looks toward -z, so +x is to the *left*)
        const world_x = self.half_width - x_offset;
        const world_y = self.half_height - y_offset;

        // using the camera transform, transform the canvas point and the origin,
        // and then compute the ray's direction. (remember that the canvas is
        // always at z=-1)
        const inv = self.transform.inverse;
        const pixel = inv.mult(Point.init(world_x, world_y, -1));
        const origin = inv.mult(Point.init(0, 0, 0));
        const direction = pixel.minus(origin).normalized();

        return Ray.init(origin, direction);
    }

    const This = @This();
};

pub const VolumePtr = union(enum) {
    sphere_idx: usize,
    plane_idx: usize,
};

pub const LightPtr = union(enum) {
    light_idx: usize,
};

const Tuple = @import("tuple.zig").Tuple;
const Point = @import("tuple.zig").Point;
const Vector = @import("tuple.zig").Vector;
const Color = @import("color.zig").Color;

const Sphere = @import("volume.zig").Sphere;
const PointLight = @import("light.zig").PointLight;

const Ray = @import("ray.zig").Ray;
const isct = @import("intersect.zig");
const light = @import("light.zig");
const Intersection = isct.Intersection;
const Intersections = isct.Intersections;
const HitData = isct.HitData;
const Matrix = @import("matrix.zig").Matrix;
const Qanvas = @import("qanvas.zig").Qanvas;

const trans = @import("transform.zig");
const rdr = @import("render.zig");

const expect = std.testing.expect;
const print = @import("u.zig").print;

fn getTestWorld(alctr: std.mem.Allocator) World {
    var w = World.init(alctr);

    { // default light
        var l = w.addLight(PointLight);
        l.ptr.position = Point.init(-10, 10, -10);
    }
    { // first sphere
        var s = w.addVolume(Sphere);
        s.ptr.material.color = Color.init(0.8, 1.0, 0.6);
        s.ptr.material.diffuse = 0.7;
        s.ptr.material.specular = 0.2;
    }
    { // second sphere
        var s = w.addVolume(Sphere);
        s.ptr.transform = trans.makeScaling(0.5, 0.5, 0.5);
    }

    return w;
}

test "World creation/destruction" {
    const alctr = std.testing.allocator;
    var world = World.init(alctr);

    try expect(world.numVolumes() == 0);
    try expect(world.numLights() == 0);

    defer world.deinit();
}

test "World: add/get sphere" {
    const alctr = std.testing.allocator;
    var world = World.init(alctr);
    defer world.deinit();

    const new = world.addVolume(Sphere);
    var sphereptr = world.getVolume(Sphere, new.handle);

    try expect(sphereptr.id != 0);
}

test "The test world" {
    // Some of the unit tests here assume that the test world meets these specs
    const w = getTestWorld(std.testing.allocator);
    defer w.deinit();

    // light
    try expect(w.numLights() == 1);

    try expect(w.lights_buf.items[0].position.equals(Point.init(-10, 10, -10)));
    try expect(w.lights_buf.items[0].intensity.equals(Color.init(1, 1, 1)));

    // spheres
    try expect(w.numVolumes() == 2);

    // sphere 1
    try expect(w.spheres_buf.items[0].material.color.equals(Color.init(0.8, 1.0, 0.6)));
    try expect(w.spheres_buf.items[0].material.diffuse == 0.7);
    try expect(w.spheres_buf.items[0].material.specular == 0.2);

    // sphere 2
    try expect(w.spheres_buf.items[1].transform.t.equals(trans.makeScaling(0.5, 0.5, 0.5).t));
}

test "Intersect a world with a ray" {
    const alctr = std.testing.allocator;

    const w = getTestWorld(alctr);
    defer w.deinit();

    const r = Ray.init(Point.init(0, 0, -5), Vector.init(0, 0, 1));

    var xs = Intersections.init(alctr, .{});
    defer xs.deinit();

    w.intersect(r, &xs);

    errdefer std.debug.print("{any}\n", .{xs.ixs.items});

    try expect(xs.ixs.items.len == 4);
    try expect(xs.ixs.items[0].t == 4.0);
    try expect(xs.ixs.items[1].t == 4.5);
    try expect(xs.ixs.items[2].t == 5.5);
    try expect(xs.ixs.items[3].t == 6.0);
}

test "Shading an intersection" {
    const alctr = std.testing.allocator;

    const w = getTestWorld(alctr);
    defer w.deinit();

    const r = Ray.init(Point.init(0, 0, -5), Vector.init(0, 0, 1));
    const s = w.spheres_buf.items[0];
    const x = Intersection.init(4.0, VolumePtr{ .sphere_idx = 0 });
    const d = HitData.init(r, x, s);

    const c = w.shadeHit(d, alctr);

    defer print(c);

    // TODO book tests are imprecise
    try expect(c.equalsTolerance(Color.init(0.38066, 0.47583, 0.2855), 100_000_000_000));
}

test "Shading an intersection from the inside" {
    const alctr = std.testing.allocator;

    const w = getTestWorld(alctr);
    defer w.deinit();

    w.lights_buf.items[0].position = Point.init(0, 0.25, 0);

    const r = Ray.init(Point.init(0, 0, 0), Vector.init(0, 0, 1));
    const s = w.spheres_buf.items[1];
    const x = Intersection.init(0.5, VolumePtr{ .sphere_idx = 1 });
    const d = HitData.init(r, x, s);

    const c = w.shadeHit(d, alctr);

    // TODO book tests are imprecise
    try expect(c.equalsTolerance(Color.init(0.90498, 0.90498, 0.90498), 100_000_000_000));
}

test "Shading an intersection in shadow" {
    const alctr = std.testing.allocator;
    var w = World.init(alctr);
    defer w.deinit();

    { // default light
        var l = w.addLight(PointLight);
        l.ptr.position = Point.init(0, 0, -10);
    }
    { // first sphere
        _ = w.addVolume(Sphere);
    }
    { // second sphere
        var s = w.addVolume(Sphere);
        s.ptr.transform = trans.makeTranslation(0, 0, 10);
    }

    const r = Ray.init(Point.init(0, 0, 5), Vector.init(0, 0, 1));
    const s = w.spheres_buf.items[1];
    const x = Intersection.init(4, VolumePtr{ .sphere_idx = 1 });
    const d = HitData.init(r, x, s);

    const c = w.shadeHit(d, alctr);

    // TODO book tests are imprecise
    try expect(c.equals(Color.init(0.1, 0.1, 0.1)));
}

test "The color when a ray misses" {
    const alctr = std.testing.allocator;

    const w = getTestWorld(alctr);
    defer w.deinit();

    const r = Ray.init(Point.init(0, 0, -5), Vector.init(0, 1, 0));

    const c = w.colorAt(r, alctr);

    try expect(c.equals(Color.init(0, 0, 0)));
}

test "The color when a ray hits" {
    const alctr = std.testing.allocator;

    const w = getTestWorld(alctr);
    defer w.deinit();

    const r = Ray.init(Point.init(0, 0, -5), Vector.init(0, 0, 1));

    const c = w.colorAt(r, alctr);

    // TODO book tests are imprecise
    try expect(c.equalsTolerance(Color.init(0.38066, 0.47583, 0.2855), 100_000_000_000));
}

test "The color with an intersections behind the ray" {
    const alctr = std.testing.allocator;

    const w = getTestWorld(alctr);
    defer w.deinit();

    w.spheres_buf.items[0].material.ambient = 1;
    w.spheres_buf.items[1].material.ambient = 1;

    const r = Ray.init(Point.init(0, 0, 0.75), Vector.init(0, 0, -1));

    const c = w.colorAt(r, alctr);

    try expect(c.equals(w.spheres_buf.items[1].material.color));
}

test "Constructing a camera" {
    const width = 160;
    const height = 120;
    const fov = std.math.pi / 2.0;

    const cam = Camera.init(width, height, fov);

    try expect(cam.width == 160);
    try expect(cam.height == 120);
    try expect(cam.fov == std.math.pi / 2.0);
    try expect(cam.transform.t.equals(Matrix(4, 4).identity()));
}

test "The pixel size for a horizontal canvas" {
    const cam = Camera.init(200, 125, std.math.pi / 2.0);

    try expect(std.math.approxEqRel(f64, cam.pixel_size, 0.01, std.math.floatEps(f64)));
}

test "The pixel size for a vertical canvas" {
    const cam = Camera.init(125, 200, std.math.pi / 2.0);

    defer print(cam.pixel_size);
    try expect(std.math.approxEqRel(f64, cam.pixel_size, 0.01, std.math.floatEps(f64)));
}

test "A ray through the center of the canvas" {
    const cam = Camera.init(201, 101, std.math.pi / 2.0);
    const ray = cam.rayForPixel(100, 50);

    try expect(ray.origin.equals(Point.init(0, 0, 0)));
    try expect(ray.direction.equals(Vector.init(0, 0, -1)));
}

test "A ray through the corner of the canvas" {
    const cam = Camera.init(201, 101, std.math.pi / 2.0);
    const ray = cam.rayForPixel(0, 0);

    try expect(ray.origin.equals(Point.init(0, 0, 0)));
    // TODO book tests are imprecise
    try expect(ray.direction.equalsTolerance(Vector.init(0.66519, 0.33259, -0.66851), 100_000_000_000));
}

test "A ray when the camera is transformed" {
    var cam = Camera.init(201, 101, std.math.pi / 2.0);
    cam.transform = trans.makeRotationY(std.math.pi / 4.0).mult(trans.makeTranslation(0, -2, 5));
    const ray = cam.rayForPixel(100, 50);

    try expect(ray.origin.equals(Point.init(0, 2, -5)));
    try expect(ray.direction.equals(Vector.init(@sqrt(2.0) / 2.0, 0, -@sqrt(2.0) / 2.0)));
}

test "Rendering a world with a camera" {
    const alctr = std.testing.allocator;
    const w = getTestWorld(alctr);
    defer w.deinit();

    var c = Camera.init(11, 11, std.math.pi / 2.0);
    const from = Point.init(0, 0, -5);
    const to = Point.init(0, 0, 0);
    const up = Vector.init(0, 1, 0);
    c.transform = trans.makeView(from, to, up);

    var qan = try Qanvas.init(alctr, @intCast(usize, c.width), @intCast(usize, c.height));
    defer qan.deinit();

    w.render(c, &qan, alctr);

    // TODO book tests are imprecise
    try expect(qan.at(5, 5).equalsTolerance(Color.init(0.38066, 0.47583, 0.2855), 100_000_000_000));
}

test "There is no shadow when nothing is collinear with point and light" {
    const alctr = std.testing.allocator;
    const w = getTestWorld(alctr);
    defer w.deinit();

    const p = Point.init(0, 10, 0);

    try expect(w.isShadowed(p, w.lights_buf.items[0], alctr) == false);
}

test "The shadow when an object is between the point and the light" {
    const alctr = std.testing.allocator;
    const w = getTestWorld(alctr);
    defer w.deinit();

    const p = Point.init(10, -10, 10);

    try expect(w.isShadowed(p, w.lights_buf.items[0], alctr) == true);
}

test "There is no shadow when an object is behind the light" {
    const alctr = std.testing.allocator;
    const w = getTestWorld(alctr);
    defer w.deinit();

    const p = Point.init(-20, 20, -20);

    try expect(w.isShadowed(p, w.lights_buf.items[0], alctr) == false);
}

test "There is no shadow when an object is behind the point" {
    const alctr = std.testing.allocator;
    const w = getTestWorld(alctr);
    defer w.deinit();

    const p = Point.init(-2, 2, -2);

    try expect(w.isShadowed(p, w.lights_buf.items[0], alctr) == false);
}
