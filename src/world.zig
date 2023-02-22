//! holds all volumes

const std = @import("std");

pub const World = struct {
    pub fn init(alctr: std.mem.Allocator) This {
        return .{
            .spheres_buf = std.ArrayList(vol.Sphere).init(alctr),
            .planes_buf = std.ArrayList(vol.Plane).init(alctr),
            .cubes_buf = std.ArrayList(vol.Cube).init(alctr),
            .cylinders_buf = std.ArrayList(vol.Cylinder).init(alctr),
            .cones_buf = std.ArrayList(vol.Cone).init(alctr),
            .triangles_buf = std.ArrayList(vol.Triangle).init(alctr),
            .lights_buf = std.ArrayList(PointLight).init(alctr),
        };
    }

    pub fn deinit(self: This) void {
        self.spheres_buf.deinit();
        self.planes_buf.deinit();
        self.cubes_buf.deinit();
        self.cylinders_buf.deinit();
        self.cones_buf.deinit();
        self.triangles_buf.deinit();
        self.lights_buf.deinit();
    }

    pub fn addVolume(self: *This, comptime T: type) struct { handle: VolumePtr, ptr: *T } {
        // create a new T and append it to the appropriate buf

        var buf = switch (T) {
            vol.Sphere => &self.spheres_buf,
            vol.Plane => &self.planes_buf,
            vol.Cube => &self.cubes_buf,
            vol.Cylinder => &self.cylinders_buf,
            vol.Cone => &self.cones_buf,
            vol.Triangle => &self.triangles_buf,
            else => unreachable,
        };

        switch (T) {
            vol.Triangle => buf.append(T.init(
                Point.init(0, 0, 0),
                Point.init(0, 1, 0),
                Point.init(1, 0, 0),
            )) catch unreachable,
            else => buf.append(T.init()) catch unreachable,
        }

        const last = @intCast(u16, buf.items.len - 1);
        const handle = switch (T) {
            vol.Sphere => VolumePtr{ .sphere_idx = last },
            vol.Plane => VolumePtr{ .plane_idx = last },
            vol.Cube => VolumePtr{ .cube_idx = last },
            vol.Cylinder => VolumePtr{ .cylinder_idx = last },
            vol.Cone => VolumePtr{ .cone_idx = last },
            vol.Triangle => VolumePtr{ .triangle_idx = last },
            else => unreachable,
        };

        return .{
            .handle = handle,
            .ptr = &buf.items[last],
        };
    }

    pub fn addLight(self: *This, comptime T: type) struct { handle: LightPtr, ptr: *T } {
        switch (T) {
            PointLight => {
                const point = @import("tuple.zig").Point.init(0, 0, 0);
                const color = @import("color.zig").Color.init(1, 1, 1);
                self.lights_buf.append(PointLight.init(point, color)) catch unreachable;
                const last = @intCast(u16, self.lights_buf.items.len - 1);
                return .{
                    .handle = LightPtr{ .light_idx = last },
                    .ptr = &self.lights_buf.items[last],
                };
            },
            else => unreachable,
        }
    }

    pub fn getVolume(self: This, vptr: VolumePtr, comptime T: type) *T {
        return switch (std.meta.activeTag(vptr)) {
            .sphere_idx => return &self.spheres_buf.items[vptr.idx()],
            .plane_idx => return &self.planes_buf.items[vptr.idx()],
            .cube_idx => return &self.cubes_buf.items[vptr.idx()],
            .cylinder_idx => return &self.cylinders.items[vptr.idx()],
            .cone_idx => return &self.cones_buf.items[vptr.idx()],
        };
    }

    pub fn getLight(self: This, comptime T: type, lptr: LightPtr) *T {
        return switch (T) {
            PointLight => return &self.lights_buf.items[lptr.light_idx],
            else => unreachable,
        };
    }

    /// Assumes that vol.Sphere defines what fields are available
    fn PropertyT(comptime name: []const u8) type {
        const fs = @typeInfo(vol.Sphere).Struct.fields;
        for (fs) |fd| {
            if (std.mem.eql(u8, name, fd.name)) {
                return *fd.type;
            }
        }
        unreachable;
    }

    /// Assumes that vol.Sphere defines what fields are available
    pub fn getProperty(
        self: This,
        volp: VolumePtr,
        comptime property: []const u8,
    ) PropertyT(property) {
        const i = switch (volp) {
            inline else => |i| i,
        };
        return switch (std.meta.activeTag(volp)) {
            .sphere_idx => return &@field(self.spheres_buf.items[i], property),
            .plane_idx => return &@field(self.planes_buf.items[i], property),
            .cube_idx => return &@field(self.cubes_buf.items[i], property),
            .cylinder_idx => return &@field(self.cylinders_buf.items[i], property),
            .cone_idx => return &@field(self.cones_buf.items[i], property),
            .triangle_idx => return &@field(self.triangles_buf.items[i], property),
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

        for (self.spheres_buf.items, 0..) |*ptr, i| {
            const vptr = VolumePtr{ .sphere_idx = @intCast(u16, i) };
            ixs.intersect(ptr.*, vptr, ray);
        }

        for (self.planes_buf.items, 0..) |*ptr, i| {
            const vptr = VolumePtr{ .plane_idx = @intCast(u16, i) };
            ixs.intersect(ptr.*, vptr, ray);
        }

        for (self.cubes_buf.items, 0..) |*ptr, i| {
            const vptr = VolumePtr{ .cube_idx = @intCast(u16, i) };
            ixs.intersect(ptr.*, vptr, ray);
        }

        for (self.cylinders_buf.items, 0..) |*ptr, i| {
            const vptr = VolumePtr{ .cylinder_idx = @intCast(u16, i) };
            ixs.intersect(ptr.*, vptr, ray);
        }

        for (self.cones_buf.items, 0..) |*ptr, i| {
            const vptr = VolumePtr{ .cone_idx = @intCast(u16, i) };
            ixs.intersect(ptr.*, vptr, ray);
        }

        for (self.triangles_buf.items, 0..) |*ptr, i| {
            const vptr = VolumePtr{ .triangle_idx = @intCast(u16, i) };
            ixs.intersect(ptr.*, vptr, ray);
        }
    }

    pub fn shadeHit(
        self: This,
        data: HitData,
        alctr: std.mem.Allocator,
        reflections_remaining: usize,
    ) Color {
        const tfm = self.getProperty(data.intersection.vptr, "transform").*;
        const mate = self.getProperty(data.intersection.vptr, "material").*;

        var color = Color.init(0, 0, 0);
        // surface color
        for (self.lights_buf.items) |l| {
            const in_shadow = self.isShadowed(data.over_point, l, alctr);
            color = color.plus(light.lighting(mate, tfm, l, data.point, data.eye_vector, data.normal_vector, in_shadow));
        }
        // reflected color
        var ref_color = self.reflectedColor(data, alctr, reflections_remaining);
        // refracted color
        var raf_color = self.refractedColor(data, alctr, reflections_remaining);

        if (mate.reflective > 0 and mate.transparency > 0) {
            // combine ref and raf with the schlick approx.
            const reflectance = data.schlick();
            return color.plus(ref_color.scaled(reflectance)).plus(raf_color.scaled(1 - reflectance));
        } else {
            return color.plus(ref_color).plus(raf_color);
        }
    }

    pub fn reflectedColor(self: This, data: HitData, alctr: std.mem.Allocator, reflections_remaining: usize) Color {
        if (reflections_remaining == 0) {
            return Color.init(0, 0, 0);
        }

        const mate = self.getProperty(data.intersection.vptr, "material");

        // material is not reflective
        if (mate.reflective == 0) return Color.init(0, 0, 0);

        // bounce ray and return what we find
        const reflected_ray = Ray.init(data.over_point, data.reflect_vector);
        const color = self.colorAt(reflected_ray, alctr, reflections_remaining - 1);

        return color.scaled(mate.reflective);
    }

    pub fn refractedColor(self: This, data: HitData, alctr: std.mem.Allocator, refractions_remaining: usize) Color {
        if (refractions_remaining == 0) {
            return Color.init(0, 0, 0);
        }

        const mate = self.getProperty(data.intersection.vptr, "material");

        // material is not refractive
        if (mate.transparency == 0) return Color.init(0, 0, 0);

        const n_ratio = data.n1 / data.n2;
        // cos(theta_i) is the same as the dot of the two vectors
        const cos_i = data.eye_vector.dot(data.normal_vector);
        // find sin(theta_t)^2 via trig identity
        const sin2_t = (n_ratio * n_ratio) * (1.0 - (cos_i * cos_i));

        // check for total internal reflection. based on snell's law (page 157)
        if (sin2_t > 1.0) return Color.init(0, 0, 0);

        // find cos(theta_t) via trig identity
        const cos_t = @sqrt(1.0 - sin2_t);

        const refr_direction = data.normal_vector.scaled(n_ratio * cos_i - cos_t).minus(data.eye_vector.scaled(n_ratio));

        const refr_ray = Ray.init(data.under_point, refr_direction);
        const color = self.colorAt(refr_ray, alctr, refractions_remaining - 1);

        // making sure to multiply by the transparency value to account for any opacity
        return color.scaled(mate.transparency);
    }

    fn getRefractiveIndicesFromBoundaries(self: This, bounds: Intersections.Boundary) struct {
        n1: f64,
        n2: f64,
    } {
        const n1 = if (bounds.lesser) |lsr|
            self.getProperty(lsr, "material").refractive_index
        else
            1.0;
        const n2 = if (bounds.greater) |gtr|
            self.getProperty(gtr, "material").refractive_index
        else
            1.0;
        return .{ .n1 = n1, .n2 = n2 };
    }

    pub fn colorAt(self: This, ray: Ray, alctr: std.mem.Allocator, reflections_remaining: usize) Color {
        var ixs = Intersections.init(alctr, .{});
        defer ixs.deinit();

        self.intersect(ray, &ixs);

        if (ixs.hit()) |hit| {
            const normal = switch (std.meta.activeTag(hit.vptr)) {
                .sphere_idx => self.spheres_buf.items[hit.idx()].normalAt(ray.position(hit.t)),
                .plane_idx => self.planes_buf.items[hit.idx()].normalAt(ray.position(hit.t)),
                .cube_idx => self.cubes_buf.items[hit.idx()].normalAt(ray.position(hit.t)),
                .cylinder_idx => self.cylinders_buf.items[hit.idx()].normalAt(ray.position(hit.t)),
                .cone_idx => self.cones_buf.items[hit.idx()].normalAt(ray.position(hit.t)),
                .triangle_idx => self.triangles_buf.items[hit.idx()].normalAt(ray.position(hit.t)),
            };

            const bounds = ixs.findBoundaryObjects(hit);
            const indices = self.getRefractiveIndicesFromBoundaries(bounds);
            const data = HitData.init(ray, hit, normal, indices.n1, indices.n2);
            return self.shadeHit(data, alctr, reflections_remaining);
        } else {
            return Color.init(0, 0, 0);
        }
    }

    pub fn render(self: This, cam: Camera, qan: *Qanvas, alctr: std.mem.Allocator) void {
        rdr.startRenderEngine(self, cam, qan, alctr);
    }

    pub fn numVolumes(self: This) usize {
        return self.spheres_buf.items.len + self.planes_buf.items.len + self.cubes_buf.items.len;
    }

    pub fn numLights(self: This) usize {
        return self.lights_buf.items.len;
    }

    spheres_buf: std.ArrayList(vol.Sphere),
    planes_buf: std.ArrayList(vol.Plane),
    cubes_buf: std.ArrayList(vol.Cube),
    cylinders_buf: std.ArrayList(vol.Cylinder),
    cones_buf: std.ArrayList(vol.Cone),
    triangles_buf: std.ArrayList(vol.Triangle),
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
    sphere_idx: u16,
    plane_idx: u16,
    cube_idx: u16,
    cylinder_idx: u16,
    cone_idx: u16,
    triangle_idx: u16,

    pub fn idx(self: VolumePtr) u16 {
        return switch (self) {
            inline else => |i| i,
        };
    }
};

pub const LightPtr = union(enum) {
    light_idx: u16,
};

fn getTestWorld(alctr: std.mem.Allocator) World {
    var w = World.init(alctr);

    { // default light
        var l = w.addLight(PointLight);
        l.ptr.position = Point.init(-10, 10, -10);
    }
    { // first sphere
        var s = w.addVolume(vol.Sphere);
        s.ptr.material.color_map = mat.FlatColor.init(Color.init(0.8, 1.0, 0.6));
        s.ptr.material.diffuse = 0.7;
        s.ptr.material.specular = 0.2;
    }
    { // second sphere
        var s = w.addVolume(vol.Sphere);
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

    _ = world.addVolume(vol.Sphere);
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
    try expect(w.spheres_buf.items[0].material.color_map.atC(0, 0, 0).equals(Color.init(0.8, 1.0, 0.6)));
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
    xs.order();
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
    const d = HitData.init(r, x, s.normalAt(r.position(x.t)), 1, 1);

    const c = w.shadeHit(d, alctr, 0);

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
    const d = HitData.init(r, x, s.normalAt(r.position(x.t)), 1, 1);

    const c = w.shadeHit(d, alctr, 0);

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
        _ = w.addVolume(vol.Sphere);
    }
    { // second sphere
        var s = w.addVolume(vol.Sphere);
        s.ptr.transform = trans.makeTranslation(0, 0, 10);
    }

    const r = Ray.init(Point.init(0, 0, 5), Vector.init(0, 0, 1));
    const s = w.spheres_buf.items[1];
    const x = Intersection.init(4, VolumePtr{ .sphere_idx = 1 });
    const d = HitData.init(r, x, s.normalAt(r.position(x.t)), 1, 1);

    const c = w.shadeHit(d, alctr, 0);

    // TODO book tests are imprecise
    try expect(c.equals(Color.init(0.1, 0.1, 0.1)));
}

test "The color when a ray misses" {
    const alctr = std.testing.allocator;

    const w = getTestWorld(alctr);
    defer w.deinit();

    const r = Ray.init(Point.init(0, 0, -5), Vector.init(0, 1, 0));

    const c = w.colorAt(r, alctr, 0);

    try expect(c.equals(Color.init(0, 0, 0)));
}

test "The color when a ray hits" {
    const alctr = std.testing.allocator;

    const w = getTestWorld(alctr);
    defer w.deinit();

    const r = Ray.init(Point.init(0, 0, -5), Vector.init(0, 0, 1));

    const c = w.colorAt(r, alctr, 0);

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

    const c = w.colorAt(r, alctr, 0);

    try expect(c.equals(w.spheres_buf.items[1].material.color_map.atC(0, 0, 0)));
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

test "The reflected color for a nonreflective material" {
    const alctr = std.testing.allocator;
    var w = getTestWorld(alctr);
    defer w.deinit();

    w.spheres_buf.items[1].material.ambient = 1;

    const i = Intersection.init(1, VolumePtr{ .sphere_idx = 1 });
    const r = Ray.init(Point.init(0, 0, 0), Vector.init(0, 0, 1));
    const data = HitData.init(r, i, w.spheres_buf.items[1].normalAt(r.position(i.t)), 1, 1);
    const color = w.reflectedColor(data, alctr, 1);

    try expect(color.equals(Color.init(0, 0, 0)));
}

test "The reflected color for a reflective material" {
    const alctr = std.testing.allocator;
    var w = getTestWorld(alctr);
    defer w.deinit();

    var pln = w.addVolume(vol.Plane);
    pln.ptr.transform = trans.makeTranslation(0, -1, 0);
    pln.ptr.material.reflective = 0.5;

    const r = Ray.init(
        Point.init(0, 0, -3),
        Vector.init(0, -@sqrt(2.0) / 2.0, @sqrt(2.0) / 2.0),
    );
    const i = Intersection.init(@sqrt(2.0), pln.handle);
    const data = HitData.init(r, i, pln.ptr.normalAt(r.position(i.t)), 1, 1);
    const color = w.reflectedColor(data, alctr, 1);

    // TODO tests in book are imprecise
    try expect(color.equalsTolerance(Color.init(0.19032, 0.2379, 0.14274), 100_000_000_000));
}

test "shadeHit color for a reflective material" {
    const alctr = std.testing.allocator;
    var w = getTestWorld(alctr);
    defer w.deinit();

    var pln = w.addVolume(vol.Plane);
    pln.ptr.transform = trans.makeTranslation(0, -1, 0);
    pln.ptr.material.reflective = 0.5;

    const r = Ray.init(
        Point.init(0, 0, -3),
        Vector.init(0, -@sqrt(2.0) / 2.0, @sqrt(2.0) / 2.0),
    );
    const i = Intersection.init(@sqrt(2.0), pln.handle);
    const data = HitData.init(r, i, pln.ptr.normalAt(r.position(i.t)), 1, 1);
    const color = w.shadeHit(data, alctr, 1);

    // TODO tests in book are imprecise
    try expect(color.equalsTolerance(Color.init(0.87677, 0.92436, 0.82918), 100_000_000_000));
}

test "colorAt doesn't get stuck in infinite recursion" {
    const alctr = std.testing.allocator;
    var w = getTestWorld(alctr);
    defer w.deinit();

    _ = w.addLight(PointLight);

    var pln1 = w.addVolume(vol.Plane);
    pln1.ptr.material.reflective = 1.0;
    pln1.ptr.transform = trans.makeTranslation(0, -1, 0);

    var pln2 = w.addVolume(vol.Plane);
    pln2.ptr.material.reflective = 1.0;
    pln2.ptr.transform = trans.makeTranslation(0, 1, 0);

    var ray = Ray.initC(0, 0, 0, 0, 1, 0);

    // test passes if this returns
    _ = w.colorAt(ray, alctr, 4);
}

test "The reflected color at the maximum recursive depth" {
    const alctr = std.testing.allocator;
    var w = getTestWorld(alctr);
    defer w.deinit();

    var pln = w.addVolume(vol.Plane);
    pln.ptr.transform = trans.makeTranslation(0, -1, 0);
    pln.ptr.material.reflective = 0.5;

    const r = Ray.init(
        Point.init(0, 0, -3),
        Vector.init(0, -@sqrt(2.0) / 2.0, @sqrt(2.0) / 2.0),
    );
    const i = Intersection.init(@sqrt(2.0), pln.handle);
    const data = HitData.init(r, i, pln.ptr.normalAt(r.position(i.t)), 1, 1);
    const color = w.reflectedColor(data, alctr, 0);

    try expect(color.equals(Color.init(0, 0, 0)));
}

test "The refracted color for a nonrefractive material" {
    const alctr = std.testing.allocator;
    var w = getTestWorld(alctr);
    defer w.deinit();

    const r = Ray.init(Point.init(0, 0, -5), Vector.init(0, 0, 1));
    const i_1 = Intersection.init(4, VolumePtr{ .sphere_idx = 0 });
    const i_2 = Intersection.init(6, VolumePtr{ .sphere_idx = 0 });

    var xs = Intersections.init(alctr, .{ i_1, i_2 });
    defer xs.deinit();

    const bnds = xs.findBoundaryObjects(i_1);

    const n1 = if (bnds.lesser) |lsr|
        w.getProperty(lsr, "material").refractive_index
    else
        1.0;

    const n2 = if (bnds.greater) |gtr|
        w.getProperty(gtr, "material").refractive_index
    else
        1.0;

    const data = HitData.init(
        r,
        i_1,
        w.spheres_buf.items[0].normalAt(r.position(i_1.t)),
        n1,
        n2,
    );

    const color = w.refractedColor(data, alctr, 1);

    try expect(color.equals(Color.init(0, 0, 0)));
}

test "The refracted color at max recursion" {
    const alctr = std.testing.allocator;
    var w = getTestWorld(alctr);
    defer w.deinit();

    prefab.toGlass(&w.spheres_buf.items[0].material);

    const r = Ray.init(Point.init(0, 0, -5), Vector.init(0, 0, 1));
    const i_1 = Intersection.init(4, VolumePtr{ .sphere_idx = 0 });
    const i_2 = Intersection.init(6, VolumePtr{ .sphere_idx = 0 });

    var xs = Intersections.init(alctr, .{ i_1, i_2 });
    defer xs.deinit();

    const bnds = xs.findBoundaryObjects(i_1);

    const n1 = if (bnds.lesser) |lsr|
        w.getProperty(lsr, "material").refractive_index
    else
        1.0;

    const n2 = if (bnds.greater) |gtr|
        w.getProperty(gtr, "material").refractive_index
    else
        1.0;

    const data = HitData.init(
        r,
        i_1,
        w.spheres_buf.items[0].normalAt(r.position(i_1.t)),
        n1,
        n2,
    );

    const color = w.refractedColor(data, alctr, 0);

    try expect(color.equals(Color.init(0, 0, 0)));
}

test "The refracted color under total internal reflection" {
    const alctr = std.testing.allocator;
    var w = getTestWorld(alctr);
    defer w.deinit();

    prefab.toGlass(&w.spheres_buf.items[0].material);

    const r = Ray.init(Point.init(0, 0, @sqrt(2.0) / 2.0), Vector.init(0, 1, 0));
    const i_1 = Intersection.init(-@sqrt(2.0) / 2.0, VolumePtr{ .sphere_idx = 0 });
    const i_2 = Intersection.init(@sqrt(2.0) / 2.0, VolumePtr{ .sphere_idx = 0 });

    var xs = Intersections.init(alctr, .{ i_1, i_2 });
    defer xs.deinit();

    // Since we're inside of the sphere, we look at i_2, not i_1
    const bnds = xs.findBoundaryObjects(i_2);

    const n1 = if (bnds.lesser) |lsr|
        w.getProperty(lsr, "material").refractive_index
    else
        1.0;

    const n2 = if (bnds.greater) |gtr|
        w.getProperty(gtr, "material").refractive_index
    else
        1.0;

    const data = HitData.init(
        r,
        i_2, // i_2 here too
        w.spheres_buf.items[0].normalAt(r.position(i_2.t)),
        n1,
        n2,
    );

    const color = w.refractedColor(data, alctr, 5);

    try expect(color.equals(Color.init(0, 0, 0)));
}

test "The refracted color" {
    const alctr = std.testing.allocator;
    var w = getTestWorld(alctr);
    defer w.deinit();

    w.spheres_buf.items[0].material.ambient = 1;
    w.spheres_buf.items[0].material.color_map = mat.TestColor.init();

    prefab.toGlass(&w.spheres_buf.items[1].material);

    const r = Ray.init(Point.init(0, 0, 0.1), Vector.init(0, 1, 0));

    const i_1 = Intersection.init(-0.9899, VolumePtr{ .sphere_idx = 0 });
    const i_2 = Intersection.init(-0.4899, VolumePtr{ .sphere_idx = 1 });
    const i_3 = Intersection.init(0.4899, VolumePtr{ .sphere_idx = 1 });
    const i_4 = Intersection.init(0.9899, VolumePtr{ .sphere_idx = 0 });

    var xs = Intersections.init(alctr, .{ i_1, i_2, i_3, i_4 });
    defer xs.deinit();

    const bnds = xs.findBoundaryObjects(i_3);

    const n1 = if (bnds.lesser) |lsr|
        w.getProperty(lsr, "material").refractive_index
    else
        1.0;

    const n2 = if (bnds.greater) |gtr|
        w.getProperty(gtr, "material").refractive_index
    else
        1.0;

    const data = HitData.init(
        r,
        i_3,
        w.spheres_buf.items[1].normalAt(r.position(i_3.t)),
        n1,
        n2,
    );

    const color = w.refractedColor(data, alctr, 5);

    // TODO book tests are imprecise
    errdefer print(color);
    try expect(color.equalsTolerance(Color.init(0, 0.99888, 0.04722), 100_000_000_000));
}

test "Shade hit handles refraction with a transparent material" {
    const alctr = std.testing.allocator;
    var w = getTestWorld(alctr);
    defer w.deinit();

    var flr = w.addVolume(vol.Plane);
    flr.ptr.material.transparency = 0.5;
    flr.ptr.material.refractive_index = 1.5;
    flr.ptr.transform = trans.makeTranslation(0, -1, 0);

    var ball = w.addVolume(vol.Sphere);
    ball.ptr.material.color_map = mat.FlatColor.init(Color.init(1, 0, 0));
    ball.ptr.material.ambient = 0.5;
    ball.ptr.transform = trans.makeTranslation(0, -3.5, -0.5);

    const r = Ray.init(Point.init(0, 0, -3), Vector.init(0, -@sqrt(2.0) / 2.0, @sqrt(2.0) / 2.0));

    const i_1 = Intersection.init(@sqrt(2.0), VolumePtr{ .plane_idx = 0 });

    var xs = Intersections.init(alctr, .{i_1});
    defer xs.deinit();

    const bnds = xs.findBoundaryObjects(i_1);

    const n1 = if (bnds.lesser) |lsr|
        w.getProperty(lsr, "material").refractive_index
    else
        1.0;

    const n2 = if (bnds.greater) |gtr|
        w.getProperty(gtr, "material").refractive_index
    else
        1.0;

    const data = HitData.init(
        r,
        i_1,
        w.planes_buf.items[0].normalAt(r.position(i_1.t)),
        n1,
        n2,
    );

    const color = w.shadeHit(data, alctr, 5);

    // TODO book tests are imprecise
    errdefer print(color);
    try expect(color.equalsTolerance(Color.init(0.93643, 0.68643, 0.68643), 100_000_000_000));
}

test "Shade hit handles a reflective and refractive surface and takes reflectance into account" {
    const alctr = std.testing.allocator;
    var w = getTestWorld(alctr);
    defer w.deinit();

    var flr = w.addVolume(vol.Plane);
    flr.ptr.material.transparency = 0.5;
    flr.ptr.material.refractive_index = 1.5;
    flr.ptr.material.reflective = 0.5;
    flr.ptr.transform = trans.makeTranslation(0, -1, 0);

    var ball = w.addVolume(vol.Sphere);
    ball.ptr.material.color_map = mat.FlatColor.init(Color.init(1, 0, 0));
    ball.ptr.material.ambient = 0.5;
    ball.ptr.transform = trans.makeTranslation(0, -3.5, -0.5);

    const r = Ray.init(Point.init(0, 0, -3), Vector.init(0, -@sqrt(2.0) / 2.0, @sqrt(2.0) / 2.0));

    const i_1 = Intersection.init(@sqrt(2.0), VolumePtr{ .plane_idx = 0 });

    var xs = Intersections.init(alctr, .{i_1});
    defer xs.deinit();

    const bnds = xs.findBoundaryObjects(i_1);

    const n1 = if (bnds.lesser) |lsr|
        w.getProperty(lsr, "material").refractive_index
    else
        1.0;

    const n2 = if (bnds.greater) |gtr|
        w.getProperty(gtr, "material").refractive_index
    else
        1.0;

    const data = HitData.init(
        r,
        i_1,
        w.planes_buf.items[0].normalAt(r.position(i_1.t)),
        n1,
        n2,
    );

    const color = w.shadeHit(data, alctr, 5);

    // TODO book tests are imprecise
    errdefer print(color);
    try expect(color.equalsTolerance(Color.init(0.93392, 0.69643, 0.69243), 100_000_000_000));
}

test "getProperty" {
    const alctr = std.testing.allocator;
    var w = World.init(alctr);
    defer w.deinit();

    var sphere = w.addVolume(vol.Sphere);
    sphere.ptr.transform = trans.makeTranslation(0, 10, 0);
    sphere.ptr.material.ambient = 0.75;

    var plane = w.addVolume(vol.Plane);
    plane.ptr.transform = trans.makeTranslation(0, -10, 0);
    plane.ptr.material.ambient = 0.25;

    try expect(w.getProperty(sphere.handle, "material").ambient == 0.75);
    try expect(w.getProperty(sphere.handle, "transform").t.equals(trans.makeTranslation(0, 10, 0).t));
    try expect(w.getProperty(plane.handle, "material").ambient == 0.25);
    try expect(w.getProperty(plane.handle, "transform").t.equals(trans.makeTranslation(0, -10, 0).t));
}

test "Finding entry and exit volumes at various intersections" {
    const alctr = std.testing.allocator;
    var w = World.init(alctr);
    defer w.deinit();

    // diagram on page 151
    const s1 = blk: {
        var s = w.addVolume(vol.Sphere);
        prefab.toGlass(&s.ptr.material);
        s.ptr.material.refractive_index = 1.5;
        s.ptr.transform = trans.makeScaling(2, 2, 2);
        break :blk s.handle;
    };
    const s2 = blk: {
        var s = w.addVolume(vol.Sphere);
        prefab.toGlass(&s.ptr.material);
        s.ptr.material.refractive_index = 2.0;
        s.ptr.transform = trans.makeTranslation(0, 0, -0.25);
        break :blk s.handle;
    };
    const s3 = blk: {
        var s = w.addVolume(vol.Sphere);
        prefab.toGlass(&s.ptr.material);
        s.ptr.material.refractive_index = 2.5;
        s.ptr.transform = trans.makeTranslation(0, 0, 0.25);
        break :blk s.handle;
    };
    var xs = Intersections.init(alctr, .{
        Intersection.init(4.75, s2),
        Intersection.init(5.25, s3),
        Intersection.init(6, s1),
        Intersection.init(2, s1),
        Intersection.init(2.75, s2),
        Intersection.init(3.25, s3),
    });
    defer xs.deinit();

    {
        const boundaries = xs.findBoundaryObjects(xs.ixs.items[0]);
        errdefer print(boundaries);
        try expect(boundaries.lesser == null);
        try expect(boundaries.greater != null);
        try expect(std.meta.eql(boundaries.greater.?, s1));
    }
    {
        const boundaries = xs.findBoundaryObjects(xs.ixs.items[1]);
        try expect(std.meta.eql(boundaries.lesser.?, s1));
        try expect(std.meta.eql(boundaries.greater.?, s2));

        const n1 = w.getProperty(boundaries.lesser.?, "material").refractive_index;
        const n2 = w.getProperty(boundaries.greater.?, "material").refractive_index;
        try expect(n1 == 1.5);
        try expect(n2 == 2.0);
    }
    {
        const boundaries = xs.findBoundaryObjects(xs.ixs.items[2]);
        try expect(std.meta.eql(boundaries.lesser.?, s2));
        try expect(std.meta.eql(boundaries.greater.?, s3));
    }
    {
        // note this result: we are passing through the edge of s2.
        // however, before we reached this edge, we entered s3 (index 2 above).
        // s3 completely surrounds this exit point on both inside and outside,
        // so lesser and greater *both* are s3.
        const boundaries = xs.findBoundaryObjects(xs.ixs.items[3]);
        errdefer print(boundaries);
        try expect(std.meta.eql(boundaries.lesser.?, s3));
        try expect(std.meta.eql(boundaries.greater.?, s3));
    }
    {
        const boundaries = xs.findBoundaryObjects(xs.ixs.items[4]);
        try expect(std.meta.eql(boundaries.lesser.?, s3));
        try expect(std.meta.eql(boundaries.greater.?, s1));
    }
    {
        const boundaries = xs.findBoundaryObjects(xs.ixs.items[5]);
        try expect(std.meta.eql(boundaries.lesser.?, s1));
        try expect(boundaries.greater == null);
    }
}

const Tuple = @import("tuple.zig").Tuple;
const Point = @import("tuple.zig").Point;
const Vector = @import("tuple.zig").Vector;
const Color = @import("color.zig").Color;

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
const vol = @import("volume.zig");
const mat = @import("material.zig");
const prefab = @import("prefab.zig");

const expect = std.testing.expect;
const print = @import("u.zig").print;
