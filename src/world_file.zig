//! create and populate a world from a text file

pub var error_message: ?[]const u8 = null;
pub var error_offending_line: ?[]const u8 = null;

const FileContents = struct {
    world: World,
    camera: Camera,
};

pub fn parseWorldFile(filename: []const u8, alctr: std.mem.Allocator) !FileContents {
    var file = try std.fs.cwd().openFile(filename, .{});
    defer file.close();

    // a gigabyte
    var txt = try file.readToEndAlloc(alctr, 1024 * 1024 * 1024);
    defer alctr.free(txt);

    return try parseWorldText(txt, alctr);
}

fn findSectionBounds(txt: []const u8) struct { min: ?usize, max: usize } {
    const keys = [_][]const u8{
        "CAMERA",
        "POINTLIGHT",
        "SPHERE",
        "PLANE",
        "CUBE",
        "CYLINDER",
        "CONE",
    };

    var keys_idxs = [_]?usize{null} ** keys.len;

    const min = blk: {
        for (keys) |key, i| {
            keys_idxs[i] = std.mem.indexOf(u8, txt, key);
        }

        const lt = struct {
            fn lt(_: void, a: ?usize, b: ?usize) bool {
                if (a == null) return false;
                return a.? < b orelse std.math.maxInt(usize);
            }
        }.lt;

        std.sort.insertionSort(?usize, &keys_idxs, {}, lt);

        break :blk keys_idxs[0].?;
    };

    const max = blk: {
        for (keys) |key, i| {
            keys_idxs[i] = std.mem.indexOf(u8, txt[min + 1 ..], key);
            if (keys_idxs[i] != null) keys_idxs[i].? += 1;
        }

        const lt = struct {
            fn lt(_: void, a: ?usize, b: ?usize) bool {
                if (a == null) return false;
                return a.? < b orelse std.math.maxInt(usize);
            }
        }.lt;

        std.sort.insertionSort(?usize, &keys_idxs, {}, lt);

        break :blk keys_idxs[0] orelse txt.len;
    };

    return .{ .min = min, .max = max };
}

pub fn parseWorldText(txt: []const u8, alctr: std.mem.Allocator) !FileContents {
    var world = World.init(alctr);
    var camera: ?Camera = null;

    var remaining = txt;
    var chunk_idxs = findSectionBounds(remaining);

    while (chunk_idxs.min) |min| {
        const chunk = remaining[min..chunk_idxs.max];
        const header = std.mem.sliceTo(chunk, '\n');
        if (std.mem.startsWith(u8, header, "CAMERA")) {
            camera = try parseCamera(chunk);
        } else if (std.mem.startsWith(u8, header, "POINTLIGHT")) {
            var world_light = world.addLight(PointLight);
            var parsed_light = try parsePointLight(chunk);
            world_light.ptr.* = parsed_light;
        } else if (std.mem.startsWith(u8, header, "CONE")) {
            var world_cone = world.addVolume(Cone);
            var parsed_cone = try parseCone(chunk);
            world_cone.ptr.* = parsed_cone;
        } else if (std.mem.startsWith(u8, header, "SPHERE")) {
            var world_sphere = world.addVolume(Sphere);
            var parsed_sphere = try parseSphere(chunk);
            world_sphere.ptr.* = parsed_sphere;
        } else if (std.mem.startsWith(u8, header, "PLANE")) {
            var world_plane = world.addVolume(Plane);
            var parsed_plane = try parsePlane(chunk);
            world_plane.ptr.* = parsed_plane;
        } else if (std.mem.startsWith(u8, header, "CUBE")) {
            var world_cube = world.addVolume(Cube);
            var parsed_cube = try parseCube(chunk);
            world_cube.ptr.* = parsed_cube;
        } else if (std.mem.startsWith(u8, header, "CYLINDER")) {
            var world_cylinder = world.addVolume(Cylinder);
            var parsed_cylinder = try parseCylinder(chunk);
            world_cylinder.ptr.* = parsed_cylinder;
        }
        std.debug.print("+ {s}\n", .{header});

        if (chunk_idxs.max == remaining.len) break;

        remaining = remaining[chunk_idxs.max..];
        chunk_idxs = findSectionBounds(remaining);
    }

    if (camera == null) return unimplementedError();

    return .{ .world = world, .camera = camera.? };
}

pub const WorldParseError = error{
    IndescribableOffense,
};

fn unimplementedError() WorldParseError {
    error_message = "Unknown";
    error_offending_line = "Unknown";
    return WorldParseError.IndescribableOffense;
}

fn setAndReturnError(e: WorldParseError, message: []const u8, offending: []const u8) WorldParseError {
    error_message = message;
    error_offending_line = offending;
    return e;
}

fn parseKeyValue(comptime KeyType: type, txt: []const u8) !KeyType {
    var tokens = std.mem.tokenize(u8, txt, " \t");
    _ = tokens.next(); // skip key
    const value = tokens.rest();

    // structs
    switch (KeyType) {
        @"3Tuple" => return try parse3Tuple(value),
        Angle => return try parseAngle(value),
        else => {},
    }

    // builtin types
    const tinfo = @typeInfo(KeyType);
    switch (tinfo) {
        .Int => return try std.fmt.parseInt(@Type(tinfo), value, 10),
        .Float => return try std.fmt.parseFloat(@Type(tinfo), value),
        .Bool => {
            if (std.mem.eql(u8, value, "true")) return true;
            if (std.mem.eql(u8, value, "false")) return false;
            unreachable;
        },
        else => unreachable,
    }
}

const Angle = struct {
    radians: f64,
};

fn parseAngle(txt: []const u8) !Angle {
    // angles can be floats or a special pi/<float> syntax
    const radians =
        std.fmt.parseFloat(f64, txt) catch |e| switch (e) {
        error.InvalidCharacter => blk: {
            if (std.mem.startsWith(u8, txt, "pi/")) {
                const number = txt[3..];
                const denominator = std.fmt.parseFloat(f64, number) catch return unimplementedError();
                break :blk PI / denominator;
            } else {
                return unimplementedError();
            }
        },
    };

    return Angle{ .radians = radians };
}

const @"3Tuple" = struct {
    a: f64,
    b: f64,
    c: f64,
};

fn parse3Tuple(txt: []const u8) !@"3Tuple" {
    // like (1.0 2 -2.3)

    errdefer tprint(txt);
    const closep = std.mem.indexOf(u8, txt, ")") orelse return unimplementedError();
    var tup_str = txt[0 .. closep + 1];
    if (tup_str[0] != '(') return unimplementedError();
    if (tup_str[tup_str.len - 1] != ')') return unimplementedError();

    tup_str = tup_str[1 .. tup_str.len - 1];
    if (tup_str.len == 0) return unimplementedError();

    var values = std.mem.tokenize(u8, tup_str, ",");
    const a_str = values.next() orelse return unimplementedError();
    const b_str = values.next() orelse return unimplementedError();
    const c_str = values.next() orelse return unimplementedError();

    return @"3Tuple"{
        .a = try std.fmt.parseFloat(f64, a_str),
        .b = try std.fmt.parseFloat(f64, b_str),
        .c = try std.fmt.parseFloat(f64, c_str),
    };
}

fn parseCamera(txt: []const u8) !Camera {
    var scale: ?f64 = null;
    var width: ?i64 = null;
    var height: ?i64 = null;
    var fov: ?f64 = null;
    var from: ?@"3Tuple" = null;
    var to: ?@"3Tuple" = null;
    var up: ?@"3Tuple" = null;

    {
        var lines = std.mem.tokenize(u8, txt, "\n");
        while (lines.next()) |line| {
            if (std.mem.startsWith(u8, line, "scale")) {
                scale = try parseKeyValue(f64, line);
            }
            if (std.mem.startsWith(u8, line, "width")) {
                width = try parseKeyValue(i64, line);
            }
            if (std.mem.startsWith(u8, line, "height")) {
                height = try parseKeyValue(i64, line);
            }
            if (std.mem.startsWith(u8, line, "fov")) {
                const angle = try parseKeyValue(Angle, line);
                fov = angle.radians;
            }
            if (std.mem.startsWith(u8, line, "from")) {
                from = try parseKeyValue(@"3Tuple", line);
            }
            if (std.mem.startsWith(u8, line, "to")) {
                to = try parseKeyValue(@"3Tuple", line);
            }
            if (std.mem.startsWith(u8, line, "up")) {
                up = try parseKeyValue(@"3Tuple", line);
            }
        }
    }

    if (width == null or height == null or fov == null or
        from == null or to == null or up == null)
        return unimplementedError();

    scale = scale orelse 1;
    const fwidth = @intToFloat(f64, width.?) * scale.?;
    const fheight = @intToFloat(f64, height.?) * scale.?;

    var camera = Camera.init(@floatToInt(i64, fwidth), @floatToInt(i64, fheight), fov.?);
    camera.transform = trans.makeView(
        Point.init(from.?.a, from.?.b, from.?.c),
        Point.init(to.?.a, to.?.b, to.?.c),
        Vector.init(up.?.a, up.?.b, up.?.c),
    );

    return camera;
}

fn parsePointLight(txt: []const u8) !PointLight {
    var position: ?@"3Tuple" = null;
    var intensity: ?@"3Tuple" = null;

    {
        var lines = std.mem.tokenize(u8, txt, "\n");
        while (lines.next()) |line| {
            if (std.mem.startsWith(u8, line, "position")) {
                position = try parseKeyValue(@"3Tuple", line);
            }
            if (std.mem.startsWith(u8, line, "intensity")) {
                intensity = try parseKeyValue(@"3Tuple", line);
            }
        }
    }

    if (position == null or intensity == null) return unimplementedError();

    return PointLight.init(
        Point.init(position.?.a, position.?.b, position.?.c),
        Color.init(intensity.?.a, intensity.?.b, intensity.?.c),
    );
}

fn parseAndApplyMaterial(txt: []const u8, material: *Material) !void {
    var tokens = std.mem.tokenize(u8, txt, " \t");
    _ = tokens.next(); // consume "material"

    const first = tokens.next() orelse return unimplementedError();

    if (std.mem.eql(u8, first, "3dchecker")) {
        const a_value_string = tokens.next() orelse return unimplementedError();
        const b_value_string = tokens.next() orelse return unimplementedError();
        const a_value = try parse3Tuple(a_value_string);
        const b_value = try parse3Tuple(b_value_string);

        material.color_map = mate.ThreeDCheckedColor.init(
            Color.init(a_value.a, a_value.b, a_value.c),
            Color.init(b_value.a, b_value.b, b_value.c),
        );
    } else if (std.mem.eql(u8, first, "ring")) {
        const a_value_string = tokens.next() orelse return unimplementedError();
        const b_value_string = tokens.next() orelse return unimplementedError();
        const a_value = try parse3Tuple(a_value_string);
        const b_value = try parse3Tuple(b_value_string);

        material.color_map = mate.RingColor.init(
            Color.init(a_value.a, a_value.b, a_value.c),
            Color.init(b_value.a, b_value.b, b_value.c),
        );
    } else if (std.mem.eql(u8, first, "flat")) {
        const a_value_string = tokens.next() orelse return unimplementedError();
        const a_value = try parse3Tuple(a_value_string);

        material.color_map = mate.FlatColor.init(
            Color.init(a_value.a, a_value.b, a_value.c),
        );
    } else if (std.mem.eql(u8, first, "diffuse")) {
        const value_string = tokens.next() orelse return unimplementedError();
        const value = try std.fmt.parseFloat(f64, value_string);

        material.diffuse = value;
    } else if (std.mem.eql(u8, first, "ambient")) {
        const value_string = tokens.next() orelse return unimplementedError();
        const value = try std.fmt.parseFloat(f64, value_string);

        material.ambient = value;
    } else if (std.mem.eql(u8, first, "specular")) {
        const value_string = tokens.next() orelse return unimplementedError();
        const value = try std.fmt.parseFloat(f64, value_string);

        material.specular = value;
    } else if (std.mem.eql(u8, first, "reflective")) {
        const value_string = tokens.next() orelse return unimplementedError();
        const value = try std.fmt.parseFloat(f64, value_string);

        material.reflective = value;
    } else if (std.mem.eql(u8, first, "transparency")) {
        const value_string = tokens.next() orelse return unimplementedError();
        const value = try std.fmt.parseFloat(f64, value_string);

        material.transparency = value;
    } else if (std.mem.eql(u8, first, "refractive")) {
        const value_string = tokens.next() orelse return unimplementedError();
        const value = try std.fmt.parseFloat(f64, value_string);

        material.refractive_index = value;
    } else if (std.mem.eql(u8, first, "transform")) {
        const tform_idx = std.mem.indexOf(u8, txt, "transform") orelse return unimplementedError();
        const substring = txt[tform_idx..];
        try parseAndApplyTransform(substring, &material.transform);
    } else {
        return unimplementedError();
    }
}

fn parseAndApplyTransform(txt: []const u8, transform: *Transform) !void {
    var tokens = std.mem.tokenize(u8, txt, " \t");
    _ = tokens.next(); // consume "transform"

    const first = tokens.next() orelse return unimplementedError();

    if (std.mem.eql(u8, first, "translate")) {
        const value_string = tokens.next() orelse return unimplementedError();
        const value = try parse3Tuple(value_string);

        const modifier = trans.makeTranslation(value.a, value.b, value.c);
        transform.* = modifier.mult(transform.*);
    } else if (std.mem.eql(u8, first, "scale")) {
        const value_string = tokens.next() orelse return unimplementedError();
        const value = try parse3Tuple(value_string);

        const modifier = trans.makeScaling(value.a, value.b, value.c);
        transform.* = modifier.mult(transform.*);
    } else if (std.mem.eql(u8, first, "rotate")) {
        const axis_string = tokens.next() orelse return unimplementedError();
        const value_string = tokens.next() orelse return unimplementedError();

        const modifier = blk: {
            if (std.mem.eql(u8, axis_string, "x")) {
                const angle = try parseAngle(value_string);
                const rads = angle.radians;
                break :blk trans.makeRotationX(rads);
            } else if (std.mem.eql(u8, axis_string, "y")) {
                const angle = try parseAngle(value_string);
                const rads = angle.radians;
                break :blk trans.makeRotationY(rads);
            } else if (std.mem.eql(u8, axis_string, "z")) {
                const angle = try parseAngle(value_string);
                const rads = angle.radians;
                break :blk trans.makeRotationZ(rads);
            } else {
                return unimplementedError();
            }
        };

        transform.* = modifier.mult(transform.*);
    }
}

fn parseCube(txt: []const u8) !Cube {
    var cube = Cube.init();

    {
        var lines = std.mem.tokenize(u8, txt, "\n");
        while (lines.next()) |line| {
            if (std.mem.startsWith(u8, line, "material")) {
                try parseAndApplyMaterial(line, &cube.material);
            } else if (std.mem.startsWith(u8, line, "transform")) {
                try parseAndApplyTransform(line, &cube.transform);
            }
        }
    }

    return cube;
}

fn parsePlane(txt: []const u8) !Plane {
    var plane = Plane.init();

    {
        var lines = std.mem.tokenize(u8, txt, "\n");
        while (lines.next()) |line| {
            if (std.mem.startsWith(u8, line, "material")) {
                try parseAndApplyMaterial(line, &plane.material);
            } else if (std.mem.startsWith(u8, line, "transform")) {
                try parseAndApplyTransform(line, &plane.transform);
            }
        }
    }

    return plane;
}

fn parseSphere(txt: []const u8) !Sphere {
    var sphere = Sphere.init();

    {
        var lines = std.mem.tokenize(u8, txt, "\n");
        while (lines.next()) |line| {
            if (std.mem.startsWith(u8, line, "material")) {
                try parseAndApplyMaterial(line, &sphere.material);
            } else if (std.mem.startsWith(u8, line, "transform")) {
                try parseAndApplyTransform(line, &sphere.transform);
            }
        }
    }

    return sphere;
}

fn parseCylinder(txt: []const u8) !Cylinder {
    var cylinder = Cylinder.init();

    {
        var lines = std.mem.tokenize(u8, txt, "\n");
        while (lines.next()) |line| {
            if (std.mem.startsWith(u8, line, "material")) {
                try parseAndApplyMaterial(line, &cylinder.material);
            } else if (std.mem.startsWith(u8, line, "transform")) {
                try parseAndApplyTransform(line, &cylinder.transform);
            } else if (std.mem.startsWith(u8, line, "closed")) {
                const value = try parseKeyValue(bool, line);
                cylinder.closed = value;
            } else if (std.mem.startsWith(u8, line, "length")) {
                const value = try parseKeyValue(f64, line);
                cylinder.length = value;
            }
        }
    }

    return cylinder;
}

fn parseCone(txt: []const u8) !Cone {
    var cone = Cone.init();

    {
        var lines = std.mem.tokenize(u8, txt, "\n");
        while (lines.next()) |line| {
            if (std.mem.startsWith(u8, line, "material")) {
                try parseAndApplyMaterial(line, &cone.material);
            } else if (std.mem.startsWith(u8, line, "transform")) {
                try parseAndApplyTransform(line, &cone.transform);
            } else if (std.mem.startsWith(u8, line, "closed")) {
                const value = try parseKeyValue(bool, line);
                cone.closed = value;
            } else if (std.mem.startsWith(u8, line, "min")) {
                const value = try parseKeyValue(f64, line);
                cone.min = value;
            } else if (std.mem.startsWith(u8, line, "max")) {
                const value = try parseKeyValue(f64, line);
                cone.max = value;
            }
        }
    }

    return cone;
}

test "parse camera" {
    const txt =
        \\CAMERA
        \\scale 2
        \\width 2400
        \\height 1200
        \\fov pi/2.5
        \\from (-0.2,3.5,-10)
        \\to (-0.2,2,0)
        \\up (0,1,0)
    ;

    const camera = try parseCamera(txt);

    try exEq(@as(i64, 4800), camera.width);
    try exEq(@as(i64, 2400), camera.height);
    try exEq(@as(f64, PI / 2.5), camera.fov);

    // transform is derived from from, to, and up
    errdefer tprint(camera.transform.t);
    const view =
        trans.makeView(
        Point.init(-0.2, 3.5, -10),
        Point.init(-0.2, 2, 0),
        Vector.init(0, 1, 0),
    ).t;
    errdefer tprint(view);

    try ex(camera.transform.t.equals(view));
}

test "parse point light" {
    const txt =
        \\POINTLIGHT
        \\position (-4,10,-5)
        \\intensity (1,1,1)
    ;

    const light = try parsePointLight(txt);

    try ex(light.position.equals(Point.init(-4, 10, -5)));
    try ex(light.intensity.equals(Color.init(1, 1, 1)));
}

test "parse cube (1)" {
    const txt =
        \\CUBE
        \\material 3dchecker (0.6,0.2,0.2) (0.2,0.2,0.6)
        \\material diffuse 0.2
        \\material ambient 0.1
        \\material reflective 0.9
        \\material transparency 0.9
        \\material refractive 1.52
        \\transform rotate y pi/-12
        \\transform scale (0.7,2,0.7)
        \\transform translate (3,2,0)
    ;

    const cube = try parseCube(txt);

    try ex(std.meta.activeTag(cube.material.color_map) == .threedchecked);
    try ex(cube.material.color_map.threedchecked.a.equals(Color.init(0.6, 0.2, 0.2)));
    try ex(cube.material.color_map.threedchecked.b.equals(Color.init(0.2, 0.2, 0.6)));
    try exEq(@as(f64, 0.2), cube.material.diffuse);
    try exEq(@as(f64, 0.1), cube.material.ambient);
    try exEq(@as(f64, 0.9), cube.material.reflective);
    try exEq(@as(f64, 0.9), cube.material.transparency);
    try exEq(@as(f64, 1.52), cube.material.refractive_index);

    // chain applies transforms in reverse order
    var t = (trans.Transform{}).chain(.{
        trans.makeTranslation(3, 2, 0),
        trans.makeScaling(0.7, 2, 0.7),
        trans.makeRotationY(PI / -12.0),
    });
    try ex(cube.transform.t.equals(t.t));
}

test "parse cube (2)" {
    const txt =
        \\CUBE
        \\material 3dchecker (0.6,0.2,0.2) (0.2,0.2,0.6)
        \\material diffuse 0.1
        \\material ambient 0
        \\material reflective 1
        \\transform rotate y pi/4.0
        \\transform translate (0,1,0)
    ;

    const cube = try parseCube(txt);

    try ex(std.meta.activeTag(cube.material.color_map) == .threedchecked);
    try ex(cube.material.color_map.threedchecked.a.equals(Color.init(0.6, 0.2, 0.2)));
    try ex(cube.material.color_map.threedchecked.b.equals(Color.init(0.2, 0.2, 0.6)));
    try exEq(@as(f64, 0.1), cube.material.diffuse);
    try exEq(@as(f64, 0), cube.material.ambient);
    try exEq(@as(f64, 1), cube.material.reflective);

    // chain applies transforms in reverse order
    var t = (trans.Transform{}).chain(.{
        trans.makeTranslation(0, 1, 0),
        trans.makeRotationY(PI / 4.0),
    });
    try ex(cube.transform.t.equals(t.t));
}

test "parse cube (3)" {
    // intentional extra whitespace on this one
    const txt =
        \\CUBE
        \\material 3dchecker (0.6,0.2,0.2) (0.2,0.2,0.6)
        \\transform scale (1.2,0.7,1.2)
        \\transform translate (-3,0.7,0)
        \\
        \\
        \\
    ;

    const cube = try parseCube(txt);

    try ex(std.meta.activeTag(cube.material.color_map) == .threedchecked);
    try ex(cube.material.color_map.threedchecked.a.equals(Color.init(0.6, 0.2, 0.2)));
    try ex(cube.material.color_map.threedchecked.b.equals(Color.init(0.2, 0.2, 0.6)));

    // chain applies transforms in reverse order
    var t = (trans.Transform{}).chain(.{
        trans.makeTranslation(-3, 0.7, 0),
        trans.makeScaling(1.2, 0.7, 1.2),
    });
    try ex(cube.transform.t.equals(t.t));
}

test "parse cube (4)" {
    // intentional blank line in this one
    const txt =
        \\CUBE
        \\material ring (0.6,0.2,0.2) (0.2,0.2,0.6)
        \\transform scale (0.3,0.3,0.3)
        \\
        \\transform rotate y 42
        \\transform translate (-3.3,1.7,0.3)
    ;

    const cube = try parseCube(txt);

    try ex(std.meta.activeTag(cube.material.color_map) == .ring);
    try ex(cube.material.color_map.ring.a.equals(Color.init(0.6, 0.2, 0.2)));
    try ex(cube.material.color_map.ring.b.equals(Color.init(0.2, 0.2, 0.6)));

    // chain applies transforms in reverse order
    var t = (trans.Transform{}).chain(.{
        trans.makeTranslation(-3.3, 1.7, 0.3),
        trans.makeRotationY(42),
        trans.makeScaling(0.3, 0.3, 0.3),
    });
    try ex(cube.transform.t.equals(t.t));
}

test "parse plane" {
    const txt =
        \\PLANE
        \\material 3dchecker (0.5,0.5,0.5) (0.25,0.25,0.25)
        \\material specular 0.1
        \\material diffuse 1
        \\material reflective 0.12
        \\material transform translate (0.5,0,0.5)
        \\material transform scale (8,8,8)
    ;

    const plane = try parsePlane(txt);

    try ex(std.meta.activeTag(plane.material.color_map) == .threedchecked);
    try ex(plane.material.color_map.threedchecked.a.equals(Color.init(0.5, 0.5, 0.5)));
    try ex(plane.material.color_map.threedchecked.b.equals(Color.init(0.25, 0.25, 0.25)));
    try exEq(@as(f64, 0.1), plane.material.specular);
    try exEq(@as(f64, 1), plane.material.diffuse);
    try exEq(@as(f64, 0.12), plane.material.reflective);

    // chain applies transforms in reverse order
    const mt = (trans.Transform{}).chain(.{
        trans.makeScaling(8, 8, 8),
        trans.makeTranslation(0.5, 0, 0.5),
    });
    try ex(plane.material.transform.t.equals(mt.t));
}

test "parse sphere (1)" {
    const txt =
        \\SPHERE
        \\material flat (0.2,0.2,0.7)
        \\transform translate (0,3,0)
    ;

    const sphere = try parseSphere(txt);

    try ex(std.meta.activeTag(sphere.material.color_map) == .flat);
    try ex(sphere.material.color_map.flat.color.equals(Color.init(0.2, 0.2, 0.7)));
    try ex(sphere.transform.t.equals(trans.makeTranslation(0, 3, 0).t));
}

test "parse sphere (2)" {
    const txt =
        \\SPHERE
        \\material flat (0.2,0.2,0.7)
        \\transform translate (3,1,-2.5)
    ;

    const sphere = try parseSphere(txt);

    try ex(std.meta.activeTag(sphere.material.color_map) == .flat);
    try ex(sphere.material.color_map.flat.color.equals(Color.init(0.2, 0.2, 0.7)));
    try ex(sphere.transform.t.equals(trans.makeTranslation(3, 1, -2.5).t));
}

test "parse cylinder (1)" {
    const txt =
        \\CYLINDER
        \\material 3dchecker (0.6,0.2,0.2) (0.2,0.2,0.6)
        \\closed true
        \\length 1.5
        \\transform scale (0.3,1,0.3)
        \\transform rotate y 32
        \\transform translate (-2.2,2.15,-0.8)
    ;

    const cylinder = try parseCylinder(txt);

    try ex(std.meta.activeTag(cylinder.material.color_map) == .threedchecked);
    try ex(cylinder.material.color_map.threedchecked.a.equals(Color.init(0.6, 0.2, 0.2)));
    try ex(cylinder.material.color_map.threedchecked.b.equals(Color.init(0.2, 0.2, 0.6)));
    try exEq(@as(bool, true), cylinder.closed);
    try exEq(@as(f64, 1.5), cylinder.length);

    // chain applies transforms in reverse order
    var t = (trans.Transform{}).chain(.{
        trans.makeTranslation(-2.2, 2.15, -0.8),
        trans.makeRotationY(32),
        trans.makeScaling(0.3, 1, 0.3),
    });
    try ex(cylinder.transform.t.equals(t.t));
}

test "parse cylinder (2)" {
    const txt =
        \\CYLINDER
        \\material flat (0.1,0.1,0.1)
        \\material diffuse 0.2
        \\material ambient 0.1
        \\material reflective 0.9
        \\material transparency 0.9
        \\material refractive 1.52
        \\closed true
        \\length 2
        \\transform scale (0.25,2,0.25)
        \\transform rotate x pi/2.0
        \\transform rotate y pi/4.0
        \\transform translate (-3,1.7,-0.5)
    ;

    const cylinder = try parseCylinder(txt);

    try ex(std.meta.activeTag(cylinder.material.color_map) == .flat);
    try ex(cylinder.material.color_map.flat.color.equals(Color.init(0.1, 0.1, 0.1)));
    try exEq(@as(bool, true), cylinder.closed);
    try exEq(@as(f64, 2), cylinder.length);
    try exEq(@as(f64, 0.2), cylinder.material.diffuse);
    try exEq(@as(f64, 0.1), cylinder.material.ambient);
    try exEq(@as(f64, 0.9), cylinder.material.reflective);
    try exEq(@as(f64, 0.9), cylinder.material.transparency);
    try exEq(@as(f64, 1.52), cylinder.material.refractive_index);

    // chain applies transforms in reverse order
    var t = (trans.Transform{}).chain(.{
        trans.makeTranslation(-3, 1.7, -0.5),
        trans.makeRotationY(PI / 4.0),
        trans.makeRotationX(PI / 2.0),
        trans.makeScaling(0.25, 2, 0.25),
    });
    try ex(cylinder.transform.t.equals(t.t));
}

test "parse cylinder (3)" {
    const txt =
        \\CYLINDER
        \\material 3dchecker (0.6,0.2,0.2) (0.2,0.2,0.6)
        \\length 3
        \\transform scale (0.8,0.4,0.8)
        \\transform rotate x pi/2
        \\transform rotate y -1.1
        \\transform translate (-3.5,0.8,-3.2)
    ;

    const cylinder = try parseCylinder(txt);

    try ex(std.meta.activeTag(cylinder.material.color_map) == .threedchecked);
    try ex(cylinder.material.color_map.threedchecked.a.equals(Color.init(0.6, 0.2, 0.2)));
    try ex(cylinder.material.color_map.threedchecked.b.equals(Color.init(0.2, 0.2, 0.6)));
    try exEq(@as(f64, 3), cylinder.length);

    // chain applies transforms in reverse order
    var t = (trans.Transform{}).chain(.{
        trans.makeTranslation(-3.5, 0.8, -3.2),
        trans.makeRotationY(-1.1),
        trans.makeRotationX(PI / 2.0),
        trans.makeScaling(0.8, 0.4, 0.8),
    });
    try ex(cylinder.transform.t.equals(t.t));
}

test "parse cone (1)" {
    const txt =
        \\CONE
        \\material 3dchecker (0.6,0.2,0.2) (0.2,0.2,0.6)
        \\max 0
        \\min -1
        \\closed false
        \\transform scale (0.3,1,0.3)
        \\transform rotate y 32
        \\transform translate (-2.2,3.9,-0.8)
    ;

    const cone = try parseCone(txt);

    try ex(std.meta.activeTag(cone.material.color_map) == .threedchecked);
    try ex(cone.material.color_map.threedchecked.a.equals(Color.init(0.6, 0.2, 0.2)));
    try ex(cone.material.color_map.threedchecked.b.equals(Color.init(0.2, 0.2, 0.6)));
    try exEq(@as(f64, 0), cone.max);
    try exEq(@as(f64, -1), cone.min);
    try exEq(@as(bool, false), cone.closed);

    // chain applies transforms in reverse order
    const t = (trans.Transform{}).chain(.{
        trans.makeTranslation(-2.2, 3.9, -0.8),
        trans.makeRotationY(32),
        trans.makeScaling(0.3, 1, 0.3),
    });
    try ex(cone.transform.t.equals(t.t));
}

test "parse cone (2)" {
    const txt =
        \\CONE
        \\material flat (0.2,0.2,0.2)
        \\max 1
        \\min 0
        \\closed true
        \\material diffuse 0.1
        \\material ambient 0.1
        \\material reflective 0.9
        \\transform translate (3,4,0)
    ;

    const cone = try parseCone(txt);

    try ex(std.meta.activeTag(cone.material.color_map) == .flat);
    try ex(cone.material.color_map.flat.color.equals(Color.init(0.2, 0.2, 0.2)));
    try exEq(@as(f64, 1), cone.max);
    try exEq(@as(f64, 0), cone.min);
    try exEq(@as(bool, true), cone.closed);
    try exEq(@as(f64, 0.1), cone.material.diffuse);
    try exEq(@as(f64, 0.1), cone.material.ambient);
    try exEq(@as(f64, 0.9), cone.material.reflective);
    try ex(cone.transform.t.equals(trans.makeTranslation(3, 4, 0).t));
}

test "parse tuples" {
    try exEq(
        @"3Tuple"{ .a = 1, .b = 2, .c = 3 },
        try parse3Tuple("(1,2,3)"),
    );
    try exEq(
        @"3Tuple"{ .a = 1, .b = -2, .c = 3 },
        try parse3Tuple("(1,-2,3)"),
    );
    try exEq(
        @"3Tuple"{ .a = 1, .b = 2, .c = 3.3 },
        try parse3Tuple("(1,2,3.3)"),
    );
    try exEq(
        @"3Tuple"{ .a = -1, .b = 2, .c = 3 },
        try parse3Tuple("(-1,2,3.0)"),
    );
    try exEq(
        @"3Tuple"{ .a = 1, .b = 2, .c = 3 },
        try parse3Tuple("(1,2,3)"),
    );
}

test "parse world (1)" {
    const txt =
        \\POINTLIGHT
        \\position (0,8,-1)
        \\intensity (1,1,1)
        \\
        \\CAMERA
        \\width 300
        \\height 300
        \\fov pi/2.5
        \\from (0,1,-4)
        \\to (0,1,0)
        \\up (0,1,0)
        \\
        \\CONE
        \\max 0
        \\min -1
        \\closed true
        \\tranform translate(0,1,0)
    ;

    const alctr = std.testing.allocator;
    const file = try parseWorldText(txt, alctr);
    const world = &file.world;
    const camera = &file.camera;
    defer world.deinit();

    try exEq(@as(usize, 1), world.lights_buf.items.len);
    try exEq(@as(usize, 1), world.cones_buf.items.len);
    try exEq(@as(i64, 300), camera.width);
    try exEq(@as(i64, 300), camera.height);
}

test "parse world (2)" {
    const txt =
        \\POINTLIGHT
        \\position (0,8,-1)
        \\intensity (1,1,1)
        \\
        \\CAMERA
        \\width 300
        \\height 300
        \\fov pi/2.5
        \\from (0,1,-4)
        \\to (0,1,0)
        \\up (0,1,0)
        \\
        \\SPHERE
        \\PLANE
        \\CUBE
        \\CYLINDER
        \\CONE
    ;

    const alctr = std.testing.allocator;
    const file = try parseWorldText(txt, alctr);
    const world = &file.world;
    const camera = &file.camera;
    defer world.deinit();

    try exEq(@as(usize, 1), world.lights_buf.items.len);
    try exEq(@as(usize, 1), world.spheres_buf.items.len);
    try exEq(@as(usize, 1), world.planes_buf.items.len);
    try exEq(@as(usize, 1), world.cubes_buf.items.len);
    try exEq(@as(usize, 1), world.cylinders_buf.items.len);
    try exEq(@as(usize, 1), world.cones_buf.items.len);
    try exEq(@as(i64, 300), camera.width);
    try exEq(@as(i64, 300), camera.height);
}

test "parse world (3)" {
    const txt =
        \\POINTLIGHT
        \\position (0,8,-1)
        \\intensity (1,1,1)
        \\
        \\CAMERA
        \\width 300
        \\height 300
        \\fov pi/2.5
        \\from (0,1,-4)
        \\to (0,1,0)
        \\up (0,1,0)
        \\
        \\CUBE
        \\CUBE
        \\CUBE
    ;

    const alctr = std.testing.allocator;
    const file = try parseWorldText(txt, alctr);
    const world = &file.world;
    const camera = &file.camera;
    defer world.deinit();

    try exEq(@as(usize, 3), world.cubes_buf.items.len);
    try exEq(@as(i64, 300), camera.width);
    try exEq(@as(i64, 300), camera.height);
}

const std = @import("std");
const trans = @import("transform.zig");
const mate = @import("material.zig");

const PI = std.math.pi;
const exEq = std.testing.expectEqual;
const ex = std.testing.expect;
const tprint = @import("u.zig").print;

const Point = @import("tuple.zig").Point;
const Vector = @import("tuple.zig").Vector;
const Color = @import("color.zig").Color;
const Transform = @import("transform.zig").Transform;
const Material = @import("material.zig").Material;

const World = @import("world.zig").World;
const Camera = @import("world.zig").Camera;
const PointLight = @import("light.zig").PointLight;
const Cube = @import("volume.zig").Cube;
const Plane = @import("volume.zig").Plane;
const Sphere = @import("volume.zig").Sphere;
const Cylinder = @import("volume.zig").Cylinder;
const Cone = @import("volume.zig").Cone;
