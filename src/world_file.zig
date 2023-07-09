//! create and populate a world from a text file

pub var error_message: ?[]const u8 = null;
pub var error_offending_line: ?[]const u8 = null;

pub const WorldParseError = error{
    IndescribableOffense,
};

const FileSection = struct {
    /// Reference to section of original text. Not owned by this object
    text: []const u8,
    header: []const u8,
    name: []const u8,
    parent_name: ?[]const u8,
    vptr: ?VolumePtr = null,
    data: Data,
    number: usize,

    const Data = union(enum) {
        camera: Camera,
        point_light: PointLight,
        cone: Cone,
        sphere: Sphere,
        plane: Plane,
        cube: Cube,
        cylinder: Cylinder,
        triangle: Triangle,
        smooth_triangle: SmoothTriangle,
    };

    pub fn inits(text: []const u8, alctr: std.mem.Allocator, number: usize) ![]This {
        const header = std.mem.sliceTo(text, '\n');

        const name = blk: {
            var split = std.mem.split(u8, header, " ");
            const shape = split.first();
            const name_text = split.rest();

            if (name_text.len > 0) break :blk name_text;
            break :blk shape;
        };

        const parent_name = blk: {
            var lines = std.mem.tokenize(u8, text, "\n");
            while (lines.next()) |line| {
                if (std.mem.startsWith(u8, line, "parent")) {
                    var split = std.mem.split(u8, line, " ");
                    _ = split.next();
                    break :blk split.rest();
                }
            }
            break :blk null;
        };

        const data = blk: {
            var result_data = std.ArrayList(Data).init(alctr);

            if (std.mem.startsWith(u8, header, "CAMERA")) {
                const datum = try parseCamera(text);
                try result_data.append(.{ .camera = datum });
            } else if (std.mem.startsWith(u8, header, "POINTLIGHT")) {
                const datum = try parsePointLight(text);
                try result_data.append(.{ .point_light = datum });
            } else if (std.mem.startsWith(u8, header, "CONE")) {
                const datum = try parseCone(text);
                try result_data.append(.{ .cone = datum });
            } else if (std.mem.startsWith(u8, header, "SPHERE")) {
                const datum = try parseSphere(text);
                try result_data.append(.{ .sphere = datum });
            } else if (std.mem.startsWith(u8, header, "PLANE")) {
                const datum = try parsePlane(text);
                try result_data.append(.{ .plane = datum });
            } else if (std.mem.startsWith(u8, header, "CUBE")) {
                const datum = try parseCube(text);
                try result_data.append(.{ .cube = datum });
            } else if (std.mem.startsWith(u8, header, "CYLINDER")) {
                const datum = try parseCylinder(text);
                try result_data.append(.{ .cylinder = datum });
            } else if (std.mem.startsWith(u8, header, "OBJ")) {
                // special case: OBJ sections should add a bunch of triange sections
                var tokens = std.mem.tokenize(u8, header, " ");
                _ = tokens.next(); // skip "OBJ"
                const filename = tokens.next().?;

                // get obj data
                var obj = try parseObjFile(filename, alctr);
                defer obj.deinit();

                // find materials or transforms
                var mat = Material.init();
                var tran = Transform{};
                var flat = false;
                {
                    var lines = std.mem.tokenize(u8, text, "\n");
                    while (lines.next()) |line| {
                        if (std.mem.startsWith(u8, line, "material")) {
                            try parseAndApplyMaterial(line, &mat);
                        } else if (std.mem.startsWith(u8, line, "transform")) {
                            try parseAndApplyTransform(line, &tran);
                        } else if (std.mem.startsWith(u8, line, "flat")) {
                            flat = try parseKeyValue(bool, line);
                        }
                    }
                }

                // XXX: some corrupt obj files can have bad geometry, e.g. a triangle
                // where two of its verts are the same. this messes up the renderer down
                // the line. we should filter those out here.

                for (obj.triangles) |*tri| {
                    tri.material = mat;
                    tri.transform = tran;
                    try result_data.append(.{ .triangle = tri.* });
                }

                for (obj.smooth_triangles) |*stri| {
                    if (flat) {
                        // if flat is specified, add all smooth triangles as flat
                        // triangles
                        var flat_tri = Triangle.init(stri.p1, stri.p2, stri.p3);
                        flat_tri.material = mat;
                        flat_tri.transform = tran;
                        try result_data.append(.{ .triangle = flat_tri });
                    } else {
                        stri.material = mat;
                        stri.transform = tran;
                        try result_data.append(.{ .smooth_triangle = stri.* });
                    }
                }
            } else {
                unreachable;
            }

            break :blk result_data.toOwnedSlice() catch unreachable;
        };
        defer alctr.free(data);

        var results = alctr.alloc(This, data.len) catch unreachable;

        for (data, 0..) |datum, i| {
            results[i] = This{
                .text = text,
                .header = header,
                .name = name,
                .parent_name = parent_name,
                .data = datum,
                .number = number,
            };
        }

        return results;
    }

    pub fn addObjectToWorld(self: *This, world: *World, camera: *?Camera) !void {
        switch (self.data) {
            .camera => |d| camera.* = d,
            .point_light => |d| {
                var new_light = world.pool.addLight(PointLight);
                new_light.ptr.* = d;
            },
            .cone => |c| {
                var new_cone = world.pool.addVolume(Cone);
                new_cone.ptr.* = c;
                self.vptr = new_cone.handle;
            },
            .sphere => |sph| {
                var new_sphere = world.pool.addVolume(Sphere);
                new_sphere.ptr.* = sph;
                self.vptr = new_sphere.handle;
            },
            .plane => |pl| {
                var new_plane = world.pool.addVolume(Plane);
                new_plane.ptr.* = pl;
                self.vptr = new_plane.handle;
            },
            .cube => |c| {
                var new_cube = world.pool.addVolume(Cube);
                new_cube.ptr.* = c;
                self.vptr = new_cube.handle;
            },
            .cylinder => |c| {
                var new_cylinder = world.pool.addVolume(Cylinder);
                new_cylinder.ptr.* = c;
                self.vptr = new_cylinder.handle;
            },
            .triangle => |t| {
                var new_triangle = world.pool.addVolume(Triangle);
                new_triangle.ptr.* = t;
                self.vptr = new_triangle.handle;
            },
            .smooth_triangle => |st| {
                var new_striangle = world.pool.addVolume(SmoothTriangle);
                new_striangle.ptr.* = st;
                self.vptr = new_striangle.handle;
            },
        }
    }

    const This = @This();
};

const FileContents = struct {
    world: World,
    camera: Camera,
    alctr: std.mem.Allocator,
    text: []const u8,
    sections: []FileSection,
};

const ParsedObj = struct {
    ignored_lines: usize,
    vertices: []Tuple,
    vertex_normals: []Tuple,
    triangles: []Triangle,
    smooth_triangles: []SmoothTriangle,
    alctr: std.mem.Allocator,

    pub fn init(txt: []const u8, alctr: std.mem.Allocator) !This {
        var lines = std.mem.tokenize(u8, txt, "\n\r");

        var ignored_lines: usize = 0;
        var verts = std.ArrayList(Tuple).init(alctr);
        var vert_normals = std.ArrayList(Tuple).init(alctr);
        var tris = std.ArrayList(Triangle).init(alctr);
        var smooth_tris = std.ArrayList(SmoothTriangle).init(alctr);

        while (lines.next()) |line| {
            // std.debug.print("/{s}/\n", .{line});
            if (std.mem.startsWith(u8, line, "v ")) {
                var nums = std.mem.tokenize(u8, line, " ");
                _ = nums.next(); // "v"

                const p1_txt = nums.next() orelse return unimplementedError();
                const p2_txt = nums.next() orelse return unimplementedError();
                const p3_txt = nums.next() orelse return unimplementedError();

                // std.debug.print("/{s}/{s}/{s}/\n", .{ p1_txt, p2_txt, p3_txt });
                const p1 = try std.fmt.parseFloat(f64, p1_txt);
                const p2 = try std.fmt.parseFloat(f64, p2_txt);
                const p3 = try std.fmt.parseFloat(f64, p3_txt);

                try verts.append(Point.init(p1, p2, p3));
            } else if (std.mem.startsWith(u8, line, "vn ")) {
                var nums = std.mem.tokenize(u8, line, " ");
                _ = nums.next(); // "vn"

                const n1_txt = nums.next() orelse return unimplementedError();
                const n2_txt = nums.next() orelse return unimplementedError();
                const n3_txt = nums.next() orelse return unimplementedError();

                // std.debug.print("/{s}/{s}/{s}/\n", .{ p1_txt, p2_txt, p3_txt });
                const n1 = try std.fmt.parseFloat(f64, n1_txt);
                const n2 = try std.fmt.parseFloat(f64, n2_txt);
                const n3 = try std.fmt.parseFloat(f64, n3_txt);

                try vert_normals.append(Vector.init(n1, n2, n3));
            } else if (std.mem.startsWith(u8, line, "f ")) {
                const data = line[2..];
                var nums = std.mem.tokenize(u8, data, " ");

                // f commands can have between 3 and inf vertices (where > 3 is describing
                // a more complex polygon). we fan triangulate them for the renderer.

                const p1_idx_txt = nums.next() orelse return unimplementedError();
                const p1_idc = try parseVertex(p1_idx_txt);

                var p2_idx_txt = nums.next() orelse return unimplementedError();
                var p2_idc = try parseVertex(p2_idx_txt);

                var p3_idx_txt_maybe = nums.next();

                while (p3_idx_txt_maybe) |p3_idx_txt| {
                    var p3_idc = try parseVertex(p3_idx_txt);

                    if (p1_idc.normal_idx != null and
                        p2_idc.normal_idx != null and
                        p3_idc.normal_idx != null)
                    {
                        // smooth triangle
                        try smooth_tris.append(SmoothTriangle.init(
                            verts.items[p1_idc.vertex_idx],
                            verts.items[p2_idc.vertex_idx],
                            verts.items[p3_idc.vertex_idx],
                            vert_normals.items[p1_idc.normal_idx.?],
                            vert_normals.items[p2_idc.normal_idx.?],
                            vert_normals.items[p3_idc.normal_idx.?],
                        ));
                    } else {
                        // flat triangle
                        try tris.append(Triangle.init(
                            verts.items[p1_idc.vertex_idx],
                            verts.items[p2_idc.vertex_idx],
                            verts.items[p3_idc.vertex_idx],
                        ));
                    }

                    // get next idx in the case of > 3 points
                    p2_idc = p3_idc;
                    p3_idx_txt_maybe = nums.next();
                }
            } else {
                // std.debug.print("Ignored line: \"{s}\"", .{line});
                ignored_lines += 1;
            }
        }

        return This{
            .ignored_lines = ignored_lines,
            .vertices = try verts.toOwnedSlice(),
            .vertex_normals = try vert_normals.toOwnedSlice(),
            .triangles = try tris.toOwnedSlice(),
            .smooth_triangles = try smooth_tris.toOwnedSlice(),
            .alctr = alctr,
        };
    }

    fn parseVertex(txt: []const u8) !struct {
        vertex_idx: usize,
        normal_idx: ?usize,
    } {
        // the data in f commands can come in a few forms:
        //   f 1 2 3
        //   f 1/2/3 4/5/6 7/8/9
        //   f 1//3 4//6 7//9
        //
        //   in the slash cases, the first number is the vertex index,
        //   the second is the texture vertex index (unused in this renderer),
        //   the third is the vertex normal index (TODO).

        var items = std.mem.split(u8, txt, "/");

        const item_1_txt = items.next() orelse return unimplementedError();
        _ = items.next(); // second item unused in this renderer
        const item_3_txt = items.next();

        // .obj files are 1-indexed, so we subtract 1 from the parsed idx
        const vert_idx = try std.fmt.parseInt(usize, item_1_txt, 10) - 1;
        const normal_idx = if (item_3_txt) |i3t| try std.fmt.parseInt(usize, i3t, 10) - 1 else null;

        return .{
            .vertex_idx = vert_idx,
            .normal_idx = normal_idx,
        };
    }

    pub fn deinit(self: *This) void {
        self.alctr.free(self.vertices);
        self.alctr.free(self.vertex_normals);
        self.alctr.free(self.triangles);
        self.alctr.free(self.smooth_triangles);
    }

    const This = @This();
};

pub fn parseWorldFile(filename: []const u8, alctr: std.mem.Allocator) !FileContents {
    var file = try std.fs.cwd().openFile(filename, .{});
    defer file.close();

    // a gigabyte
    var txt = try file.readToEndAlloc(alctr, 1024 * 1024 * 1024);

    // intentionally no defer alctr.free(txt), txt is owned by resulting FileContents

    return try parseWorldText(txt, alctr);
}

pub fn parseObjFile(filename: []const u8, alctr: std.mem.Allocator) !ParsedObj {
    var file = try std.fs.cwd().openFile(filename, .{});
    defer file.close();

    // a gigabyte
    var txt = try file.readToEndAlloc(alctr, 1024 * 1024 * 1024);
    defer alctr.free(txt);

    return try parseObj(txt, alctr);
}

pub fn parseWorldText(txt: []const u8, alctr: std.mem.Allocator) !FileContents {
    var world = World.init(alctr);
    var camera: ?Camera = null;

    var sections = try getTextSections(txt, alctr);
    defer sections.deinit();

    // populate world
    for (sections.items) |*section| {
        try section.addObjectToWorld(&world, &camera);
    }

    // link parents/children
    var tree = getVolumeTree(sections.items, alctr);
    defer tree.deinit();
    tree.applyTransforms(&world);

    // get triangle blocks
    var tri_blocks = getTriangleBlocks(sections.items, alctr);
    defer alctr.free(tri_blocks);

    // create AABBs
    for (tri_blocks) |tblock| {
        var bs = AABB.Bounds{
            .min_x = std.math.floatMax(f64),
            .min_y = std.math.floatMax(f64),
            .min_z = std.math.floatMax(f64),
            .max_x = std.math.floatMin(f64),
            .max_y = std.math.floatMin(f64),
            .max_z = std.math.floatMin(f64),
        };
        var transform: trans.Transform = undefined;
        switch (tblock.buffer) {
            .flat => {
                // find extents in object space
                for (world.pool.triangles_buf.items[tblock.start..tblock.end]) |t| {
                    bs.min_x = @min(t.p1.x(), t.p2.x(), t.p3.x(), bs.min_x);
                    bs.min_y = @min(t.p1.y(), t.p2.y(), t.p3.y(), bs.min_y);
                    bs.min_z = @min(t.p1.z(), t.p2.z(), t.p3.z(), bs.min_z);
                    bs.max_x = @max(t.p1.x(), t.p2.x(), t.p3.x(), bs.max_x);
                    bs.max_y = @max(t.p1.y(), t.p2.y(), t.p3.y(), bs.max_y);
                    bs.max_z = @max(t.p1.z(), t.p2.z(), t.p3.z(), bs.max_z);
                }
                // assuming all tforms are the same
                transform = world.pool.triangles_buf.items[tblock.start].transform;
            },
            .smooth => {
                // find extents in object space
                for (world.pool.smooth_triangles_buf.items[tblock.start..tblock.end]) |t| {
                    bs.min_x = @min(t.p1.x(), t.p2.x(), t.p3.x(), bs.min_x);
                    bs.min_y = @min(t.p1.y(), t.p2.y(), t.p3.y(), bs.min_y);
                    bs.min_z = @min(t.p1.z(), t.p2.z(), t.p3.z(), bs.min_z);
                    bs.max_x = @max(t.p1.x(), t.p2.x(), t.p3.x(), bs.max_x);
                    bs.max_y = @max(t.p1.y(), t.p2.y(), t.p3.y(), bs.max_y);
                    bs.max_z = @max(t.p1.z(), t.p2.z(), t.p3.z(), bs.max_z);
                }
                // assuming all tforms are the same
                transform = world.pool.smooth_triangles_buf.items[tblock.start].transform;
            },
        }

        var aabb = world.pool.addVolume(AABB);
        aabb.ptr.bounds = bs;
        aabb.ptr.transform = transform;
        aabb.ptr.first_idx = tblock.start;
        aabb.ptr.range = switch (tblock.buffer) {
            .flat => .{ .flat = world.pool.triangles_buf.items[tblock.start..tblock.end] },
            .smooth => .{ .smooth = world.pool.smooth_triangles_buf.items[tblock.start..tblock.end] },
        };
    }

    if (camera == null) return unimplementedError();

    return .{
        .world = world,
        .camera = camera.?,
        .alctr = alctr,
        .text = txt,
        .sections = sections.toOwnedSlice() catch unreachable,
    };
}

fn getTextSections(txt: []const u8, alctr: std.mem.Allocator) !std.ArrayList(FileSection) {
    var sections = std.ArrayList(FileSection).init(alctr);

    var remaining = txt;
    var chunk_idxs = findSectionBounds(remaining);

    var section_number: usize = 0;

    while (chunk_idxs.min) |min| {
        const chunk = remaining[min..chunk_idxs.max];

        const parsed_sections = try FileSection.inits(chunk, alctr, section_number);
        section_number += 1;
        defer alctr.free(parsed_sections);

        for (parsed_sections) |section| {
            sections.append(section) catch unreachable;
        }

        // std.debug.print("p {s}\n", .{sections.items[sections.items.len - 1].header()});
        if (chunk_idxs.max == remaining.len) break;

        remaining = remaining[chunk_idxs.max..];
        chunk_idxs = findSectionBounds(remaining);
    }

    return sections;
}

const VolumeTree = struct {
    nodes: []VolumeTreeNode,
    alctr: std.mem.Allocator,

    pub fn deinit(self: *VolumeTree) void {
        for (self.nodes) |*n| {
            n.children_idxs.deinit();
        }
        self.alctr.free(self.nodes);
    }

    fn printDepthFirst(self: VolumeTree, idx: usize, level: usize) void {
        {
            // indent
            var i: usize = 0;
            while (i < level) : (i += 1) {
                std.debug.print("  ", .{});
            }
        }

        if (self.nodes[idx].vptr) |vptr| {
            // print node
            std.debug.print("{}\n", .{vptr});

            // print children
            for (self.nodes[idx].children_idxs.items) |c_idx| {
                self.printDepthFirst(c_idx, level + 1);
            }
        }
    }

    pub fn print(self: VolumeTree) void {
        for (self.nodes, 0..) |n, i| {
            if (n.parent_idx == null) {
                self.printDepthFirst(i, 0);
            }
        }
    }

    fn applyTransformsDepthFirst(self: VolumeTree, world: *World, idx: usize, parent_tform: Transform) void {
        const vptr = self.nodes[idx].vptr.?;
        var vol_transform = world.pool.getProperty(vptr, "transform");

        vol_transform.* = parent_tform.mult(vol_transform.*);

        for (self.nodes[idx].children_idxs.items) |c_idx| {
            self.applyTransformsDepthFirst(world, c_idx, vol_transform.*);
        }
    }

    pub fn applyTransforms(self: VolumeTree, world: *World) void {
        for (self.nodes, 0..) |n, i| {
            if (n.parent_idx == null and n.vptr != null) {
                self.applyTransformsDepthFirst(world, i, Transform{});
            }
        }
    }
};

const VolumeTreeNode = struct {
    vptr: ?VolumePtr,
    parent_idx: ?usize,
    children_idxs: std.ArrayList(usize),
};

fn getVolumeTree(sections: []FileSection, alctr: std.mem.Allocator) VolumeTree {
    var tree_nodes = alctr.alloc(VolumeTreeNode, sections.len) catch unreachable;

    // initialize flat tree
    for (sections, 0..) |s, i| {
        const pidx: ?usize = blk: {
            if (s.parent_name) |pname| {
                for (sections, 0..) |ps, j| {
                    if (std.mem.eql(u8, ps.name, pname))
                        break :blk j;
                }
                unreachable;
            }
            break :blk null;
        };

        tree_nodes[i] = VolumeTreeNode{
            .vptr = s.vptr,
            .parent_idx = pidx,
            .children_idxs = std.ArrayList(usize).init(alctr),
        };
    }

    // populate children lists
    for (tree_nodes, 0..) |tn, i| {
        if (tn.parent_idx) |pi| {
            tree_nodes[pi].children_idxs.append(i) catch unreachable;
        }
    }

    return VolumeTree{
        .nodes = tree_nodes,
        .alctr = alctr,
    };
}

const TriangleBlock = struct {
    const Buffer = enum {
        flat,
        smooth,
    };

    start: usize,
    end: usize,
    buffer: Buffer,
};

// scan the file section for blocks of triangles or smooth triangles and produce a
// list of AABBs enclosing each group of obj file triangles.
fn getTriangleBlocks(sections: []FileSection, alctr: std.mem.Allocator) []TriangleBlock {
    var result = std.ArrayList(TriangleBlock).init(alctr);

    var start_it: usize = 0;
    var end_it: usize = 0;

    var next_flat_idx: usize = 0;
    var next_smooth_idx: usize = 0;

    while (start_it < sections.len) : (start_it = end_it) {
        // find block bounds. scan forward until off the list or different section number
        while (end_it < sections.len and sections[start_it].number == sections[end_it].number)
            end_it += 1;

        var block = sections[start_it..end_it];

        // if it's a triangle, figure out which buffer it will be in and make a block
        switch (sections[start_it].data) {
            .triangle => {
                result.append(TriangleBlock{
                    .start = next_flat_idx,
                    .end = next_flat_idx + (block.len),
                    .buffer = .flat,
                }) catch unreachable;
                next_flat_idx += (block.len);
            },
            .smooth_triangle => {
                result.append(TriangleBlock{
                    .start = next_smooth_idx,
                    .end = next_smooth_idx + (block.len),
                    .buffer = .smooth,
                }) catch unreachable;
                next_smooth_idx += (block.len);
            },
            else => {}, // skip non-triangle
        }
    }

    return result.toOwnedSlice() catch unreachable;
}

fn findSectionBounds(txt: []const u8) struct { min: ?usize, max: usize } {
    const keys = [_][]const u8{
        "CAMERA", "POINTLIGHT",
        "SPHERE", "PLANE",
        "CUBE",   "CYLINDER",
        "CONE",   "OBJ",
    };

    var keys_idxs = [_]?usize{null} ** keys.len;

    const min = blk: {
        for (keys, 0..) |key, i| {
            keys_idxs[i] = std.mem.indexOf(u8, txt, key);
        }

        const lt = struct {
            fn lt(_: void, a: ?usize, b: ?usize) bool {
                if (a == null) return false;
                return a.? < b orelse std.math.maxInt(usize);
            }
        }.lt;

        std.sort.insertion(?usize, &keys_idxs, {}, lt);

        break :blk keys_idxs[0].?;
    };

    const max = blk: {
        for (keys, 0..) |key, i| {
            keys_idxs[i] = std.mem.indexOf(u8, txt[min + 1 ..], key);
            if (keys_idxs[i] != null) keys_idxs[i].? += 1;
        }

        const lt = struct {
            fn lt(_: void, a: ?usize, b: ?usize) bool {
                if (a == null) return false;
                return a.? < b orelse std.math.maxInt(usize);
            }
        }.lt;

        std.sort.insertion(?usize, &keys_idxs, {}, lt);

        break :blk keys_idxs[0] orelse txt.len;
    };

    return .{ .min = min, .max = max };
}

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
    const fwidth = @as(f64, @floatFromInt(width.?)) * scale.?;
    const fheight = @as(f64, @floatFromInt(height.?)) * scale.?;

    var camera = Camera.init(@as(i64, @intFromFloat(fwidth)), @as(i64, @intFromFloat(fheight)), fov.?);
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

fn parseObj(txt: []const u8, alctr: std.mem.Allocator) !ParsedObj {
    return try ParsedObj.init(txt, alctr);
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
    defer alctr.free(file.sections);
    const world = &file.world;
    const camera = &file.camera;
    defer world.deinit();

    try exEq(@as(usize, 1), world.pool.lights_buf.items.len);
    try exEq(@as(usize, 1), world.pool.cones_buf.items.len);
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
    defer alctr.free(file.sections);
    const world = &file.world;
    const camera = &file.camera;
    defer world.deinit();

    try exEq(@as(usize, 1), world.pool.lights_buf.items.len);
    try exEq(@as(usize, 1), world.pool.spheres_buf.items.len);
    try exEq(@as(usize, 1), world.pool.planes_buf.items.len);
    try exEq(@as(usize, 1), world.pool.cubes_buf.items.len);
    try exEq(@as(usize, 1), world.pool.cylinders_buf.items.len);
    try exEq(@as(usize, 1), world.pool.cones_buf.items.len);
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
    defer alctr.free(file.sections);
    const world = &file.world;
    const camera = &file.camera;
    defer world.deinit();

    try exEq(@as(usize, 3), world.pool.cubes_buf.items.len);
    try exEq(@as(i64, 300), camera.width);
    try exEq(@as(i64, 300), camera.height);
}

test "parse world (4)" {
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
        \\OBJ obj/test.obj
    ;

    const alctr = std.testing.allocator;
    const file = try parseWorldText(txt, alctr);
    defer alctr.free(file.sections);
    const world = &file.world;
    defer world.deinit();

    // test.obj has 3 triangles in it
    try exEq(@as(usize, 3), world.pool.triangles_buf.items.len);
    try ex(world.pool.triangles_buf.items[0].equals(Triangle.init(
        Point.init(-1, 1, 0),
        Point.init(-1, 0, 0),
        Point.init(1, 0, 0),
    )));
    try ex(world.pool.triangles_buf.items[1].equals(Triangle.init(
        Point.init(-1, 1, 0),
        Point.init(1, 0, 0),
        Point.init(1, 1, 0),
    )));
    try ex(world.pool.triangles_buf.items[2].equals(Triangle.init(
        Point.init(-1, 1, 0),
        Point.init(1, 1, 0),
        Point.init(0, 2, 0),
    )));
}

test "parse parents: parser parses names and parent names" {
    const txt =
        \\CAMERA
        \\width 300
        \\height 300
        \\fov pi/2.5
        \\from (0,1,-4)
        \\to (0,1,0)
        \\up (0,1,0)
        \\
        \\CONE name 1
        \\max 1
        \\min 0
        \\closed false
        \\
        \\CONE
        \\max 1
        \\min 0
        \\closed false
        \\
        \\SPHERE midlevel
        \\parent name 1
        \\
        \\SPHERE
        \\parent midlevel
    ;

    const alctr = std.testing.allocator;
    const file = try parseWorldText(txt, alctr);
    defer alctr.free(file.sections);
    const world = &file.world;
    defer world.deinit();

    try exEq(@as(usize, 5), file.sections.len);

    try exEqStr("CAMERA", std.mem.sliceTo(file.sections[0].name, 0));
    try exEq(@as(?[]const u8, null), file.sections[0].parent_name);

    try exEqStr("name 1", std.mem.sliceTo(file.sections[1].name, 0));
    try exEq(@as(?[]const u8, null), file.sections[1].parent_name);

    try exEqStr("CONE", std.mem.sliceTo(file.sections[2].name, 0));
    try exEq(@as(?[]const u8, null), file.sections[2].parent_name);

    try exEqStr("midlevel", std.mem.sliceTo(file.sections[3].name, 0));
    try ex(file.sections[3].parent_name != null);
    try exEqStr("name 1", std.mem.sliceTo(file.sections[3].parent_name.?, 0));

    try exEqStr("SPHERE", std.mem.sliceTo(file.sections[4].name, 0));
    try ex(file.sections[4].parent_name != null);
    try exEqStr("midlevel", std.mem.sliceTo(file.sections[4].parent_name.?, 0));
}

test "parse parents: parser parses vptrs" {
    const txt =
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
    defer alctr.free(file.sections);
    const world = &file.world;
    defer world.deinit();

    try exEq(@as(usize, 6), file.sections.len);

    try exEq(VolumePtr{ .sphere_idx = 0 }, file.sections[1].vptr.?);
    try exEq(VolumePtr{ .plane_idx = 0 }, file.sections[2].vptr.?);
    try exEq(VolumePtr{ .cube_idx = 0 }, file.sections[3].vptr.?);
    try exEq(VolumePtr{ .cylinder_idx = 0 }, file.sections[4].vptr.?);
    try exEq(VolumePtr{ .cone_idx = 0 }, file.sections[5].vptr.?);
}

test "parse parents: parser applies transforms: translation" {
    const txt =
        \\CAMERA
        \\width 300
        \\height 300
        \\fov pi/2.5
        \\from (0,1,-4)
        \\to (0,1,0)
        \\up (0,1,0)
        \\
        \\CONE gem bottom
        \\max 1
        \\min 0
        \\closed false
        \\transform translate (1,0,0)
        \\
        \\CONE gem top
        \\max 1
        \\min 0
        \\closed false
        \\transform translate (0,2,0)
        \\parent gem bottom
    ;

    const alctr = std.testing.allocator;
    const file = try parseWorldText(txt, alctr);
    defer alctr.free(file.sections);
    const world = &file.world;
    defer world.deinit();

    try exEq(@as(usize, 2), world.pool.cones_buf.items.len);

    // assuming that they are added to the world in order...

    errdefer tprint(world.pool.cones_buf.items[1].transform.t);
    errdefer tprint(world.pool.cones_buf.items[0].transform.t);

    try ex(world.pool.cones_buf.items[0].transform.t.equals(trans.makeTranslation(1, 0, 0).t));
    // test that the top was moved along with the bottom
    try ex(world.pool.cones_buf.items[1].transform.t.equals(trans.makeTranslation(1, 2, 0).t));
}

test "parse parents: parser applies transforms: rotation" {
    const txt =
        \\CAMERA
        \\scale 1
        \\width 12
        \\height 12
        \\fov pi/2.5
        \\from (-0.2,3.25,-10)
        \\to (-0.2,1.75,0)
        \\up (0,1,0)
        \\
        \\SPHERE origin
        \\transform scale (0.5,0.5,0.5)
        \\transform rotate z pi/2
        \\
        \\SPHERE right
        \\transform translate (3,0,0)
        \\parent origin
        \\
        \\SPHERE
        \\transform translate (0,3,0)
        \\parent origin
        \\
        \\SPHERE
        \\transform translate (0,3,0)
        \\parent right
    ;

    const alctr = std.testing.allocator;
    const file = try parseWorldText(txt, alctr);
    defer alctr.free(file.sections);
    const world = &file.world;
    defer world.deinit();

    try exEq(@as(usize, 4), world.pool.spheres_buf.items.len);

    // assuming that they are added to the world in order...

    errdefer {
        tprint(world.pool.spheres_buf.items[0].transform.t);
        tprint(world.pool.spheres_buf.items[1].transform.t);
        tprint(world.pool.spheres_buf.items[2].transform.t);
        tprint(world.pool.spheres_buf.items[3].transform.t);
    }

    const identity = trans.Transform{};

    {
        const test_t_1 =
            identity.chain(.{
            trans.makeRotationZ(PI / 2.0),
            trans.makeScaling(0.5, 0.5, 0.5),
        }).t;
        try ex(world.pool.spheres_buf.items[0].transform.t.equals(test_t_1));
    }

    {
        const test_t_2 = identity.chain(.{
            trans.makeTranslation(0, 1.5, 0),
            trans.makeScaling(0.5, 0.5, 0.5),
            trans.makeRotationZ(PI / 2.0),
        }).t;
        try ex(world.pool.spheres_buf.items[1].transform.t.equals(test_t_2));
    }

    {
        const test_t_3 = identity.chain(.{
            trans.makeTranslation(-1.5, 0, 0),
            trans.makeRotationZ(PI / 2.0),
            trans.makeScaling(0.5, 0.5, 0.5),
        }).t;
        try ex(world.pool.spheres_buf.items[2].transform.t.equals(test_t_3));
    }

    {
        const test_t_4 = identity.chain(.{
            trans.makeTranslation(-1.5, 1.5, 0),
            trans.makeRotationZ(PI / 2.0),
            trans.makeScaling(0.5, 0.5, 0.5),
        }).t;
        try ex(world.pool.spheres_buf.items[3].transform.t.equals(test_t_4));
    }
}

test "Obj: Ignore unrecognized lines" {
    const txt =
        \\This obj file
        \\is not formatted
        \\in the right format.
        \\The parser should
        \\ignore all of it.
    ;

    const alctr = std.testing.allocator;
    var obj = try parseObj(txt, alctr);
    defer obj.deinit();

    try exEq(@as(usize, 5), obj.ignored_lines);
}

test "Obj: parse vertex records" {
    const txt =
        \\v -1 1 0
        \\v -1.0000 0.5000 0.0000
        \\v 1 0 0
        \\v 1 1 0
    ;

    const alctr = std.testing.allocator;
    var obj = try parseObj(txt, alctr);
    defer obj.deinit();

    try exEq(@as(usize, 4), obj.vertices.len);
    try ex(obj.vertices[0].equals(Point.init(-1, 1, 0)));
    try ex(obj.vertices[1].equals(Point.init(-1, 0.5, 0)));
    try ex(obj.vertices[2].equals(Point.init(1, 0, 0)));
    try ex(obj.vertices[3].equals(Point.init(1, 1, 0)));
}

test "Obj: parse triangles" {
    const txt =
        \\v -1 1 0
        \\v -1 0 0
        \\v 1 0 0
        \\v 1 1 0
        \\
        \\f 1 2 3
        \\f 1 3 4
    ;

    const alctr = std.testing.allocator;
    var obj = try parseObj(txt, alctr);
    defer obj.deinit();

    try exEq(@as(usize, 4), obj.vertices.len);
    try ex(obj.vertices[0].equals(Point.init(-1, 1, 0)));
    try ex(obj.vertices[1].equals(Point.init(-1, 0, 0)));
    try ex(obj.vertices[2].equals(Point.init(1, 0, 0)));
    try ex(obj.vertices[3].equals(Point.init(1, 1, 0)));

    try exEq(@as(usize, 2), obj.triangles.len);
    try ex(obj.triangles[0].equals(Triangle.init(
        obj.vertices[0],
        obj.vertices[1],
        obj.vertices[2],
    )));
    try ex(obj.triangles[1].equals(Triangle.init(
        obj.vertices[0],
        obj.vertices[2],
        obj.vertices[3],
    )));
}

test "Obj: parse polygons" {
    const txt =
        \\v -1 1 0
        \\v -1 0 0
        \\v 1 0 0
        \\v 1 1 0
        \\v 0 2 0
        \\
        \\f 1 2 3 4 5
    ;

    const alctr = std.testing.allocator;
    var obj = try parseObj(txt, alctr);
    defer obj.deinit();

    try exEq(@as(usize, 5), obj.vertices.len);
    try ex(obj.vertices[0].equals(Point.init(-1, 1, 0)));
    try ex(obj.vertices[1].equals(Point.init(-1, 0, 0)));
    try ex(obj.vertices[2].equals(Point.init(1, 0, 0)));
    try ex(obj.vertices[3].equals(Point.init(1, 1, 0)));
    try ex(obj.vertices[4].equals(Point.init(0, 2, 0)));

    // parser will triangulate the polygon
    try exEq(@as(usize, 3), obj.triangles.len);
    try ex(obj.triangles[0].equals(Triangle.init(
        obj.vertices[0],
        obj.vertices[1],
        obj.vertices[2],
    )));
    try ex(obj.triangles[1].equals(Triangle.init(
        obj.vertices[0],
        obj.vertices[2],
        obj.vertices[3],
    )));
    try ex(obj.triangles[2].equals(Triangle.init(
        obj.vertices[0],
        obj.vertices[3],
        obj.vertices[4],
    )));
}

test "Obj: parse polygons: complex faces" {
    const txt =
        \\v -1 1 0
        \\v -1 0 0
        \\v 1 0 0
        \\v 1 1 0
        \\v 0 2 0
        \\
        \\vn 0 2 0
        \\vn 0 2 0
        \\vn 0 2 0
        \\vn 0 2 0
        \\
        \\f 1//1 2//2 3/3/3 4/4/4 5
    ;

    const alctr = std.testing.allocator;
    var obj = try parseObj(txt, alctr);
    defer obj.deinit();

    try exEq(@as(usize, 5), obj.vertices.len);
    try ex(obj.vertices[0].equals(Point.init(-1, 1, 0)));
    try ex(obj.vertices[1].equals(Point.init(-1, 0, 0)));
    try ex(obj.vertices[2].equals(Point.init(1, 0, 0)));
    try ex(obj.vertices[3].equals(Point.init(1, 1, 0)));
    try ex(obj.vertices[4].equals(Point.init(0, 2, 0)));

    // parser will triangulate the polygon
    try exEq(@as(usize, 2), obj.smooth_triangles.len);
    try ex(obj.smooth_triangles[0].equals(SmoothTriangle.init(
        obj.vertices[0],
        obj.vertices[1],
        obj.vertices[2],
        obj.vertex_normals[0],
        obj.vertex_normals[1],
        obj.vertex_normals[2],
    )));
    try ex(obj.smooth_triangles[1].equals(SmoothTriangle.init(
        obj.vertices[0],
        obj.vertices[2],
        obj.vertices[3],
        obj.vertex_normals[0],
        obj.vertex_normals[2],
        obj.vertex_normals[3],
    )));

    try exEq(@as(usize, 1), obj.triangles.len);
    try ex(obj.triangles[0].equals(Triangle.init(
        obj.vertices[0],
        obj.vertices[3],
        obj.vertices[4],
    )));
}

test "Obj: Can parse suzanne without exploding" {
    var file = try std.fs.cwd().openFile("obj/suzanne.obj", .{});
    defer file.close();

    const alctr = std.testing.allocator;

    var txt = try file.readToEndAlloc(alctr, 1024 * 1024);
    defer alctr.free(txt);

    var obj = try parseObj(txt, alctr);
    defer obj.deinit();
}

test "Obj: Can parse teapot-low without exploding" {
    var file = try std.fs.cwd().openFile("obj/teapot-low.obj", .{});
    defer file.close();

    const alctr = std.testing.allocator;

    var txt = try file.readToEndAlloc(alctr, 1024 * 1024);
    defer alctr.free(txt);

    var obj = try parseObj(txt, alctr);
    defer obj.deinit();
}

test "Obj: apply material" {
    const txt =
        \\CAMERA
        \\scale 1
        \\width 12
        \\height 12
        \\fov pi/2.5
        \\from (-0.2,3.25,-10)
        \\to (-0.2,1.75,0)
        \\up (0,1,0)
        \\
        \\OBJ obj/test.obj
        \\material flat (1,0,0)
    ;

    const alctr = std.testing.allocator;
    const file = try parseWorldText(txt, alctr);
    defer alctr.free(file.sections);
    const world = &file.world;
    defer world.deinit();

    try exEq(@as(usize, 3), world.pool.triangles_buf.items.len);
    try ex(world.pool.triangles_buf.items[0].material.color_map.flat.color.equals(Color.init(1, 0, 0)));
}

test "Obj: apply transform" {
    const txt =
        \\CAMERA
        \\scale 1
        \\width 12
        \\height 12
        \\fov pi/2.5
        \\from (-0.2,3.25,-10)
        \\to (-0.2,1.75,0)
        \\up (0,1,0)
        \\
        \\OBJ obj/test.obj
        \\transform translate (0,1,0)
    ;

    const alctr = std.testing.allocator;
    const file = try parseWorldText(txt, alctr);
    defer alctr.free(file.sections);
    const world = &file.world;
    defer world.deinit();

    try exEq(@as(usize, 3), world.pool.triangles_buf.items.len);
    try ex(world.pool.triangles_buf.items[0].transform.t.equals(trans.makeTranslation(0, 1, 0).t));
    try ex(world.pool.triangles_buf.items[1].transform.t.equals(trans.makeTranslation(0, 1, 0).t));
    try ex(world.pool.triangles_buf.items[2].transform.t.equals(trans.makeTranslation(0, 1, 0).t));
}

test "obj: parent" {
    const txt =
        \\CAMERA
        \\width 300
        \\height 300
        \\fov pi/2.5
        \\from (0,1,-4)
        \\to (0,1,0)
        \\up (0,1,0)
        \\
        \\SPHERE middle
        \\transform translate (0,1,0)
        \\
        \\OBJ obj/test.obj
        \\parent middle
    ;

    const alctr = std.testing.allocator;
    const file = try parseWorldText(txt, alctr);
    defer alctr.free(file.sections);
    const world = &file.world;
    defer world.deinit();

    try exEq(@as(usize, 3), world.pool.triangles_buf.items.len);
    try ex(world.pool.triangles_buf.items[0].transform.t.equals(trans.makeTranslation(0, 1, 0).t));
    try ex(world.pool.triangles_buf.items[1].transform.t.equals(trans.makeTranslation(0, 1, 0).t));
    try ex(world.pool.triangles_buf.items[2].transform.t.equals(trans.makeTranslation(0, 1, 0).t));
}

test "obj: parse vertex normals" {
    const txt =
        \\vn 0 0 1
        \\vn 0.707 0 -0.707
        \\vn 1 2 3
    ;

    const alctr = std.testing.allocator;
    var obj = try parseObj(txt, alctr);
    defer obj.deinit();

    try exEq(@as(usize, 3), obj.vertex_normals.len);
    try ex(obj.vertex_normals[0].equals(Vector.init(0, 0, 1)));
    try ex(obj.vertex_normals[1].equals(Vector.init(0.707, 0, -0.707)));
    try ex(obj.vertex_normals[2].equals(Vector.init(1, 2, 3)));
}

test "obj: parse vertex normals 2" {
    const txt =
        \\v 0 1 0
        \\v -1 0 0
        \\v 1 0 0
        \\
        \\vn -1 0 0
        \\vn 1 0 0
        \\vn 0 1 0
        \\
        \\f 1//3 2//1 3//2
        \\f 1/0/3 2/102/1 3/14/2
    ;

    const alctr = std.testing.allocator;
    var obj = try parseObj(txt, alctr);
    defer obj.deinit();

    try exEq(@as(usize, 3), obj.vertices.len);
    try ex(obj.vertices[0].equals(Point.init(0, 1, 0)));
    try ex(obj.vertices[1].equals(Point.init(-1, 0, 0)));
    try ex(obj.vertices[2].equals(Point.init(1, 0, 0)));

    try exEq(@as(usize, 3), obj.vertex_normals.len);
    try ex(obj.vertex_normals[0].equals(Vector.init(-1, 0, 0)));
    try ex(obj.vertex_normals[1].equals(Vector.init(1, 0, 0)));
    try ex(obj.vertex_normals[2].equals(Vector.init(0, 1, 0)));

    try exEq(@as(usize, 2), obj.smooth_triangles.len);
    try ex(obj.smooth_triangles[0].equals(SmoothTriangle.init(
        obj.vertices[0],
        obj.vertices[1],
        obj.vertices[2],
        obj.vertex_normals[2],
        obj.vertex_normals[0],
        obj.vertex_normals[1],
    )));
    try ex(obj.smooth_triangles[1].equals(SmoothTriangle.init(
        obj.vertices[0],
        obj.vertices[1],
        obj.vertices[2],
        obj.vertex_normals[2],
        obj.vertex_normals[0],
        obj.vertex_normals[1],
    )));
}

const std = @import("std");
const trans = @import("transform.zig");
const mate = @import("material.zig");

const PI = std.math.pi;
const exEqStr = std.testing.expectEqualStrings;
const exEq = std.testing.expectEqual;
const ex = std.testing.expect;
const tprint = @import("u.zig").print;

const Tuple = @import("tuple.zig").Tuple;
const Point = @import("tuple.zig").Point;
const Vector = @import("tuple.zig").Vector;
const Color = @import("color.zig").Color;
const Transform = @import("transform.zig").Transform;
const Material = @import("material.zig").Material;

const World = @import("world.zig").World;
const Camera = @import("world.zig").Camera;
const PointLight = @import("light.zig").PointLight;
const Cube = @import("volume.zig").Cube;
const AABB = @import("volume.zig").AABB;
const Plane = @import("volume.zig").Plane;
const Sphere = @import("volume.zig").Sphere;
const Cylinder = @import("volume.zig").Cylinder;
const Cone = @import("volume.zig").Cone;
const Triangle = @import("volume.zig").Triangle;
const SmoothTriangle = @import("volume.zig").SmoothTriangle;

const VolumePtr = @import("volume.zig").VolumePool.VolumePtr;
