//! render threads

const std = @import("std");
const intx = @import("intersect.zig");
const light = @import("light.zig");

const Tuple = @import("tuple.zig").Tuple;
const Point = @import("tuple.zig").Point;
const Vector = @import("tuple.zig").Vector;
const Ray = @import("ray.zig").Ray;
const Color = @import("color.zig").Color;

const World = @import("world.zig").World;
const Sphere = @import("sphere.zig").Sphere;
const VolPtr = @import("world.zig").VolumePtr;

const Intersections = @import("intersect.zig").Intersections;

const Qanvas = @import("qanvas.zig").Qanvas;

const Chunk = struct {
    start_x: i64,
    end_x: i64,
    start_y: i64,
    end_y: i64,
    done: bool,
};

fn getChunks(w: i64, qan: Qanvas, alctr: std.mem.Allocator) []Chunk {
    var chunks = std.ArrayList(Chunk).init(alctr);

    const y_w = w;

    const chunks_x = @divTrunc(@intCast(i64, qan.width) -| 1, w) + 1;
    const chunks_y = @divTrunc(@intCast(i64, qan.height) -| 1, y_w) + 1;

    var y: i64 = 0;
    while (y < chunks_y) : (y += 1) {
        var x: i64 = 0;
        while (x < chunks_x) : (x += 1) {
            chunks.append(.{
                .start_x = w * x,
                .end_x = @min(@intCast(i64, qan.width), w * (x + 1)),
                .start_y = y_w * y,
                .end_y = @min(@intCast(i64, qan.height), y_w * (y + 1)),
                .done = false,
            }) catch unreachable;
        }
    }

    return chunks.toOwnedSlice() catch unreachable;
}

pub fn startRenderEngine(world: World, window_center: Tuple, heading: Tuple, qan: *Qanvas, alctr: std.mem.Allocator) void {
    std.debug.print("Starting render.\n", .{});
    var prog_ctx = std.Progress{};
    var prog = prog_ctx.start("Pixels", qan.width * qan.height);

    // TODO base this size on cache size?
    var chunks = getChunks(32, qan.*, alctr);
    defer alctr.free(chunks);

    var prng = std.rand.DefaultPrng.init(0);
    var random = prng.random();
    random.shuffle(Chunk, chunks);

    const cpus = std.Thread.getCpuCount() catch unreachable;
    const num_threads = @divTrunc(cpus * 3, 4) + 1;
    std.debug.print("Using {} threads.\n", .{@min(num_threads, 16)});

    var threads_idle_buf = [_]bool{true} ** 16;
    var threads_idle = threads_idle_buf[0..num_threads];

    var timer = std.time.Timer.start() catch unreachable;

    while (true) {
        // find chunk that needs rendering
        const mchunk: ?*Chunk = blk: {
            for (chunks) |*chunk| {
                if (!chunk.done) break :blk chunk;
            }
            break :blk null;
        };

        // start thread on chunk
        if (mchunk) |chunk| {
            const mthread_i: ?usize = blk: {
                for (threads_idle) |ti, i| {
                    if (ti == true) {
                        break :blk i;
                    }
                }
                break :blk null;
            };

            if (mthread_i) |ti| {
                // start thread
                // zig fmt: off
                threads_idle[ti] = false;
                chunk.done = true;
                var thr = std.Thread.spawn(.{}, render,
                    .{world, window_center, heading, qan, prog, alctr,
                    chunk.start_x, chunk.end_x, chunk.start_y, chunk.end_y,
                    &threads_idle[ti]}
                ) catch unreachable;
                thr.detach();
                // zig fmt: on
            } else {
                // wait and try again
                std.time.sleep(1_000);
            }
        } else {
            break;
        }
    }

    // wait for threads to finish
    while (true) {
        const workth_i: ?usize = blk: {
            for (threads_idle) |ti, i| {
                if (ti == false) {
                    break :blk i;
                }
            }
            break :blk null;
        };

        if (workth_i != null) {
            // wait and try again
            std.time.sleep(1_000);
        } else {
            break;
        }
    }

    prog.end();
    std.debug.print("Done in {d:.3} seconds.\n", .{@intToFloat(f64, timer.read()) / 1_000_000_000});
}

// zig fmt: off
fn render(
    world: World,
    window_center: Tuple,
    heading: Tuple,
    qan: *Qanvas,
    prog: *std.Progress.Node,
    alctr: std.mem.Allocator,
    start_x: i64, end_x: i64,
    start_y: i64, end_y: i64,
    done_flag: *bool
) void {
// zig fmt: on
    defer done_flag.* = true;
    const cam = window_center.plus(heading.normalized().scaled(-2));

    {
        var pix_y: i64 = start_y;
        while (pix_y < end_y) : (pix_y += 1) {
            var pix_x: i64 = start_x;
            while (pix_x < end_x) : (pix_x += 1) {
                qan.write(Color.init(1, 0, 1), pix_x, pix_y);
            }
        }
    }

    var intersections = Intersections.init(alctr, .{});
    defer intersections.deinit();

    var pix_y: i64 = start_y;
    while (pix_y < end_y) : (pix_y += 1) {
        var pix_x: i64 = start_x;
        while (pix_x < end_x) : (pix_x += 1) {
            const pixels_per_unit = 256;
            // get point on window
            const pix_x_w = @intToFloat(f64, pix_x - @divTrunc(@intCast(i64, qan.width), 2)) / pixels_per_unit;
            const pix_y_w = @intToFloat(f64, pix_y - @divTrunc(@intCast(i64, qan.height), 2)) / pixels_per_unit;
            // note the -pix_y_w. pixel space y increases downwards
            // but object space y increases upwards
            const aim_point = window_center.plus(Vector.init(pix_x_w, -pix_y_w, 0));
            const aim_ray = Ray.init(cam, aim_point.minus(cam).normalized());

            for (world.spheres_buf.items) |*sphptr, i| {
                const sph = VolPtr{ .sphere_idx = i };
                intersections.intersect(sph, sphptr.*, aim_ray);
            }
            defer intersections.clear();

            if (intersections.hit()) |x| {
                const p = aim_ray.position(x.t);
                const volume = world.get(Sphere, x.vptr);
                const n = volume.normalAt(p);
                const e = aim_ray.direction.scaled(-1);
                const l = world.get(light.PointLight, VolPtr{ .light_idx = 0 });

                const color = light.lighting(volume.material, l.*, p, e, n);
                qan.write(color, pix_x, pix_y);
            } else {
                qan.write(qan.backgroundColor(pix_x, pix_y), pix_x, pix_y);
            }

            prog.completeOne();
        }
    }
}
