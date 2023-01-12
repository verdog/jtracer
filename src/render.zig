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
const Camera = @import("world.zig").Camera;
const Sphere = @import("volume.zig").Sphere;
const VolPtr = @import("world.zig").VolumePtr;

const Qanvas = @import("qanvas.zig").Qanvas;

const Intersections = @import("intersect.zig").Intersections;

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

pub fn startRenderEngine(world: World, cam: Camera, qan: *Qanvas, alctr: std.mem.Allocator) void {
    std.debug.print("Starting render.\n", .{});
    var prog_ctx = std.Progress{};
    var prog = prog_ctx.start("Pixels", qan.width * qan.height);

    // TODO base this size on cache size?
    var chunks = getChunks(128, qan.*, alctr);
    defer alctr.free(chunks);

    var prng = std.rand.DefaultPrng.init(0);
    var random = prng.random();
    random.shuffle(Chunk, chunks);

    const cpus = std.Thread.getCpuCount() catch unreachable;
    const num_threads = @divTrunc(cpus * 3, 4) + 1;
    std.debug.print("Using {} threads.\n", .{@min(num_threads, 16)});

    var threads_idle_buf = [_]bool{true} ** 16;
    var threads_idle = threads_idle_buf[0..num_threads];

    var threads_working_mem_buf = [_][]u8{undefined} ** 16;
    var threads_working_mem = threads_working_mem_buf[0..num_threads];

    for (threads_working_mem) |*slc| {
        slc.* = alctr.alloc(u8, 1024 * 1024 * 256) catch unreachable;
    }

    defer for (threads_working_mem) |slc| {
        alctr.free(slc);
    };

    var threads_alctrs_buf = [_]std.heap.FixedBufferAllocator{undefined} ** 16;
    var threads_alctrs = threads_alctrs_buf[0..num_threads];

    for (threads_alctrs) |*alc, i| {
        alc.* = std.heap.FixedBufferAllocator.init(threads_working_mem[i]);
    }

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
                threads_idle[ti] = false;
                chunk.done = true;
                var thr = std.Thread.spawn(.{}, render, .{
                    world,             cam,                            qan,
                    prog,              threads_alctrs[ti].allocator(), chunk.start_x,
                    chunk.end_x,       chunk.start_y,                  chunk.end_y,
                    &threads_idle[ti],
                }) catch unreachable;
                thr.detach();
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

fn render(
    world: World,
    cam: Camera,
    qan: *Qanvas,
    prog: *std.Progress.Node,
    alctr: std.mem.Allocator,
    start_x: i64,
    end_x: i64,
    start_y: i64,
    end_y: i64,
    done_flag: *bool,
) void {
    defer done_flag.* = true;

    {
        var pix_y: i64 = start_y;
        while (pix_y < end_y) : (pix_y += 1) {
            var pix_x: i64 = start_x;
            while (pix_x < end_x) : (pix_x += 1) {
                qan.write(Color.init(1, 0, 1), pix_x, pix_y);
            }
        }
    }

    var pix_y: i64 = start_y;
    while (pix_y < end_y) : (pix_y += 1) {
        var pix_x: i64 = start_x;
        while (pix_x < end_x) : (pix_x += 1) {
            // 8 reflections
            qan.write(world.colorAt(cam.rayForPixel(pix_x, pix_y), alctr, 8), pix_x, pix_y);
            prog.completeOne();
        }
    }
}
