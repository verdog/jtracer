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
    pixel_size: i64,
};

fn getChunks(w: i64, qan: Qanvas, alctr: std.mem.Allocator) []Chunk {
    var chunks = std.ArrayList(Chunk).init(alctr);

    const y_w = w;
    std.debug.assert(@popCount(w) == 1);

    const chunks_x = @divTrunc(@as(i64, @intCast(qan.width)) -| 1, w) + 1;
    const chunks_y = @divTrunc(@as(i64, @intCast(qan.height)) -| 1, y_w) + 1;

    var y: i64 = 0;
    while (y < chunks_y) : (y += 1) {
        var x: i64 = 0;
        while (x < chunks_x) : (x += 1) {
            chunks.append(.{
                .start_x = w * x,
                .end_x = @min(@as(i64, @intCast(qan.width)), w * (x + 1)),
                .start_y = y_w * y,
                .end_y = @min(@as(i64, @intCast(qan.height)), y_w * (y + 1)),
                .pixel_size = w,
            }) catch unreachable;
        }
    }

    const closeToCenter = struct {
        fn f(qan_ctx: Qanvas, lhs: Chunk, rhs: Chunk) bool {
            const center_left = .{
                .x = @divTrunc(lhs.start_x + lhs.end_x, 2),
                .y = @divTrunc(lhs.start_y + lhs.end_y, 2),
            };
            const center_right = .{
                .x = @divTrunc(rhs.start_x + rhs.end_x, 2),
                .y = @divTrunc(rhs.start_y + rhs.end_y, 2),
            };
            const center_canvas = .{
                .x = @as(i64, @intCast(@divTrunc(qan_ctx.width, 2))),
                .y = @as(i64, @intCast(@divTrunc(qan_ctx.height, 2))),
            };

            const left_distance = blk: {
                const dx = @as(f64, @floatFromInt(center_left.x - center_canvas.x));
                const dy = @as(f64, @floatFromInt(center_left.y - center_canvas.y));
                break :blk @sqrt(dx * dx + dy * dy);
            };
            const right_distance = blk: {
                const dx = @as(f64, @floatFromInt(center_right.x - center_canvas.x));
                const dy = @as(f64, @floatFromInt(center_right.y - center_canvas.y));
                break :blk @sqrt(dx * dx + dy * dy);
            };

            return left_distance < right_distance;
        }
    }.f;

    // sort chunks by distance from center to attempt to render the focus of the image first
    std.sort.block(Chunk, chunks.items, qan, closeToCenter);

    // TODO add this as a setting
    // var prng = std.rand.DefaultPrng.init(0);
    // var random = prng.random();
    // random.shuffle(Chunk, chunks.items);

    return chunks.toOwnedSlice() catch unreachable;
}

pub fn startRenderEngine(world: World, cam: Camera, qan: *Qanvas, alctr: std.mem.Allocator) void {
    std.debug.print("Starting render.\n", .{});
    var chunks = getChunks(32, qan.*, alctr);
    defer alctr.free(chunks);

    var prog_ctx = std.Progress{};
    var prog = blk: {
        var total_chunk_passes: u64 = 0;

        // assume all chunks have the same pixel size
        var working_size = chunks[0].pixel_size;
        while (working_size > 0) {
            total_chunk_passes += chunks.len;
            working_size = @divTrunc(working_size, 2);
        }

        break :blk prog_ctx.start("Chunks", total_chunk_passes);
    };

    const cpus = std.Thread.getCpuCount() catch unreachable;
    // const num_threads = @max(1, @divTrunc(cpus * 3, 4));
    const num_threads = cpus;
    std.debug.print("Using {} threads.\n", .{num_threads});

    const threads_buf_size = 32;
    var threads_idle_buf = [_]bool{true} ** threads_buf_size;
    std.debug.assert(threads_idle_buf.len >= num_threads);
    var threads_idle = threads_idle_buf[0..num_threads];

    var threads_working_mem_buf = [_][]u8{undefined} ** threads_buf_size;
    var threads_working_mem = threads_working_mem_buf[0..num_threads];

    for (threads_working_mem) |*slc| {
        slc.* = alctr.alloc(u8, 1024 * 1024 * 64) catch unreachable;
    }

    defer for (threads_working_mem) |slc| {
        alctr.free(slc);
    };

    var threads_alctrs_buf = [_]std.heap.FixedBufferAllocator{undefined} ** threads_buf_size;
    var threads_alctrs = threads_alctrs_buf[0..num_threads];

    for (threads_alctrs, 0..) |*alc, i| {
        alc.* = std.heap.FixedBufferAllocator.init(threads_working_mem[i]);
    }

    var timer = std.time.Timer.start() catch unreachable;

    while (true) {
        // find chunk that needs rendering
        const mchunk: ?*Chunk = blk: {
            var biggest_chunk: ?*Chunk = null;
            var biggest_pixel_size: i64 = 0;

            for (chunks) |*chunk| {
                if (chunk.pixel_size > biggest_pixel_size) {
                    biggest_chunk = chunk;
                    biggest_pixel_size = chunk.pixel_size;
                }
            }

            break :blk biggest_chunk;
        };

        // start thread on chunk
        if (mchunk) |chunk| {
            const mthread_i: ?usize = blk: {
                for (threads_idle, 0..) |ti, i| {
                    if (ti == true) {
                        break :blk i;
                    }
                }
                break :blk null;
            };

            if (mthread_i) |ti| {
                // start thread
                threads_idle[ti] = false;
                threads_alctrs[ti].reset();
                var thr = std.Thread.spawn(.{}, render, .{
                    world,            cam,                            qan,
                    prog,             threads_alctrs[ti].allocator(), chunk.start_x,
                    chunk.end_x,      chunk.start_y,                  chunk.end_y,
                    chunk.pixel_size, &threads_idle[ti],
                }) catch unreachable;
                thr.detach();

                chunk.pixel_size = chunk.pixel_size >> 1;
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
            for (threads_idle, 0..) |ti, i| {
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
    std.debug.print("Done in {d:.3} seconds.\n", .{@as(f64, @floatFromInt(timer.read())) / 1_000_000_000});
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
    pixel_size: i64,
    thread_done_flag: *bool,
) void {
    defer thread_done_flag.* = true;

    // draw dots
    const magenta = Color.init(1, 0, 1);
    qan.fill(start_x, start_y, 1, end_y - start_y, magenta);
    qan.fill(start_x, end_y - 1, end_x - start_x, 1, magenta);

    var pix_y: i64 = start_y;
    while (pix_y < end_y) : (pix_y += pixel_size) {
        var pix_x: i64 = start_x;
        while (pix_x < end_x) : (pix_x += pixel_size) {
            // 5 reflections
            const color = world.colorAt(cam.rayForPixel(pix_x, pix_y), alctr, 5);
            qan.fill(pix_x, pix_y, pixel_size, pixel_size, color);
        }
    }

    prog.completeOne();
}
