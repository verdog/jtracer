const std = @import("std");
const SDLSdk = @import("external/SDL.zig/Sdk.zig");

pub fn build(b: *std.build.Builder) void {
    const optimize = b.standardOptimizeOption(.{
        .preferred_optimize_mode = .ReleaseSafe,
    });

    const SDL = SDLSdk.init(b, null);

    const exe = b.addExecutable(.{
        .name = "jtracer",
        .root_source_file = .{ .path = "src/main.zig" },
        .optimize = optimize,
    });
    exe.emit_docs = .emit;
    exe.setMainPkgPath("./src");

    SDL.link(exe, .dynamic);
    exe.addModule("sdl2", SDL.getWrapperModule());
    exe.linkSystemLibrary("sdl2_image");

    b.installArtifact(exe);

    const run_cmd = b.addRunArtifact(exe);
    run_cmd.step.dependOn(b.getInstallStep());
    if (b.args) |args| {
        run_cmd.addArgs(args);
    }

    const run_step = b.step("run", "Run the app");
    run_step.dependOn(&run_cmd.step);

    const test_step = b.step("test", "Run unit tests");
    var ts = b.addTest(.{
        .root_source_file = .{ .path = "src/main.zig" },
        .optimize = optimize,
    });
    test_step.dependOn(&b.addRunArtifact(ts).step);
}
