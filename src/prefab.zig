//! Functions to apply common properties to objects. Because objects are closely tied
//! to the World, these functions modify objects instead of creating and returning them.

const vol = @import("volume.zig");
const mat = @import("material.zig");

pub fn toGlass(mate: *mat.Material) void {
    mate.transparency = 1.0;
    mate.refractive_index = 1.5;
}
