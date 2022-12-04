//! project specific math definitions

const std = @import("std");

/// tolerance for two floating point numbers to be considered the same
pub const floatTolerance = 512 * std.math.floatEps(f64);
