#include "hit.h"

const int A = 2;
const int A_squared = A * A;
const float y_z_max_module = 0.649519052838329;
const float AXIS_RANGE[] = {0, 2, -y_z_max_module, y_z_max_module, -y_z_max_module, y_z_max_module};

bool hit_test(float x, float y, float z) {
    return (x * x * x * (x - A) <= -A_squared * (y * y + z * z));
}

const float* get_axis_range() {
    return AXIS_RANGE;
} 