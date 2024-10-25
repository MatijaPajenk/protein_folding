#include "point.h"
#include <cmath>

float point::distance(const point &p) const {
	const float dx = x - p.x;
	const float dy = y - p.y;
	const float dz = z - p.z;
	return sqrtf(dx * dx + dy * dy + dz * dz);
}
