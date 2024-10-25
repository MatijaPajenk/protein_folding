#pragma once
class point {
public:
	float x, y, z;
	point() = default;
	point(const float x, const float y, const float z) : x(x), y(y), z(z) {}

	[[nodiscard]] float distance(const point &p) const;
};

