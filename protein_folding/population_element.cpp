#include "amino_acid.h"
#include "point.h"
#include "population_element.h"
#include <functional>
#include <vector>

constexpr double pi = 3.141592653589793;

std::vector<float> population_element::generate_x(const int d, const std::function<float()> &randomizer) {
	std::vector<float> x(d);
	for (int i = 0; i < d; ++i) {
		x[i] = static_cast<float>(-pi + 2 * pi * randomizer());
		if (x[i] == 0.f) x[i] = 1e-6f;
	}
	return x;
}

std::vector<point> population_element::generate_p(const int l, const std::vector<float> &x) {
	std::vector<point> p(l);
	p[0] = point(0, 0, 0);
	p[1] = point(0, 1, 0);
	p[2] = point(cosf(x[0]), 1.0f + sinf(x[0]), 0);

	const int beta = l - 2;
	for (int i = 3; i < l; ++i) {
		p[i] = point(
			p[i - 1].x + cosf(x[i - 2]) * cosf(x[beta + i - 3]),
			p[i - 1].y + sinf(x[i - 2]) * cosf(x[beta + i - 3]),
			p[i - 1].z + sinf(x[beta + i - 3])
		);
	}
	return p;
}

float population_element::c(const amino_acid &a, const amino_acid &b) {
	if (a == amino_acid::a && b == amino_acid::a) return 1.0f;
	if (a == amino_acid::b && b == amino_acid::b) return 0.5f;
	return -0.5f;
}

double population_element::calculate_e(const int l, const std::vector<amino_acid> &s, const std::vector<float> &x, const std::vector<point> &p) {
	double energy_1 = 0;

	for (int i = 0; i < l - 2; ++i) energy_1 += 1 - cosf(x[i]);
	energy_1 *= 0.25;

	double energy_2 = 0;
	for (int i = 0; i < l - 2; ++i) {
		for (int j = i + 2; j < l; ++j) {
			const float d = p[i].distance(p[j]);
			const float d2 = d * d;
			const float d6 = d2 * d2 * d2;
			const float inv_d6 = 1 / d6;
			energy_2 += (inv_d6 * inv_d6 - c(s[i], s[j]) * inv_d6);
		}
	}
	return energy_1 + 4 * energy_2;
}
