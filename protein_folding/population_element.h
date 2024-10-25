#pragma once
#include <functional>
#include <vector>
#include "amino_acid.h"
#include "Point.h"
#include <cfloat>


class population_element {
public:
	std::vector<float> x;
	float f = 0.5f;
	float cr = 0.9f;
	double e = DBL_MAX;

	population_element() = default;

	explicit population_element(const std::vector<float> &x) : x(x) {}


	static unsigned int number_of_energy_calculations;

	bool operator<(const population_element &other) const {
		return e < other.e;
	}

	bool operator>(const population_element &other) const {
		return e > other.e;
	}

	static std::vector<float> generate_x(int d, const std::function<float()> &randomizer);

	static std::vector<point> generate_p(int l, const std::vector<float> &x);

	static double calculate_e(int l, const std::vector<amino_acid> &s, const std::vector<float> &x, const std::vector<point> &p);

	static float c(const amino_acid &a, const amino_acid &b);
};

