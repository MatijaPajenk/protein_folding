#pragma once
#include <algorithm>
#include <iostream>
#include <vector>

class test_angles_util {
public:
	static const std::vector<float> seed_66;
	static const std::vector<float> actual_solution;

	static std::vector<float> calculate_diff(const std::vector<float> &a, const std::vector<float> &b) {
		if (a.size() != b.size()) {
			std::cerr << "Vectors must have the same size!\n";
			return {};
		}

		std::vector<float> differences(a.size());

		std::transform(a.begin(), a.end(), b.begin(), differences.begin(), [](float a, float b) {
			const float diff = a - b; // Calculate the difference
			float normalized = std::fmod(diff + 360.0f, 360.0f); // Normalize to range [0, 360)

			// If the normalized difference is greater than 180, choose the shorter path (negative equivalent)
			if (normalized > 180.0f) {
				normalized -= 360.0f; // Make it negative for shorter distance
			}

			return normalized; }
		);

		return differences;
	}
};

const std::vector<float> test_angles_util::actual_solution = {
	43.2915f, 2.88166f, -48.728f, 0.0655009f, 12.6242f, 66.0927f, -6.40805f, 8.96332f,8.80015f, 2.23544f, 74.0763f,
	-6.62061f, 1.31798f, -104.099f, 160.341f, -177.384f, -20.6892f, 26.8003f, 127.789f, 166.27f, 10.2979f
};

const std::vector<float> test_angles_util::seed_66 = {
	24.6194f, -61.6804f, -23.269f, 26.4055f, 43.21f, 76.6225f, 7.10225f, 30.688f, 3.95914f, 60.5764f, -21.8314f,
	129.854f, 17.9836f, 126.889f, -134.559f, -177.388f, 42.5306f, -14.9708f, -106.605f, -166.419f, 40.5914f
};



