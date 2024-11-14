// ReSharper disable CppTooWideScopeInitStatement
#include "amino_acid.h"
#include "population_element.h"
#include <algorithm>
#include <chrono>
#include <exception>
#include <functional>
#include <iostream>
#include <map>
#include <random>
#include <string>
#include <vector>
#include <cmath>

constexpr double pi = 3.141592653589793;

namespace {
	std::vector<amino_acid> parse_amino_acid_string(const std::string &s) {
		std::vector<amino_acid> amino_acids(s.size());
		for (int i = 0; i < static_cast<int>(s.size()); ++i) {
			switch (s[i]) {
				case 'A':
					amino_acids[i] = amino_acid::a;
					break;
				case 'B':
					amino_acids[i] = amino_acid::b;
					break;
				default:
					std::cerr << "Error: Invalid amino acid character." << "\n";
					return {};
			}
		}
		return amino_acids;
	}

	std::map<std::string, std::string> parse_arguments(const int argc, char *argv[]) {
		std::map<std::string, std::string> args;
		for (int i = 2; i < argc; i += 2) {
			std::string key = argv[i];

			// Ensure there's a value after the key
			if (i + 1 < argc) {
				const std::string value = argv[i + 1];
				args[key] = value;
			} else {
				std::cerr << "Error: No value provided for " << key << '\n';
				return {};
			}
		}
		return args;
	}

	std::vector<population_element> initialize_population(const unsigned int np, const int l, const std::function<float()> &randomizer) {
		std::vector<population_element> population(np);
		for (int i = 0; i < static_cast<int>(np); ++i) {
			population[i] = population_element(population_element::generate_x(2 * l - 5, randomizer));
		}
		return population;
	}

	bool check_argument_count(const int argc, char *argv[]) {
		if (argc != 12) {
			std::cerr << "Usage: " << *argv
				<< "<program> <vector<'a'|'b'>> -seed <unsigned int> -target <float> -nfesLmt <unsigned int> -runtimeLmt <unsigned int> -Np <unsigned int>"
				<< "\n";
			return false;
		}
		return true;
	}
}

int main(const int argc, char *argv[]) {
	if (!check_argument_count(argc, argv)) return 1;

	const std::string amino_acid_string = argv[1];
	const int l = static_cast<int>(amino_acid_string.size());
	const int d = 2 * l - 5;

	const std::vector<amino_acid> s = parse_amino_acid_string(amino_acid_string);
	if (s.empty()) return 1;

	std::map<std::string, std::string> args = parse_arguments(argc, argv);
	if (args.empty()) return 1;

	bool target_reached = false, nfes_limit_reached = false, runtime_limit_reached = false;
	unsigned int seed = std::stoul(args["-seed"]);
	const double target = std::stof(args["-target"]);
	const unsigned int nfes_lmt = std::stoul(args["-nfesLmt"]);
	const unsigned int runtime_lmt = std::stoul(args["-runtimeLmt"]);
	const unsigned int np = std::stoul(args["-Np"]);

	std::mt19937 gen(seed);

	std::uniform_real_distribution<float> float_dist(0.0f, 1.0f);
	auto float_randomizer = [&]() { return float_dist(gen); };

	std::uniform_int_distribution<int> int_dist_np(0, static_cast<int>(np - 1));
	auto int_randomizer_np = [&]() { return int_dist_np(gen); };

	std::uniform_int_distribution<int> int_dist_d(0, d - 1);
	auto int_randomizer_d = [&]() { return int_dist_d(gen); };

	std::vector<population_element> population = initialize_population(np, l, float_randomizer);

	int best_index = static_cast<int>(std::distance(population.begin(), std::min_element(population.begin(), population.end())));

	int r1, r2;

	const auto start = std::chrono::high_resolution_clock::now();
	const auto end = start + std::chrono::seconds(runtime_lmt);
	unsigned int number_of_energy_calculations = np;

	while (!nfes_limit_reached && !target_reached && !runtime_limit_reached) {
		for (int i = 0; i < static_cast<int>(np); ++i) {
			const float f = float_randomizer() < 0.1f ? 0.1f + 0.9f * float_randomizer() : population[i].f;
			const float cr = float_randomizer() < 0.1f ? float_randomizer() : population[i].cr;

			do { r1 = int_randomizer_np(); } while (r1 == i || r1 == best_index);
			do { r2 = int_randomizer_np(); } while (r2 == i || r2 == r1 || r2 == best_index);

			const int j_rand = int_randomizer_d();

			std::vector<float> u(d);
			for (int j = 0; j < d; ++j) {
				if (j == j_rand || float_randomizer() < cr) {
					u[j] = population[best_index].x[j] + f * (population[r1].x[j] - population[r2].x[j]);
					if (u[j] <= -pi) u[j] += static_cast<float>(2 * pi);
					if (u[j] > pi) u[j] += static_cast<float>(2 * -pi);
				} else {
					u[j] = population[i].x[j];
				}
			}

			const double eu = population_element::calculate_e(l, s, u, population_element::generate_p(l, u));
			number_of_energy_calculations++;

			if (eu <= population[i].e) {
				std::vector<float> us(d);
				for (int j = 0; j < d; ++j) {
					us[j] = population[best_index].x[j] + 0.5f * (u[j] - population[i].x[j]);
					if (us[j] <= -pi) us[j] += static_cast<float>(2 * pi);
					if (us[j] > pi) us[j] += static_cast<float>(2 * -pi);
				}

				const double eus = population_element::calculate_e(l, s, us, population_element::generate_p(l, us));
				number_of_energy_calculations++;
				population[i].x = eus <= eu ? us : u;
				population[i].e = eus <= eu ? eus : eu;
				population[i].f = f;
				population[i].cr = cr;
			}

			if (population[i].e < population[best_index].e) best_index = i;

			if (number_of_energy_calculations >= nfes_lmt) nfes_limit_reached = true;
			if (end < std::chrono::high_resolution_clock::now()) runtime_limit_reached = true;
			if (population[i].e <= target + 1e-6) target_reached = true;
		}
	}
	const auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start).count();

	std::cout << seed << ";" << duration << ";" << number_of_energy_calculations / duration * 1000 << ";" << population[best_index].e << ";" << number_of_energy_calculations << "\n";

	//std::ranges::for_each(population[best_index].x.begin(), population[best_index].x.end(), [](const float x) {std::cout << x * 180 / pi << " "; });

	return 0;
}
