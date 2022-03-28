

#include "cxxopts.hpp"
#define ANKERL_NANOBENCH_IMPLEMENT
#include "nanobench.h"

#include <pagmo/algorithm.hpp>
#include <pagmo/algorithms/nsga2.hpp>
#include <pagmo/batch_evaluators/thread_bfe.hpp>
#include <pagmo/population.hpp>
#include <pagmo/problem.hpp>
#include <pagmo/problems/dtlz.hpp>
#include <pagmo/utils/nondominated_sorting.hpp>

#include <tbb/global_control.h>
#include <tbb/parallel_for.h>

#include <iostream>
#include <thread>
#include <unordered_map>

using namespace pagmo;

namespace nb = ankerl::nanobench;

int main(int argc, char** argv)
{
    cxxopts::Options opts("nds", "empirical tests of non-dominated sorting algorithms using the NSGA-II algorithm.");

    opts.add_options()
        ("t,threads", "the number of threads to use", cxxopts::value<size_t>()->default_value(std::to_string(std::thread::hardware_concurrency())))
        ("n,population-size", "the population size", cxxopts::value<size_t>()->default_value("1000"))
        ("m,objectives", "the number of objectives", cxxopts::value<size_t>()->default_value("2"))
        ("d,problem-dimension", "the problem dimension", cxxopts::value<size_t>()->default_value("10"))
        ("e,epochs", "the number of epochs (runs) during benchmarking", cxxopts::value<size_t>()->default_value("5"))
        ("g,generations", "the number of generations", cxxopts::value<size_t>()->default_value("500"))
        ("s,seed", "the random seed", cxxopts::value<size_t>())
        ("r,reps", "number of repetions of a benchmark", cxxopts::value<size_t>()->default_value("10"))
        ("a,algorithm", "the non-dominated sorting algorithm to use", cxxopts::value<std::string>()->default_value("RS"))
        ("p,problem", "the id of the problem e.g. 1 for DTLZ1, 2 for DTLZ2,...", cxxopts::value<size_t>()->default_value("1"));

    const std::unordered_map<std::string, non_dominated_sorting_algorithm> algorithms {
        { "DS", non_dominated_sorting_algorithm::DS },
        { "ENS-BS", non_dominated_sorting_algorithm::ENS_BS },
        { "ENS-SS", non_dominated_sorting_algorithm::ENS_SS },
        { "HS", non_dominated_sorting_algorithm::HS },
        { "MNDS", non_dominated_sorting_algorithm::MNDS },
        { "RS", non_dominated_sorting_algorithm::RS },
        { "RO", non_dominated_sorting_algorithm::RO }
    };

    auto result = opts.parse(argc, argv);
    size_t problem_dim = result["problem-dimension"].as<size_t>();
    size_t num_threads = result["threads"].as<size_t>();
    size_t n = result["population-size"].as<size_t>();
    size_t m = result["objectives"].as<size_t>();
    size_t epochs = result["epochs"].as<size_t>();
    size_t generations = result["generations"].as<size_t>();
    size_t reps = result["reps"].as<size_t>();
    size_t seed = pagmo::random_device::next();
    if (result["seed"].count() > 0) {
        seed = result["seed"].as<size_t>();
    }

    std::string algorithm_name = result["algorithm"].as<std::string>();
    size_t problem_id = result["problem"].as<size_t>();

    auto nds_alg = algorithms.find(algorithm_name);
    if (nds_alg == algorithms.end()) {
        throw std::runtime_error("unknown algorithm");
    }

    //tbb::global_control c(tbb::global_control::max_allowed_parallelism, num_threads); // in case we use a threaded bfe
    nsga2 algo(
            generations, // generations
            0.95,        // crossover rate
            10.,         // eta_c (distribution index for crossover)
            0.1,         // mutation probability
            50.,         // eta_m (distribution index for mutation),
            seed,        // seed
            nds_alg->second, // the algorithm for non-dominated sorting
            false        // dominate on equal
            );

    //algo.set_bfe(bfe{thread_bfe{}});

    nb::Bench bench;
    problem prob{dtlz(problem_id, problem_dim, m)};
    population pop{prob, n, seed};
    std::string run_label = "DTLZ" + std::to_string(problem_id) + ", " + algorithm_name + ", " + std::to_string(n) + ", " + std::to_string(m);
    bench.output(nullptr);
    for (size_t rep = 0; rep < reps; ++rep) {
        bench.epochIterations(epochs).run(run_label, [&]() {
            algo.evolve(pop);
            //std::cout << "The population: \n" << pop;
        });
    }
    //std::cout << "seed: " << seed << "\n";
    bench.render(ankerl::nanobench::templates::csv(), std::cout);

    // 5 - Output the population
    std::cout << "The population: \n" << pop;
}
