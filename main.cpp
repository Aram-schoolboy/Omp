#include "hit.h"
#include <iostream>
#include <omp.h>
#include <random>
#include <cstdlib>
#include <ctime>


uint32_t max = std::numeric_limits<uint32_t>::max() ;
int threads = 0;
size_t total_hit = 0;
float answer = 0;

struct Options {
    int threads = 0;
    bool threads_default = false;
    const char* output_file_name;
    size_t n = 1;
};

int ParseInt(const std::string& str) {
    int num = 0;
    for (int i = 0; i < str.size(); ++i) {
        num += str[i] - '0';
        num *= 10;
    }
    return num / 10;
}

bool Parse(int argc, char* argv[], Options& options) {
    if (argc < 2 || argc > 7) {
        std::cerr << "invalid count of args";
        return false;
    }

    uint8_t arg_index = 1;
    bool omp, input_file, output_file;
    omp = input_file = output_file = false;
    while (arg_index < argc) {
        if (static_cast<std::string>(argv[arg_index]) == "--no-omp") {
            if (omp) {
                std::cerr << "omp option should be given only once";
                return false;
            }
            omp = true;
            options.threads = 0;
            ++arg_index;
        } else if (static_cast<std::string>(argv[arg_index]) == "--omp-threads") {
            if (arg_index + 1 == argc) {
                std::cerr << "write threads number";
                return false;
            }
            if (omp) {
                std::cerr << "omp option should be given only once";
                return false;
            }
            omp = true;
            ++arg_index;
            if (static_cast<std::string>(argv[arg_index]) == "default") {
                options.threads_default = true;
            } else {
                options.threads = ParseInt(static_cast<std::string>(argv[arg_index]));
            }
            ++arg_index;
        } else if (static_cast<std::string>(argv[arg_index]) == "--input") {
            if (arg_index + 1 == argc) {
                std::cerr << "write input file name";
                return false;
            }
            if (input_file) {
                std::cerr << "input file option should be given only once";
                return false;
            }
            input_file = true;
            ++arg_index;
            FILE* file = fopen(argv[arg_index], "r");
            fscanf(file, "%llu\n", &options.n);
            fclose(file);
            ++arg_index;
        } else if (static_cast<std::string>(argv[arg_index]) == "--output") {
            if (arg_index + 1 == argc) {
                std::cerr << "write output file name";
                return false;
            }
            if (output_file) {
                std::cerr << "output file option should be given only once";
                return false;
            }
            output_file = true;
            ++arg_index;
            options.output_file_name = argv[arg_index];
            ++arg_index;
        } else {
            std::cerr << "invalid argument";
            return false;
        }
    }

    return (omp && input_file && output_file);
}

void ParseAXIS_RANGE(float& x_min, float& x_max, float& y_min, float& y_max, float& z_min, float& z_max,
                     float& x_range, float& y_range, float& z_range) {
    const float* p = get_axis_range();
    x_min = p[0];
    x_max = p[1];
    y_min = p[2];
    y_max = p[3];
    z_min = p[4];
    z_max = p[5];
    x_range = x_max - x_min;
    y_range = y_max - y_min;
    z_range = z_max - z_min;
}

void WriteAnswer(const char* file_name) {
    FILE* file = fopen(file_name, "w");
    fprintf(file, "%g\n", answer);
    fclose(file);
}

void CalculateAnswer(size_t n, float x_range, float y_range, float z_range, float x_min, float y_min, float z_min, bool omp) {
    #pragma omp parallel if(omp) shared(x_range, y_range, z_range, x_min, y_min, z_min)
    {
        #pragma omp single
        threads = omp_get_num_threads();

        size_t hit = 0;
        unsigned int seed = (time(nullptr) - omp_get_thread_num()) % 2143;
        std::mt19937 generator{seed};

        #pragma omp for schedule(static)
        for (size_t i = 0; i < n; ++i) {
            if (hit_test((float) generator() / max * x_range + x_min, (float) generator() / max * y_range + y_min, (float) generator() / max * z_range + z_min)) { ++hit; }
        }
        #pragma omp critical
        total_hit += hit;
    }
    answer = x_range * y_range * z_range * total_hit / n;
}

int main(int argc, char* argv[]) {
    Options options;
    if (!Parse(argc, argv, options)) { return 1; }

    size_t n = options.n;
    threads = options.threads;
    const char* output_file_name = options.output_file_name;

    float x_min, x_max, y_min, y_max, z_min, z_max, x_range, y_range, z_range;
    ParseAXIS_RANGE(x_min, x_max, y_min, y_max, z_min, z_max, x_range, y_range, z_range);

    bool omp = false;
    if (options.threads_default) {
        omp = true;
    } else if (threads > 0) {
        omp_set_num_threads(threads);
        omp = true;
    }

    double start = omp_get_wtime();
    CalculateAnswer(n, x_range, y_range, z_range, x_min, y_min, z_min, omp);
    double end = omp_get_wtime();

    WriteAnswer(output_file_name);

    if (!omp) { threads = 0; }
    printf("Time (%i thread(s)): %g ms\n", threads, (end - start) * 1000);
}