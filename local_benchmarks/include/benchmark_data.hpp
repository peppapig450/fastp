#pragma once

#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <memory>
#include <optional>
#include <ostream>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <thread>
#include <type_traits>
#include <utility>

#include "evaluator.h"
#include "fastqreader.h"
#include "fragment_config.hpp"
#include "knownadapters.h"
#include "read.h"

namespace fs = std::filesystem;

#if defined(__linux__)
    #include <sched.h>
    #include <unistd.h>
#elif defined(__APPLE__)
    #include <sys/sysctl.h>
    #include <sys/types.h>
#elif defined(_WIN32)
    #include <windows.h>
#endif

namespace benchmark_data {
namespace cpu {

// Count set bits for affinity masks (portable-ish)
inline auto popcount_ull(std::uint64_t x) -> unsigned {
    // Simple fallback, fast enough for this usage
    unsigned c = 0;
    while (x != 0U) {
        // 'x & (x - 1)' clears the least significant bit (LSB) that is set to 1.
        // Example: 1100 (12) becomes 1000 (8).
        // The loop runs exactly as many times as there are set bits.
        x &= (x - 1);
        ++c;
    }

    return c;
}

inline auto affinity_concurrency() -> unsigned {
#if defined(__linux__)
    // Respect current process's allowed CPUs (taskset/cgroups cpuset)
    cpu_set_t set;
    CPU_ZERO(&set);

    if (sched_getaffinity(0, sizeof(set), &set) == 0) {
        unsigned logical_cores = 0U;
        for (int i = 0; i < CPU_SETSIZE; ++i) {
            if (CPU_ISSET(i, &set)) {
                ++logical_cores;
            }
        }

        if (logical_cores != 0U) {
            return logical_cores;
        }
    }

    // Fallback: online logical CPUs
    std::uint32_t logical_cores = sysconf(_SC_NPROCESSORS_ONLN);
    return (logical_cores > 0U) ? static_cast<unsigned>(logical_cores) : 0U;

#elif defined(__APPLE__)
    // Wrapper around BSD's libc method for retrieving system information
    auto get_apple_cpu_count = [](const char* name) -> unsigned {
        int         logical_cores = 0;
        std::size_t len           = sizeof(logical_cores);

        if (sysctlbyname(name, &logical_cores, &len, nullptr, 0) == 0 && logical_cores > 0) {
            return static_cast<unsigned>(logical_cores);
        }

        return 0U;
    };

    // Active logical CPUs (respects thermal/power state)
    if (auto logical_cores = get_apple_cpu_count("hw.activecpu"); logical_cores > 0) {
        return logical_cores;
    }

    // Fallback: total logical CPUs
    if (auto logical_cores = get_apple_cpu_count("hw.logicalcpu"); logical_cores > 0) {
        return logical_cores;
    }

    return 0U;

#elif defined(_WIN32)
    // Microsoft why do you have to name your C methods like they're Java methods???
    // Best effort: count bits in process affinity across all groups
    DWORD len = 0;
    if (!GetProcessGroupAffinity(GetCurrentProcess(), &len, nullptr)
        && GetLastError() == ERROR_INSUFFICIENT_BUFFER) {
        std::vector<GROUP_AFFINITY> groups(len);

        if (GetProcessGroupAffinity(GetCurrentProcess(), &len, groups.data())) {
            unsigned total_cores = 0;
            for (DWORD i = 0; i < len; ++i) {
                total_cores += popcount_ull(groups[i].Mask);
            }

            if (total_cores != 0U) {
                return total_cores;
            }
        }
    }

    // Fallback that handles processor groups
    if (WORD count = GetActiveProcessorCount(ALL_PROCESSOR_GROUPS); count != 0U) {
        return count;
    }

    // Oldest fallback
    SYSTEM_INFO si {};
    GetSystemInfo(&si);
    return si.dwNumberOfProcessors;
#else
    return 0U;
#endif
}

// Use a cached value with a safe fallback to std::thread::hardware_concurrency()
inline auto hardware_concurrency() -> unsigned {
    static const auto cached = [] {
        auto logical_cores = affinity_concurrency();
        if (logical_cores == 0U) {
            logical_cores = std::thread::hardware_concurrency();
        }

        return std::max(1U, logical_cores);
    }();

    return cached;
}
}  // namespace cpu

namespace {  // anon

#ifdef _WIN32
constexpr const char* path_separator = ";";
#else
constexpr const char* path_separator = ":";
#endif

auto is_executable(const fs::path& path) {
#ifdef _WIN32
    return fs::exists(path);
#else
    return fs::exists(path)
           && (fs::status(path).permissions() & fs::perms::owner_exec) != fs::perms::none;
#endif
}

auto is_command_available(std::string_view command) {
    const char* path_env = std::getenv("PATH");
    if (path_env == nullptr) {
        return false;
    }

    std::string_view path_str(path_env);
    std::size_t      start = 0;
    while (start < path_str.size()) {
        std::size_t end = path_str.find(path_separator, start);
        if (end == std::string_view::npos) {
            end = path_str.size();
        }

        fs::path dir = path_str.substr(start, end - start);
        if (is_executable(dir / command)) {
            return true;
        }
#ifdef _WIN32
        if (is_executable(dir / (std::string(command) + ".exe"))) {
            return true;
        }
#endif
        start = end + 1;
    }

    return false;
}

auto getAdapterSequences() -> std::vector<std::string> {
    const auto&              adapters = adapters::getKnown();
    std::vector<std::string> adaptersSeqs;
    adaptersSeqs.reserve(adapters.size());

    for (const auto& adapterPair : adapters) {
        adaptersSeqs.push_back(adapterPair.first);
    }

    return adaptersSeqs;
}

void applyAdapterContamination(std::string&     sequence,
                               std::string&     quality,
                               std::string_view adapter,
                               std::mt19937&    rng) {
    constexpr std::size_t minInsert          = 15;
    constexpr std::size_t kMinAdapterOverlap = 5;
    const auto            minAdapter = std::min<std::size_t>(kMinAdapterOverlap, adapter.size());

    if (sequence.size() > minInsert + minAdapter) {
        std::uniform_int_distribution<std::size_t> cutDist(minInsert, sequence.size() - minAdapter);
        const std::size_t                          cut = cutDist(rng);

        const auto adapterToTake = std::min<std::size_t>(sequence.size() - cut, adapter.size());

        // Keep first 'cut' bases and then append adapter
        sequence.resize(cut);
        sequence.append(adapter.substr(0, adapterToTake));

        // Keep first 'cut' quals, then append 'I' chars
        quality.resize(cut);
        quality.append(adapterToTake, 'I');
    }
}

void contaminateStream(std::istream& input,
                       std::ostream& output,
                       double        contamRate,
                       unsigned int  rngSeed = 669) {
    const auto adapterSeqs = getAdapterSequences();

    std::mt19937                           rng(rngSeed);
    std::uniform_real_distribution<double> contamDist {0.0, 1.0};

    // Only define adapterChoice if we actually have adapters (should never happen since adapters
    // are hardcoded, but just in case)
    std::optional<std::uniform_int_distribution<int>> adapterChoice;
    if (!adapterSeqs.empty()) {
        adapterChoice.emplace(0, static_cast<int>(adapterSeqs.size() - 1));
    }

    std::string h;
    std::string s;
    std::string p;
    std::string q;
    while (std::getline(input, h) && std::getline(input, s) && std::getline(input, p)
           && std::getline(input, q)) {
        if (adapterChoice && s.size() && q.size() && contamDist(rng) < contamRate) {
            applyAdapterContamination(s, q, adapterSeqs[(*adapterChoice)(rng)], rng);
        }

        output << h << '\n' << s << '\n' << p << '\n' << q << '\n';
    }
}

template <typename TDerived, typename... TBase>
constexpr bool is_any_base_of_v = (... || std::is_base_of_v<TBase, TDerived>);

template <typename FstreamType>
concept is_fstream = is_any_base_of_v<FstreamType, std::fstream, std::ifstream, std::ofstream>;

template <is_fstream StreamType>
auto openFile(const std::string& path) -> StreamType {
    StreamType file(path);
    if (!file) {
        throw std::runtime_error("Cannot open: " + path);
    }

    return file;
}

void contamOneFile(const std::string& inFastq, const std::string& outFastq, double contamRate) {
    auto input  = openFile<std::ifstream>(inFastq);
    auto output = openFile<std::ofstream>(outFastq);
    contaminateStream(input, output, contamRate);
}

// Helper function for running a command
void runCommand(const std::string& cmd, const std::string& errorMessage) {
    if (std::system(cmd.c_str()) != 0) {
        throw std::runtime_error(errorMessage + ": " + cmd);
    }
}

}  // namespace

struct MasonAvailability {
    bool simulator   = false;
    bool variator    = false;
    bool materialize = false;
};

// Generate realistic test data
class MasonDataGenerator {
public:
    MasonDataGenerator()
        : cpu_count_ {cpu::hardware_concurrency()} {
        ensureMasonAvailable();
        setupWorkDirectory();
    }

    ~MasonDataGenerator() { cleanup(); }

    // Single-end
    // 5-step pipeline (ref -> [variants] -> materialize -> simulate -> contam) -> FASTQ
    // if refFasta is empty, we synthesize a 1Mb reference
    auto generateDataset(const std::string& name,
                         std::size_t        numReads,
                         std::size_t        readLength,
                         const std::string& platform          = "illumina",
                         double             adapterContamRate = 0.1,
                         const std::string& errorModel        = "hiseq",
                         bool               withVariants      = false,
                         std::string_view   refFasta          = {}) -> std::string {
        return generateDatasetImpl<false>(name,
                                          numReads,
                                          readLength,
                                          platform,
                                          adapterContamRate,
                                          errorModel,
                                          withVariants,
                                          refFasta);
    }

    // 5-step pipeline -> pair of FASTQs (R1, R2), optional external reference
    auto generateDatasetPE(const std::string& name,
                           std::size_t        numReads,
                           std::size_t        readLength,
                           const std::string& platform          = "illumina",
                           double             adapterContamRate = 0.1,
                           const std::string& errorModel        = "hiseq",
                           bool               withVariants      = false,
                           std::string_view refFasta = {}) -> std::pair<std::string, std::string> {
        return generateDatasetImpl<true>(name,
                                         numReads,
                                         readLength,
                                         platform,
                                         adapterContamRate,
                                         errorModel,
                                         withVariants,
                                         refFasta);
    }

    auto loadReadsFromFastq(const std::string& fastqPath, std::size_t maxReads = 0)
        -> std::vector<std::unique_ptr<Read>> {
        std::vector<std::unique_ptr<Read>> reads;
        FastqReader                        reader {fastqPath};

        std::size_t count = 0;
        while (!reader.eof()) {
            auto read = std::unique_ptr<Read> {reader.read()};
            if (read == nullptr) {
                break;
            }

            reads.push_back(std::move(read));
            if (maxReads > 0 && ++count >= maxReads) {
                break;
            }
        }

        return reads;
    }

    void setFragmentConfig(std::optional<FragmentConfig> fragConfig) { fragment_cfg_ = fragConfig; }

private:
    fs::path                      workDir_;
    unsigned                      cpu_count_;
    MasonAvailability             mason_;
    std::optional<FragmentConfig> fragment_cfg_;

    template <bool IsPairedEnd>
    using DatasetReturn = std::conditional_t<IsPairedEnd,
                                             std::pair<std::string, std::string>,  // PE
                                             std::string>;                         // SE

    // Forward declare so the compiler is aware of the return type
    template <bool IsPairedEnd>
    auto generateDatasetImpl(const std::string& name,
                             std::size_t        numReads,
                             std::size_t,
                             const std::string& platform,
                             double             adapterContamRate,
                             const std::string& errorModel,
                             bool               withVariants,
                             std::string_view   refFasta) -> DatasetReturn<IsPairedEnd>;

    void ensureMasonAvailable() {
        mason_.simulator   = is_command_available("mason_simulator");
        mason_.variator    = is_command_available("mason_variator");
        mason_.materialize = is_command_available("mason_materializer");

        if (!mason_.simulator) {
            std::cerr << "Warning: mason_simulator not found.\n";
        }
        if (!mason_.variator) {
            std::cerr << "Note: mason_variator not found (variants will be skipped).\n";
        }
        if (!mason_.materialize) {
            std::cerr << "Note: mason_materializer not found (materialization will be skipped).\n";
        }
        if (!mason_.simulator && !mason_.variator && !mason_.materialize) {
            std::cerr
                << "Install via Conda:\n"
                << "  conda install -c bioconda mason\n"
                << "Or build from source: https://github.com/seqan/seqan/tree/main/apps/mason2\n";
        }
    }

    void setupWorkDirectory() {
        workDir_ = fs::temp_directory_path() / "fastp_benchmark_data";
        fs::create_directories(workDir_);
    }

    void cleanup() {
        // We could cleanup generated files here if we want
        // fs::remove_all(workDir_);
    }

    // If user passes an external FASTA, ensure it exists and looks like FASTA
    auto validateReferenceGenome(std::string_view refFasta) -> fs::path {
        fs::path path {refFasta};

        if (!fs::exists(path)) {
            throw std::runtime_error("Reference FASTA not found: " + path.string());
        }

        // very light sanity check: first char '>
        std::ifstream in(path);
        char          firstChar = 0;
        if (!(in >> firstChar) || firstChar != '>') {
            throw std::runtime_error("Reference FASTA does not appear to be in FASTA format: "
                                     + path.string());
        }

        return path;
    }

    // Synthesize a 1 Mbp contig with ~42% GC, name-tagged to avoid clashes
    auto generateReferenceGenome(const std::string& tag) -> std::string {
        constexpr int    genomeLength = 1000000;
        constexpr double gcContent    = 0.42;  // ~21% GC, 21% C
        constexpr double atContent    = 0.58;  // ~29% A, 29% T
        constexpr int    lineWrap     = 80;

        // Probabilities for each base based on the target GC/AT content
        constexpr double gProb = gcContent / 2.0;            // 0.21
        constexpr double cProb = gProb + (gcContent / 2.0);  // 0.42
        constexpr double aProb = cProb + (atContent / 2.0);  // 0.71

        constexpr unsigned int rngSeed = 420;

        auto refPath = workDir_ / (tag + "_reference.fasta");
        if (!fs::exists(refPath)) {
            std::ofstream ref(refPath);
            ref << ">chr1\n";

            std::mt19937                           rng(rngSeed);
            std::uniform_real_distribution<double> uniformDistribution(0.0, 1.0);
            for (int i = 0; i < genomeLength; ++i) {
                auto randomValue = uniformDistribution(rng);
                char base;
                if (randomValue < gProb) {
                    base = 'G';
                } else if (randomValue < cProb) {
                    base = 'C';
                } else if (randomValue < aProb) {
                    base = 'A';
                } else {
                    base = 'T';
                }

                ref << base;
                if ((i + 1) % lineWrap == 0) {
                    ref << '\n';
                }
            }

            if (genomeLength % lineWrap != 0) {
                ref << '\n';
            }
        }

        return refPath.string();
    }

    // (2) Generate variants (VCF) with mason_variator (conservative default rates)
    auto generateVariantsVCF(const std::string& refPath, const std::string& tag) -> std::string {
        if (!mason_.variator) {
            return {};
        }

        auto vcfPath   = workDir_ / (tag + ".vcf");
        auto vcfString = vcfPath.string();

        if (fs::exists(vcfPath)) {
            return vcfString;
        }

        std::ostringstream cmd;
        // clang-format off
        cmd << "mason_variator"
            << " -ir " << std::quoted(refPath) 
            << " -ov " << std::quoted(vcfString)
            << " --snp-rate 0.001 --small-indel-rate 0.0001";
        // clang-format on
        runCommand(cmd.str(), "mason_variator failed");

        return vcfString;
    }

    // (3) Apply VCF to reference genome using mason_materializer
    auto materializeGenome(const std::string& refPath,
                           const std::string& vcfPath,
                           const std::string& tag) -> std::string {
        if (!mason_.materialize) {
            return {};
        }

        auto outFa       = workDir_ / (tag + "materialized.fa");
        auto outFaString = outFa.string();

        if (fs::exists(outFa)) {
            return outFaString;
        }

        std::ostringstream cmd;
        // clang-format off
        cmd << "mason_materializer"
            << " -ir " << std::quoted(refPath)
            << " -iv " << std::quoted(vcfPath)
            << " -o "  << std::quoted(outFaString);
        // clang-format on
        runCommand(cmd.str(), "mason_materializer failed: " + cmd.str());

        return outFaString;
    }

    // Helper function for single and pair end Illumina read generation
    auto generateCleanReadsImpl(const std::string&                refPath,
                                std::size_t                       numReads,
                                std::size_t                       readLength,
                                const std::string&                output1,
                                const std::optional<std::string>& output2 = std::nullopt) -> void {
        if (mason_.simulator) {
            // pick caller specified config or use safe defaults based on read length
            FragmentConfig fragmentConfig =
                fragment_cfg_.value_or(FragmentConfig::SafeDefaults(readLength));
            fragmentConfig.validate(readLength);

            std::ostringstream cmd;
            // clang-format off
            cmd << "mason_simulator"
                << " -ir " << std::quoted(refPath)
                << " -n "  << numReads
                << " --seq-technology illumina"
                << " --illumina-read-length " << readLength
                << " --num-threads " << cpu_count_ << ' ';
            // clang-format on

            cmd << " --fragment-size-model " << fragmentConfig.model;
            if (fragmentConfig.model == FragmentModel::Normal) {
                // clang-format off
                cmd << " --fragment-mean-size "     << fragmentConfig.mean
                    << " --fragment-size-std-dev "  << fragmentConfig.stddev;
                // clang-format on
            } else if (fragmentConfig.model == FragmentModel::Uniform) {
                // clang-format off
                cmd << " --fragment-min-size " << fragmentConfig.min_size
                    << " --fragment-max-size " << fragmentConfig.max_size;
                // clang-format on
            }

            cmd << " -o " << std::quoted(output1);
            if (output2) {
                cmd << " -or " << std::quoted(*output2);
            }

            std::string mode = output2 ? "PE" : "SE";
            runCommand(cmd.str(), "mason_simulator (" + mode + ") failed: " + cmd.str());
        } else {
            generateFallbackReads(output1, numReads, readLength);
            if (output2) {
                generateFallbackReads(*output2, numReads, readLength);
            }
        }
    }

    // (4-SE) Simulate single-end Illumina reads
    auto generateCleanReads(const std::string& refPath,
                            std::size_t        numReads,
                            std::size_t        readLength,
                            const std::string& tag) -> std::string {
        auto outFq = workDir_ / (tag + "_clean_SE.fastq");
        generateCleanReadsImpl(refPath, numReads, readLength, outFq.string());

        return outFq.string();
    }

    // (4-PE) Simulate paired-end Illumina reads
    auto generateCleanReadsPE(const std::string& refPath,
                              std::size_t        numReads,
                              std::size_t        readLength,
                              const std::string& tag) -> std::pair<std::string, std::string> {
        auto outFq1 = workDir_ / (tag + "_clean_R1.fastq");
        auto outFq2 = workDir_ / (tag + "_clean_R2.fastq");
        generateCleanReadsImpl(refPath, numReads, readLength, outFq1.string(), outFq2.string());

        return {outFq1.string(), outFq2.string()};
    }

    // 5 (SE) 3' adapter contamination (cut insert, append adapter prefix, fix quals)
    auto addAdapterContamination(const std::string& cleanFastq,
                                 double             contamRate,
                                 const std::string& suffix) -> std::string {
        constexpr int randomSeed = 666;

        auto outputPath = workDir_ / ("contaminated_" + suffix + ".fastq");
        auto input      = openFile<std::ifstream>(cleanFastq);
        auto output     = openFile<std::ofstream>(outputPath);

        contaminateStream(input, output, contamRate, randomSeed);
        return outputPath.string();
    }

    // (5-PE) contam both mates (this can be extended to correlated inserts)
    auto addAdapterContaminationPE(const std::string& cleanR1,
                                   const std::string& cleanR2,
                                   double             contamRate,
                                   const std::string& suffix)
        -> std::pair<std::string, std::string> {
        auto outR1 = workDir_ / ("contaminated_" + suffix + "_R1.fastq");
        auto outR2 = workDir_ / ("contaminated_" + suffix + "_R2.fastq");

        contamOneFile(cleanR1, outR1.string(), contamRate);
        contamOneFile(cleanR2, outR2.string(), contamRate);

        return {outR1.string(), outR2.string()};
    }

    // Fallback generator (used only if mason is missing)
    void generateFallbackReads(const std::string& outputPath,
                               std::size_t        numReads,
                               std::size_t        readLength) {
        constexpr int                 randomSeed = 69420;
        constexpr std::array<char, 4> bases      = {'A', 'T', 'G', 'C'};

        std::ofstream out(outputPath);

        std::mt19937                       rng(randomSeed);
        std::uniform_int_distribution<int> baseDist {0, 3};
        std::uniform_int_distribution<int> qualDist {33, 73};

        for (std::size_t i = 0; i < numReads; ++i) {
            std::string seq(readLength, 'A');
            for (auto& base : seq) {
                base = bases[baseDist(rng)];
            }

            std::string qual(readLength, 'I');
            for (auto& quality : qual) {
                quality = static_cast<char>(qualDist(rng));
            }

            // clang-format off
            out << "@fallback_read_" << i    << "\n"
                << seq << "\n+\n"    << qual << "\n";
            // clang-format on
        }
    }
};

// Generic 5 step pipeline (reference -> [variants] -> materialize -> simulate -> contam) -> FASTQ
template <bool IsPairedEnd>
auto MasonDataGenerator::generateDatasetImpl(const std::string& name,
                                             std::size_t        numReads,
                                             std::size_t        readLength,
                                             const std::string& platform,
                                             double             adapterContamRate,
                                             const std::string& errorModel,
                                             bool               withVariants,
                                             std::string_view   refFasta)
    -> DatasetReturn<IsPairedEnd> {
    // Check if output already exists
    if constexpr (IsPairedEnd) {
        const auto r1Out = workDir_ / (name + "_R1.fastq");
        const auto r2Out = workDir_ / (name + "_R2.fastq");
        if (fs::exists(r1Out) && fs::exists(r2Out)) {
            return std::make_pair(r1Out.string(), r2Out.string());
        }
    } else {
        const auto outputPath = workDir_ / (name + ".fastq");
        if (fs::exists(outputPath)) {
            return outputPath.string();
        }
    }

    // (1) reference: provided or synthetic
    const auto baseRef = (refFasta.empty() ? fs::path(generateReferenceGenome(name))
                                           : validateReferenceGenome(refFasta));

    // (2) optional variants (VCF) tagged by `name`
    std::string vcfPath;
    if (withVariants && mason_.variator) {
        vcfPath = generateVariantsVCF(baseRef, name);
    }

    // (3) materialize modified genome (if VCF present)
    const auto simRef = (withVariants && !vcfPath.empty() && mason_.materialize)
                            ? fs::path(materializeGenome(baseRef, vcfPath, name))
                            : baseRef;

    // Steps 4-5: branch based on single-end vs paired-end
    // Both add 3' adapter contamination
    if constexpr (IsPairedEnd) {
        auto [cleanR1, cleanR2] = generateCleanReadsPE(simRef, numReads, readLength, name);
        auto [contR1, contR2] =
            addAdapterContaminationPE(cleanR1, cleanR2, adapterContamRate, name);
        return std::make_pair(contR1, contR2);
    } else {
        const auto cleanReads        = generateCleanReads(simRef, numReads, readLength, name);
        const auto contaminatedReads = addAdapterContamination(cleanReads, adapterContamRate, name);
        return contaminatedReads;
    }
}

inline auto getGenerator() -> MasonDataGenerator& {
    static MasonDataGenerator instance;
    return instance;
}
}  // namespace benchmark_data