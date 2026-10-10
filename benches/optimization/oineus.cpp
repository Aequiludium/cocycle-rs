// Native serial R/V critical-set adapter for the pinned Oineus API.
#include <oineus/params.h>
#include <oineus/decomposition.h>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <map>

using Clock = std::chrono::steady_clock;
using D = oineus::VRUDecomposition<int>;
using Targets = std::vector<std::pair<int, double>>;

struct Source {
  int n, q, repeats;
  std::vector<double> values;
  std::vector<int> dimensions;
  std::vector<std::vector<int>> boundaries;
  Targets proposals;
};

double elapsed(Clock::time_point start) {
  return std::chrono::duration<double, std::milli>(Clock::now() - start).count();
}

long memory(const std::string& name) {
  std::ifstream status("/proc/self/status");
  std::string line;
  while (std::getline(status, line)) {
    if (line.starts_with(name)) return std::stol(line.substr(name.size()));
  }
  throw std::runtime_error("missing process memory field");
}

class State {
  const Source& source;
  bool full_u;
  std::unique_ptr<D> primal, dual;
  std::map<int, D::MatrixData> primal_rows, dual_rows;

  std::unique_ptr<D> build(bool dualize) {
    D::MatrixData matrix(source.n);
    for (int j = 0; j < source.n; ++j) {
      for (int i : source.boundaries[j]) {
        if (dualize) matrix[source.n - 1 - i].push_back(source.n - 1 - j);
        else matrix[j].push_back(i);
      }
    }
    for (auto& column : matrix) std::sort(column.begin(), column.end());
    auto state = std::make_unique<D>(matrix);
    oineus::ReductionParams params;
    params.compute_v = true;
    params.use_clearing = false;
    state->reduce(params);
    if (full_u) {
      // Raw-square constructors expose one matrix block; this deliberately
      // computes the full inverse in both decompositions for the ablation.
      auto value_at = [&](int i) { return source.values[dualize ? source.n - 1 - i : i]; };
      state->compute_full_u_rows<double>(0, value_at, 1, false);
    }
    return state;
  }

public:
  State(const Source& input, bool full) : source(input), full_u(full), primal(build(false)) {}

  std::vector<std::pair<int, int>> pairs() const {
    std::map<int, int> paired;
    for (int j = 0; j < source.n; ++j)
      if (!primal->r_data[j].empty()) paired.emplace(primal->r_data[j].back(), j);
    return {paired.begin(), paired.end()};
  }

  std::vector<int> essential() const {
    std::vector<bool> paired(source.n);
    for (auto [b, d] : pairs()) paired[b] = true;
    std::vector<int> result;
    for (int i = 0; i < source.n; ++i)
      if (primal->r_data[i].empty() && !paired[i]) result.push_back(i);
    return result;
  }

  Targets targets() {
    std::map<int, double> result;
    for (auto [id, target] : source.proposals) {
      if (target == source.values[id]) continue;
      bool negative = !primal->r_data[id].empty();
      bool increase = target > source.values[id];
      bool dualize = !negative; // Benchmark inputs contain finite endpoints only.
      bool use_u = negative == increase;
      if (dualize && !dual) dual = build(true);
      auto& state = dualize ? *dual : *primal;
      int index = dualize ? source.n - 1 - id : id;
      D::IntSparseColumn support;
      if (use_u && full_u) {
        support = state.u_data_t[index];
      } else if (use_u) {
        auto& cache = dualize ? dual_rows : primal_rows;
        int dimension = source.dimensions[id];
        if (!cache.contains(dimension)) {
          D::MatrixData rows(source.n);
          for (int j = 0; j < source.n; ++j) {
            int original = dualize ? source.n - 1 - j : j;
            if (source.dimensions[original] == dimension)
              for (int i : state.v_data[j]) rows[i].push_back(j);
          }
          cache.emplace(dimension, std::move(rows));
        }
        auto value_at = [&](int i) { return dualize ? -source.values[source.n - 1 - i] : source.values[i]; };
        support = state.compute_u_row_bounded(index, cache.at(dimension), dualize ? -target : target,
            value_at, [](double value, double bound) { return value > bound; });
      } else {
        support = state.v_data[index];
      }
      for (int i : support) {
        int original = dualize ? source.n - 1 - i : i;
        double value = source.values[original];
        if (increase ? value > target : value < target) continue;
        auto previous = result.find(original);
        if (previous == result.end()) result.emplace(original, target);
        else if (std::abs(target - value) > std::abs(previous->second - value)) previous->second = target;
      }
    }
    return {result.begin(), result.end()};
  }
};

int main(int argc, char** argv) {
  if (argc != 3) throw std::runtime_error("usage: oineus MODE FIXTURE");
  std::ifstream input(argv[2]);
  Source source;
  int queries; double scale;
  input >> source.n >> source.q >> queries >> source.repeats >> scale;
  source.values.resize(source.n); source.dimensions.resize(source.n); source.boundaries.resize(source.n);
  std::map<std::vector<int>, int> ids;
  for (int j = 0; j < source.n; ++j) {
    int k; input >> k;
    std::vector<int> vertices(k);
    for (int& v : vertices) input >> v;
    input >> source.values[j]; source.dimensions[j] = k - 1;
    if (k > 1) for (int omitted = 0; omitted < k; ++omitted) {
      auto face = vertices; face.erase(face.begin() + omitted);
      source.boundaries[j].push_back(ids.at(face));
    }
    ids.emplace(vertices, j);
    std::sort(source.boundaries[j].begin(), source.boundaries[j].end());
  }
  for (int i = 0; i < queries; ++i) { int id; double target; input >> id >> target; source.proposals.emplace_back(id, target); }
  bool full = std::string(argv[1]) == "oineus_full";
  long rss = memory("VmRSS:"), hwm = memory("VmHWM:");
  auto begin = Clock::now();
  auto state = std::make_unique<State>(source, full);
  double build_ms = elapsed(begin);
  auto query = Clock::now(); auto targets = state->targets(); double cold_ms = elapsed(query);
  auto pairs = state->pairs(); auto essential = state->essential();
  auto cleanup = Clock::now(); state.reset(); double cleanup_ms = elapsed(cleanup);
  state = std::make_unique<State>(source, full);
  if (state->targets() != targets) throw std::runtime_error("unstable targets");
  auto warm = Clock::now();
  size_t checksum = 0;
  for (int i = 0; i < source.repeats; ++i) {
    auto output = state->targets();
    checksum += output.size();
    asm volatile("" : : "g"(output.data()) : "memory");
  }
  double warm_ms = elapsed(warm) / source.repeats;
  std::cout << std::setprecision(17) << "{\"protocol\":\"critical-sets-native-v1\",\"times_ms\":{\"build\":" << build_ms
      << ",\"cold_query\":" << cold_ms << ",\"cleanup\":" << cleanup_ms << ",\"cold_workflow\":" << build_ms + cold_ms + cleanup_ms
      << ",\"warm_batch\":" << warm_ms << "},\"rss_input_kib\":" << rss << ",\"hwm_input_kib\":" << hwm
      << ",\"peak_rss_kib\":" << memory("VmHWM:") << ",\"extra_payload_count\":0,\"targets\":[";
  for (size_t i = 0; i < targets.size(); ++i) { if (i) std::cout << ','; std::cout << '[' << targets[i].first << ',' << targets[i].second << ']'; }
  std::cout << "],\"pairs\":[";
  for (size_t i = 0; i < pairs.size(); ++i) { if (i) std::cout << ','; std::cout << '[' << pairs[i].first << ',' << pairs[i].second << ']'; }
  std::cout << "],\"essential\":[";
  for (size_t i = 0; i < essential.size(); ++i) { if (i) std::cout << ','; std::cout << essential[i]; }
  std::cout << "],\"intervals\":[";
  bool first = true;
  for (auto [b,d] : pairs) if (source.dimensions[b] <= source.q && source.values[b] < source.values[d]) {
    if (!first) std::cout << ','; first = false;
    std::cout << '[' << source.dimensions[b] << ',' << source.values[b] << ',' << source.values[d] << ']';
  }
  for (int b : essential) if (source.dimensions[b] <= source.q) {
    if (!first) std::cout << ','; first = false;
    std::cout << '[' << source.dimensions[b] << ',' << source.values[b] << ",null]";
  }
  std::cout << "],\"checksum\":" << checksum << "}\n";
}
