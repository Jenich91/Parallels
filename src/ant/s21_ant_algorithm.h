#ifndef sfleta_ANT_ALGORITHM
#define sfleta_ANT_ALGORITHM

#include <map>
#include <mutex>
#include <random>
#include <thread>

#include "sfleta_graph.h"

namespace sfleta {
static const constexpr size_t kInf = 1e8;
const double kAlfa = 1.0;
const double kBeta = 1.0;
const double kQ = 4.0;
const double kP = 0.5;
const double kInitPher = 0.2;
const size_t kColonies = 200;
const size_t kAnts = 500;

struct TsmResult {
  std::vector<int> vertices;
  double distance = INFINITY;
};

class AntAlgorithm {
 public:
  AntAlgorithm() = default;
  ~AntAlgorithm() = default;

  const TsmResult SolveTravelingSalesmanProblem(Graph &graph, bool mthread);
  void PrintResult(const TsmResult &);

 private:
  std::vector<std::vector<int>> adjacency_matrix_;
  std::vector<std::vector<double>> pheromones_matrix_;
  std::vector<TsmResult> ants_routes_;

  void FindRouteOfAnt(size_t graph_size);
  void ClearPheromonesMatrix();
  std::mutex mtx;

  void InitPheromonesMatrix(size_t graph_size);
  int ChoiceRandomNode(int iMin, int iMax);
  double ChoiceDoubleRand(double fMin, double fMax);
  std::map<int, double> FindNotVistedNodes(
      size_t graph_size, const std::vector<bool> &visted_nodes_,
      int current_node_);
  int ChoiceNode(std::map<int, double> &probabilitys);
  double CalculateLengthAntRoute(const std::vector<int> &ant_path_);
  void UpdatePheromonesMatrix();
  void UpdateResult(TsmResult &result);
};
}  // namespace sfleta
#endif  // sfleta_ANT_ALGORITHM