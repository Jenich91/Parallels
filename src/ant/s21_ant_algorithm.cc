#include "sfleta_ant_algorithm.h"

namespace sfleta {
const TsmResult AntAlgorithm::SolveTravelingSalesmanProblem(Graph &graph,
                                                            bool mthread) {
  TsmResult result;
  adjacency_matrix_ = graph.GetAdjacencyMatrix();
  size_t graph_size = graph.GetSize();

  InitPheromonesMatrix(graph_size);

  for (size_t step_i = 0; step_i < kColonies; step_i++) {
    for (size_t ant_i = 0; ant_i < kAnts; ant_i++) {
      if (mthread == false) {
        FindRouteOfAnt(graph_size);
      } else {
        int count_thread = std::thread::hardware_concurrency();
        std::vector<std::thread> vector_thread;
        int tmp = 0;
        for (; tmp < count_thread; ++tmp, ++ant_i) {
          vector_thread.push_back(
              std::thread(&AntAlgorithm::FindRouteOfAnt, this, graph_size));
        }
        for (int k = 0; k < tmp; ++k) {
          vector_thread[k].join();
        }
      }
    }
    UpdatePheromonesMatrix();
    if (ants_routes_.size() != 0) {
      UpdateResult(result);
      ants_routes_.clear();
    }
  }

  if (result.distance == INFINITY) {
    throw std::logic_error("Can't solve Salesman Problem\n");
  }
  return result;
}

void AntAlgorithm::FindRouteOfAnt(size_t graph_size) {
  bool loop = false;
  std::vector<int> ant_path_;
  std::vector<bool> visted_nodes_;
  int current_node_;
  double l_0 = INFINITY;

  visted_nodes_.resize(graph_size, true);
  current_node_ = ChoiceRandomNode(0, graph_size - 1);

  ant_path_.push_back(current_node_);
  visted_nodes_[current_node_] = false;

  for (size_t node_i = 0; node_i < graph_size - 1; node_i++) {
    std::map<int, double> not_visited_nodes =
        FindNotVistedNodes(graph_size, visted_nodes_, current_node_);

    std::map<int, double> probabilitys;
    double sum_p_all = 0;
    for (const auto &m : not_visited_nodes) {
      sum_p_all += pow(pheromones_matrix_[current_node_][m.first], kAlfa) *
                   pow(1 / m.second, kBeta);
    }

    double range = 0;
    for (const auto &j : not_visited_nodes) {
      double p_ij = 0;
      double value_teta =
          pow(pheromones_matrix_[current_node_][j.first], kAlfa);
      double value_eta = pow(1 / j.second, kBeta);
      p_ij = (value_teta * value_eta) / sum_p_all;
      range += p_ij;
      probabilitys.insert(std::make_pair(j.first, range));
    }
    int next_node = ChoiceNode(probabilitys);

    if (next_node != -1) {
      if (adjacency_matrix_[current_node_][next_node] != 0) {
        ant_path_.push_back(next_node);
        visted_nodes_[next_node] = false;
        current_node_ = next_node;
      }
    } else {
      loop = true;
      break;
    }
  }
  if (loop != true) {
    if (adjacency_matrix_[current_node_][ant_path_.front()] != 0) {
      ant_path_.push_back(ant_path_.front());
      if (ant_path_.size() == adjacency_matrix_.size() + 1) {
        double current_len_route = CalculateLengthAntRoute(ant_path_);
        if (current_len_route < l_0) {
          l_0 = current_len_route;
          TsmResult tmp_res{{ant_path_}, {l_0}};
          mtx.lock();
          ants_routes_.push_back(tmp_res);
          mtx.unlock();
        }
      }
    }
  }
}

void AntAlgorithm::ClearPheromonesMatrix() {
  if (pheromones_matrix_.size()) {
    for (auto &row : pheromones_matrix_) row.clear();
    pheromones_matrix_.clear();
  }
}

void AntAlgorithm::InitPheromonesMatrix(size_t graph_size) {
  ClearPheromonesMatrix();
  pheromones_matrix_.resize(graph_size, std::vector<double>(graph_size));
  for (size_t i = 0; i < graph_size; i++) {
    for (size_t j = 0; j < graph_size; j++) {
      if (adjacency_matrix_[i][j] != 0) {
        pheromones_matrix_[i][j] = kInitPher;
      }
    }
  }
}

int AntAlgorithm::ChoiceRandomNode(int iMin, int iMax) {
  std::random_device rd;
  std::default_random_engine engine(rd());
  std::uniform_int_distribution<int> node(iMin, iMax);
  return node(engine);
}

double AntAlgorithm::ChoiceDoubleRand(double fMin, double fMax) {
  std::random_device rd;
  std::default_random_engine eng(rd());
  std::uniform_real_distribution<double> random_double(fMin, fMax);
  return random_double(eng);
}

std::map<int, double> AntAlgorithm::FindNotVistedNodes(
    size_t graph_size, const std::vector<bool> &visted_nodes_,
    int current_node_) {
  std::map<int, double> list_not_visited_nodes;

  for (size_t j = 0; j < graph_size; j++) {
    if (visted_nodes_[j] && adjacency_matrix_[current_node_][j] != 0) {
      list_not_visited_nodes.insert(
          std::make_pair(j, adjacency_matrix_[current_node_][j]));
    }
  }
  return list_not_visited_nodes;
}

int AntAlgorithm::ChoiceNode(std::map<int, double> &probabilitys) {
  int next_node;
  if (probabilitys.size() == 1) {
    next_node = probabilitys.begin()->first;
  } else {
    if (probabilitys.size() != 0) {
      double random = ChoiceDoubleRand(0, 1);
      for (auto it = probabilitys.begin(); it != probabilitys.end(); ++it) {
        if (random <= it->second) {
          next_node = it->first;
          break;
        }
      }
    } else {
      next_node = -1;
    }
  }
  return next_node;
}

double AntAlgorithm::CalculateLengthAntRoute(
    const std::vector<int> &ant_path_) {
  double length_route = 0;
  for (size_t i = 0; i < ant_path_.size() - 1; i++) {
    length_route += adjacency_matrix_[ant_path_[i]][ant_path_[i + 1]];
  }
  return length_route;
}

void AntAlgorithm::UpdatePheromonesMatrix() {
  for (size_t i = 0; i < pheromones_matrix_.size(); i++) {
    for (size_t j = 0; j < pheromones_matrix_[i].size(); j++) {
      pheromones_matrix_[i][j] = pheromones_matrix_[i][j] * (1 - kP);
    }
  }

  for (size_t i = 0; i < ants_routes_.size(); i++) {
    double delta_teta = kQ / ants_routes_[i].distance;
    for (size_t j = 0; j < ants_routes_[i].vertices.size() - 1; j++) {
      pheromones_matrix_[ants_routes_[i].vertices[j]]
                        [ants_routes_[i].vertices[j + 1]] += delta_teta;
    }
  }
}

void AntAlgorithm::UpdateResult(TsmResult &result) {
  for (size_t i = 0; i < ants_routes_.size(); i++) {
    if (ants_routes_[i].distance < result.distance) {
      result.vertices = ants_routes_[i].vertices;
      result.distance = ants_routes_[i].distance;
    }
  }
}

void AntAlgorithm::PrintResult(const TsmResult &result) {
  auto [vertices, distance] = (result);
  std::cout << "Выгодный (короткий) маршрут:" << std::endl;
  for (auto &value : vertices) std::cout << value + 1 << ' ';
  std::cout << "\nДлина этого маршрута: " << distance << std::endl;
}

}  // namespace sfleta
