#ifndef sfleta_GRAPH_H
#define sfleta_GRAPH_H

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace sfleta {
class Graph {
 public:
  Graph() : adjacency_matrix_{}, size_{0} {}
  void LoadGraphFromFile(std::string);

  auto GetAdjacencyMatrix() { return adjacency_matrix_; }
  auto GetSize() { return size_; }

  void Print();

 private:
  std::vector<std::vector<int>> adjacency_matrix_;
  size_t size_;

  void Clear();
};
}  // namespace sfleta
#endif  // sfleta_GRAPH_H
