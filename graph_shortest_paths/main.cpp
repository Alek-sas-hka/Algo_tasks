#include <iostream>
#include <queue>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace graph {
template <typename T, typename WeightT>
struct DefaultEdge : std::pair<std::pair<T, T>, WeightT> {
  DefaultEdge(const T& first, const T& second, const WeightT& weight)
      : std::pair<std::pair<T, T>, WeightT>({first, second}, weight) {}
  using BaseClass = std::pair<std::pair<T, T>, WeightT>;
  const T& Start() const { return BaseClass::first.first; }
  const T& Finish() const { return BaseClass::first.second; }
  const WeightT& Weight() const { return BaseClass::second; }
};
template <typename Vertex = size_t, typename Weight = size_t,
          typename Edge = DefaultEdge<Vertex, Weight>>
class AbstractGraph {
 public:
  explicit AbstractGraph(size_t vertices_num, size_t edges_num = 0)
      : vertices_number_(vertices_num), edges_number_(edges_num) {}
  virtual std::vector<std::pair<Vertex, Weight>> GetNeighbours(
      const Vertex& vertex) = 0;
  const Vertex& GetVertexAmount() const { return vertices_number_; }
  const Vertex& GetEdgesAmount() const { return edges_number_; }

 protected:
  size_t vertices_number_ = 0;
  size_t edges_number_ = 0;
};

template <typename Vertex = size_t, typename Weight = size_t,
          typename Edge = DefaultEdge<Vertex, Weight>>
class AdjacencyListGraph : public AbstractGraph<Vertex, Weight, Edge> {
 public:
  AdjacencyListGraph(size_t vertices_num, const std::vector<Edge>& edges)
      : AbstractGraph<Vertex, Weight, Edge>(vertices_num, edges.size()) {
    for (const auto& edge : edges) {
      list[edge.Start()].push_back({edge.Finish(), edge.Weight()});
      list[edge.Finish()].push_back({edge.Start(), edge.Weight()});
    }
  }
  std::vector<std::pair<Vertex, Weight>> GetNeighbours(
      const Vertex& vertex) final {
    return list[vertex];
  }

  std::unordered_map<Vertex, std::vector<std::pair<Vertex, Weight>>> list;
};
}  // namespace graph

template <typename Vertex = size_t, typename Weight = size_t>
class Visitor {
 public:
  Visitor(size_t size, const Vertex& start, const Weight& max) {
    dist_ = std::vector<Weight>(size, max);
    prev_ = std::vector<Vertex>(size, start);
    point_ = start;
    dist_[start] = 0;
  }
  Weight GetDist(const Vertex& vert) { return dist_[vert]; }
  Vertex GetPoint() { return point_; }

  Vertex GetPrev(const Vertex& vert) { return prev_[vert]; }

  std::vector<Weight> GetAllDist() { return dist_; }

  void Relax(const Vertex& vert, const Vertex& vert2, const Weight& length) {
    dist_[vert2] = dist_[vert] + length;
    prev_[vert2] = vert;
  }

 private:
  Vertex point_;
  std::vector<Weight> dist_;
  std::vector<Vertex> prev_;
};

template <typename Vertex = size_t, typename Weight = size_t,
          typename Edge = graph::DefaultEdge<Vertex, Weight>>
void Dijkstra(Visitor<Vertex, Weight>& vis,
              graph::AbstractGraph<Vertex, Weight, Edge>& graph) {
  std::set<std::pair<Weight, Vertex>> queue;
  queue.insert({vis.GetDist(vis.GetPoint()), vis.GetPoint()});
  while (!queue.empty()) {
    Vertex nearest = queue.begin()->second;
    queue.erase(queue.begin());
    for (auto idx : graph.GetNeighbours(nearest)) {
      Vertex next = idx.first;
      if (vis.GetDist(nearest) + idx.second < vis.GetDist(next)) {
        queue.erase({vis.GetDist(next), next});
        vis.Relax(nearest, next, idx.second);
        queue.insert({vis.GetDist(next), next});
      }
    }
  }
}

int main() {
  size_t vertex_num;
  size_t edge_num;
  size_t map_num;
  std::cin >> map_num;
  const size_t kDoNotTry = 2009000999;
  for (size_t k = 0; k < map_num; ++k) {
    std::vector<graph::DefaultEdge<size_t, size_t>> edges;
    std::cin >> vertex_num >> edge_num;

    for (size_t i = 0; i < edge_num; ++i) {
      size_t edge_from;
      size_t edge_to;
      size_t edge_weight;
      std::cin >> edge_from >> edge_to >> edge_weight;
      edges.push_back({edge_from, edge_to, edge_weight});
    }
    graph::AdjacencyListGraph<size_t, size_t> graph(vertex_num, edges);

    size_t our_position;
    std::cin >> our_position;

    Visitor<size_t, size_t> vis(graph.GetVertexAmount(), our_position,
                                kDoNotTry);
    Dijkstra<size_t, size_t>(vis, graph);

    for (auto idx : vis.GetAllDist()) {
      std::cout << idx << " ";
    }
    std::cout << "\n";
  }
  return 0;
}
