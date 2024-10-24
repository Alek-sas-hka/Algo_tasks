#include <functional>
#include <iostream>
#include <limits>
#include <unordered_map>
#include <unordered_set>
#include <vector>

const int64_t kDoNotTry = std::numeric_limits<int>::max();

template <typename Vertex>
struct DefaultEdge : std::pair<Vertex, Vertex> {
  using BaseClass = std::pair<Vertex, Vertex>;

  DefaultEdge(const Vertex& first, const Vertex& second)
      : BaseClass(first, second) {}

  const Vertex& GetFrom() const { return BaseClass::first; }
  const Vertex& GetTo() const { return BaseClass::second; }
};


namespace std {
template <typename Vertex>
struct hash<DefaultEdge<Vertex>> {
  size_t operator()(const DefaultEdge<Vertex>& edge) const {
    return edge.GetFrom() ^ edge.GetTo();
  }
};
} 
template <typename Vertex>
DefaultEdge<Vertex> GetReverseEdge(const DefaultEdge<Vertex>& edge) {
  return {edge.GetTo(), edge.GetFrom()};
}

template <typename VertexType = int,
          typename EdgeType = DefaultEdge<VertexType>>
class Graph {
 public:
  using Vertex = VertexType;
  using Edge = EdgeType;
  using InternalEdgeContainer = std::unordered_set<EdgeType>;

  explicit Graph(size_t number_vertices) : number_vertices_(number_vertices) {}

  void AddEdge(const EdgeType& edge) {
    adjacency_list_[edge.GetFrom()].insert({edge.GetFrom(), edge.GetTo()});
    ++number_edges_;
  }

  void AddVertex() {
    adjacency_list_.emplace_back();
    ++number_vertices_;
  }

  size_t VerticesCount() const { return number_vertices_; }

  size_t EdgesCount() const { return number_edges_; }

  InternalEdgeContainer& EdgesFromVertex(int vertex) {
    return adjacency_list_[vertex];
  }

  size_t EdgesNumberFromVertex(int vertex) const {
    return adjacency_list_[vertex].size();
  }

 private:
  size_t number_edges_ = 0;
  size_t number_vertices_ = 0;
  std::unordered_map<Vertex, std::unordered_set<Edge>> adjacency_list_;
};

template <class Graph>
class FilteredGraph {
 public:

  FilteredGraph(Graph& graph,
                std::function<bool(typename Graph::Edge)> predicate)
      : graph_(graph), predicate_(std::move(predicate)) {}

  typename Graph::InternalEdgeContainer EdgesFromVertex(
      typename Graph::Vertex& vertex) {
    auto outgoing_edges = graph_.EdgesFromVertex(vertex);
    typename Graph::InternalEdgeContainer filtered_edges;
    for (const auto& edge : outgoing_edges) {
      if (predicate_(edge)) {
        filtered_edges.insert(edge);
      }
    }
    return filtered_edges;
  }

  size_t VerticesCount() const { return graph_.VerticesCount(); }

 private:
  Graph& graph_;
  std::function<bool(typename Graph::Edge)> predicate_;
};

template <typename FlowNetwork>
class GoldNetworkBuilder;

template <typename GraphType>
class FlowNetwork {
 public:
  using Vertex = typename GraphType::Vertex;
  using Edge = typename GraphType::Edge;
  using Graph = GraphType;

  explicit FlowNetwork(Graph graph) : graph_(std::move(graph)) {}

  void UpdateFlowByEdge(const Edge& edge, int extra_flow) {
    edges_properties_[edge].flow += extra_flow;
    const auto kReverseEdge = GetReverseEdge<Vertex>(edge);
    edges_properties_[kReverseEdge].flow -= extra_flow;
  }

  int ResidualCapacity(const Edge& edge) {
    return edges_properties_[edge].capacity - edges_properties_[edge].flow;
  }

  size_t NumberEdges() const { return graph_.EdgesCount(); }

  size_t NumberVertices() const { return graph_.VerticesCount(); }

  Vertex& Source() { return source_; }

  Vertex& Sink() { return sink_; }

  FilteredGraph<GraphType> ResidualNetworkView() {
    return FilteredGraph(graph_, [this](const Edge& edge) {
      return ResidualCapacity(edge) > 0;
    });
  }

  friend class GoldNetworkBuilder<FlowNetwork>;

 private:
  struct EdgeProperties {
    EdgeProperties() = default;
    EdgeProperties(int flow, int capacity, size_t reverse_edge_id)
        : flow(flow), capacity(capacity), reverse_edge_id(reverse_edge_id) {}
    int flow = 0;
    int capacity = 0;
    size_t reverse_edge_id = 0;
  };

  Vertex source_;
  Vertex sink_;
  std::unordered_map<Edge, EdgeProperties> edges_properties_;
  GraphType graph_;
};

template <class FlowNetwork>
class GoldNetworkBuilder {
 public:
  using Vertex = typename FlowNetwork::Vertex;
  using Edge = typename FlowNetwork::Edge;

  explicit GoldNetworkBuilder(size_t vertices_number)
      : network_(Graph(vertices_number)) {}
  void AddEdge(const Edge& edge, int capacity) {
    network_.edges_properties_[edge] = {0, capacity,
                                        network_.edges_properties_.size() + 1};
    network_.graph_.AddEdge(edge);

    const auto kReverseEdge = GetReverseEdge<Vertex>(edge);
    network_.edges_properties_[kReverseEdge] = {
        0, 0, network_.edges_properties_.size() - 1};
    network_.graph_.AddEdge(kReverseEdge);
  }

  void AssignSource(const Vertex& source) { network_.Source() = source; }

  void AssignSink(const Vertex& sink) { network_.Sink() = sink; }

  FlowNetwork GetFlowNetwork() const { return std::move(network_); }

 private:
  FlowNetwork network_;
};

template <typename Vertex = int>
class Visitor {
 public:
  Visitor() {
    res_capacity_ = kDoNotTry;
    prev_res_cap_ = res_capacity_;
  }

  bool NotBeenHere(Vertex point) {
    return (visited_.find(point) == visited_.end());
  }

  void VisitVertex(Vertex point) { visited_.insert(point); }

  void SetResCapacity(int cap) {
    prev_res_cap_ = res_capacity_;
    res_capacity_ = cap;
  }

  int ResCapacity() { return res_capacity_; }

  int PrevCap() { return prev_res_cap_; }

 private:
  int res_capacity_;
  int prev_res_cap_;
  std::unordered_set<Vertex> visited_;
};

template <class Vertex = int>
FlowNetwork<Graph<Vertex>> Creation(
    int country_b, GoldNetworkBuilder<FlowNetwork<Graph<Vertex>>>& builder) {
  for (int i = 0; i < country_b; i++) {
    Vertex from;
    Vertex dst;
    int capacity;
    std::cin >> from >> dst >> capacity;
    builder.AddEdge({from - 1, dst - 1}, capacity);
  }
  return builder.GetFlowNetwork();
}

template <class Network, class Vertex>
void Dfs(Vertex start, Vertex target, Network& network, Visitor<Vertex>& vis) {
  vis.VisitVertex(start);
  if (start == target) {
    return;
  }

  for (auto tmp : network.ResidualNetworkView().EdgesFromVertex(start)) {
    if (vis.NotBeenHere(tmp.second)) {
      vis.SetResCapacity(
          std::min(vis.ResCapacity(),
                   network.ResidualCapacity({tmp.first, tmp.second})));
      int save = vis.ResCapacity();
      Dfs(tmp.second, target, network, vis);
      // vis.SetResCapacity(t);
      if (vis.ResCapacity() <= 0) {
        vis.SetResCapacity(save);
        continue;
      }
      network.UpdateFlowByEdge({start, tmp.second}, vis.ResCapacity());
      return;
    }
  }
  vis.SetResCapacity(0);
}

template <class Network>
int MaxFlow(Network& network) {
  int flow = 0;

  while (true) {
    Visitor<int> vis;
    Dfs(network.Source(), network.Sink(), network, vis);
    int res = vis.ResCapacity();
    if (res <= 0) {
      break;
    }
    flow += res;
  }
  return flow;
}

int main() {
  int country_a;
  int country_b;
  std::cin >> country_a >> country_b;
  GoldNetworkBuilder<FlowNetwork<Graph<int>>> builder(country_a - 1);
  builder.AssignSink(country_a - 1);
  builder.AssignSource(0);

  auto network = Creation(country_b, builder);

  std::cout << MaxFlow(network) << "\n";

  return 0;
}
