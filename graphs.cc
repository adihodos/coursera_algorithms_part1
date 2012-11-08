#include <cassert>
#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <cctype>

#include <algorithm>
#include <functional>
#include <iterator>
#include <list>
#include <memory>
#include <queue>
#include <stack>
#include <vector>
#include <unordered_map>

#define LAST_ERRNO(func_name) \
    do { \
        auto err_code = errno; \
        fprintf(stderr, "\n[%s, %d] function %s failed, errno %d", \
                __FILE__, __LINE__, #func_name, err_code); \
    } while (0)

template<typename T, size_t U>
inline size_t constexpr array_size(T (&)[U]) noexcept {
    return U;
}

struct crt_file_deleter {
    void operator()(FILE* fp) const noexcept {
        if (fp)
            fclose(fp);
    }
};

struct graph_arc {
    int32_t destination;
    int32_t weight;
};

struct graph_edge {
    int32_t source;
    int32_t destination;
    int32_t weight;
};

enum class arc_direction : char {
    incoming = 0,
    outgoing
};


struct graph_vertex {
    int32_t id;
    std::list<graph_arc> edges[2];

    typedef std::list<graph_arc> arc_list;
    typedef arc_list::iterator arc_list_iterator;
    typedef arc_list::const_iterator arc_list_const_iterator;
    typedef arc_list::reverse_iterator arc_list_reverse_iterator;
    typedef arc_list::const_reverse_iterator arc_list_const_reverse_iterator;

    arc_list_const_iterator cbegin() const {
        return edges[(int) arc_direction::outgoing].cbegin();
    }

    arc_list_const_iterator cend() const {
        return edges[(int) arc_direction::outgoing].cend();
    }

    arc_list_const_iterator crbegin() const {
        return edges[(int) arc_direction::incoming].cbegin();
    }

    arc_list_const_iterator crend() const {
        return edges[(int) arc_direction::incoming].cend();
    }
};

class graph {
private :
    typedef std::unordered_map<int32_t, graph_vertex>   vertex_table_t;

    vertex_table_t  vertices_;
    bool            directed_;

    void add_vertex(int32_t source, int32_t dest);

    void add_arc(int32_t source, int32_t dest, arc_direction type);

    inline bool check_valid_vertex_id(int32_t id) const noexcept {
        return id > 0 && id < (int32_t) vertices_.size();
    }

public :
    typedef vertex_table_t::iterator            vertex_iterator;
    typedef vertex_table_t::const_iterator      vertex_const_iterator;

    graph(bool directed = false);

    bool is_directed() const {
        return directed_;
    }

    bool read_from_file(const char* file_name);

    size_t inner_degree(int32_t vertex) noexcept {
        assert(check_valid_vertex_id(vertex));
        return vertices_[vertex].edges[(int) arc_direction::incoming].size();
    }

    size_t outer_degree(int32_t vertex) noexcept {
        assert(check_valid_vertex_id(vertex));
        return vertices_[vertex].edges[(int) arc_direction::outgoing].size();
    }

    size_t get_vertex_count() const noexcept {
        return vertices_.size() - 1;
    }

    const graph_vertex& get_vertex(int32_t id) const {
        assert(check_valid_vertex_id(id));
        return vertices_.find(id)->second;
    }

    void debug_print() const noexcept;

    vertex_iterator begin() {
        return vertices_.begin();
    }

    vertex_iterator end() {
        return vertices_.end();
    }

    vertex_const_iterator cbegin() {
        return vertices_.cbegin();
    }

    vertex_const_iterator cend() {
        return vertices_.cend();
    }
};

graph::graph(bool directed)
    :   vertices_(),
        directed_(directed) {
    vertices_[0] = graph_vertex();
}

bool graph::read_from_file(const char *file_name) {
    std::unique_ptr<FILE, crt_file_deleter> fp(fopen(file_name, "r"));
    if (!fp) {
        LAST_ERRNO(fopen());
        return false;
    }

    vertices_.clear();
    vertices_[0] = graph_vertex();

    char buff_line[1024];
    while (fgets(buff_line, array_size(buff_line), fp.get())) {
        const char* p = buff_line;
        int32_t source_vertex = 0;

        for (;;) {
            while (*p && isspace(*p))
                ++p;

            if (!*p)
                break;

            char tmp_buff[64];
            int i = 0;
            while (*p && isdigit(*p)) {
                tmp_buff[i++] = *p++;
            }

            if (!i)
                break;

            tmp_buff[i] = '\0';
            int32_t new_vertex = atoi(tmp_buff);
            if (!new_vertex)
                break;

            if (!source_vertex) {
                source_vertex = new_vertex;
                add_vertex(source_vertex, -1);
            } else {
                //
                // add vertex
                add_vertex(source_vertex, new_vertex);
            }
        }
    }

    return true;
}

void graph::add_vertex(int32_t source, int32_t dest) {
    if (dest == -1) {
        graph_vertex& gvert = vertices_[source];
        gvert.id = source;
        return;
    }

    add_arc(source, dest, arc_direction::outgoing);
    add_arc(dest, source, arc_direction::incoming);

    if (!directed_) {
        add_arc(source, dest, arc_direction::incoming);
        add_arc(dest, source, arc_direction::outgoing);
    }
}

void graph::add_arc(int32_t source, int32_t dest, arc_direction type) {
    //assert(check_valid_vertex_id(source));
    //assert(check_valid_vertex_id(dest));

    graph_vertex& vertex = vertices_[source];
    vertex.id = source;
    vertex.edges[(int) type].push_back({ dest, -1 });
}

void graph::debug_print() const noexcept {
    auto itr_vertex = vertices_.cbegin();
    auto itr_vertex_end = vertices_.cend();

    for (++itr_vertex; itr_vertex != itr_vertex_end; ++itr_vertex) {
        printf("\nVertex %d", itr_vertex->first);
        const char* const kArcDescription[] = { "Incoming", "Outgoing" };
        
        for (int j = 0; j < (int) array_size(itr_vertex->second.edges); ++j) {
            printf("\n%s edges : \n", kArcDescription[j]);

            auto itr_edge = itr_vertex->second.edges[j].cbegin();
            auto itr_end = itr_vertex->second.edges[j].cend();

            while (itr_edge != itr_end) {
                printf(" (%d -> %d) ", itr_vertex->first, itr_edge->destination);
                ++itr_edge;
            }
        }
    }
}

enum class node_color : char {
    white = 0,
    gray,
    black
};

class graph_search_base {
public :
    typedef std::vector<int32_t>    distance_table_t;
    typedef std::vector<int32_t>    predecessor_table_t;
    
    const distance_table_t& get_distance_table() const {
        return distances_;
    }

    const predecessor_table_t& get_predecessors_table() const {
        return predecessors_;
    }

protected :
    graph&                  graph_;
    std::vector<node_color> colors_;
    distance_table_t        distances_;
    predecessor_table_t     predecessors_;
    
    void initialize() {
        const size_t vertex_count = graph_.get_vertex_count() + 1;
        colors_.assign(vertex_count, node_color::white);
        distances_.assign(vertex_count, -1);
        predecessors_.assign(vertex_count, -1);
    }

    graph_search_base(graph& g) : graph_(g) {}

    ~graph_search_base() {}

    int32_t get_distance(int32_t vertex_id) const {
        return distances_[vertex_id];
    }

    int32_t get_predecessor(int32_t vertex_id) const {
        return predecessors_[vertex_id];
    }

    node_color get_color(int32_t vertex_id) const {
        return colors_[vertex_id];
    }
};

class bfs_search {
public :
    typedef std::vector<int32_t>    distance_table_t;
    typedef std::vector<int32_t>    predecessor_table_t;
private :
    graph&                  graph_;
    std::vector<node_color> color_table_;
    distance_table_t        distances_;
    predecessor_table_t     predecessors_;

    void initialize() {
        const size_t vertex_count = graph_.get_vertex_count() + 1;
        color_table_.assign(vertex_count, node_color::white);
        distances_.assign(vertex_count, -1);
        predecessors_.assign(vertex_count, -1);
    }

public :
    bfs_search(graph& g)
        : graph_(g) {
    }

    void execute(int32_t start_vertex, bool reverse = false);

    const distance_table_t& get_distance_table() const {
        return distances_;
    }

    const predecessor_table_t& get_predecessors_table() const {
        return predecessors_;
    }
};

void bfs_search::execute(int32_t start_vertex, bool reverse) {
    initialize();

    std::queue<int32_t> vertex_queue;
    vertex_queue.push(start_vertex);
    color_table_[start_vertex] = node_color::gray;
    distances_[start_vertex] = 0;

    typedef graph_vertex::arc_list_const_iterator (graph_vertex::*const_itr_fn)() const;
    const_itr_fn itr_fn[4] = { 
        &graph_vertex::cbegin, &graph_vertex::cend,
        &graph_vertex::crbegin, &graph_vertex::crend,
    };

    while (!vertex_queue.empty()) {
        auto src_vertex = graph_.get_vertex(vertex_queue.front());

        auto adj_itr = (src_vertex.*itr_fn[reverse * 2])();
        auto adj_end = (src_vertex.*itr_fn[reverse * 2 + 1])();

        while (adj_itr != adj_end) {
            printf("\n(%d ~> %d)", src_vertex.id, adj_itr->destination);

            if (color_table_[adj_itr->destination] == node_color::white) {
                color_table_[adj_itr->destination] = node_color::gray;
                predecessors_[adj_itr->destination] = src_vertex.id;
                distances_[adj_itr->destination] = distances_[src_vertex.id] + 1;
                vertex_queue.push(adj_itr->destination);
            }
            ++adj_itr;
        }

        color_table_[src_vertex.id] = node_color::black;
        vertex_queue.pop();
    }
}

struct null_edge_discovery_fn {
    void operator()(int32_t, int32_t) const {}
};

struct null_vertex_fn {
    void operator()(int32_t) const {}
};

//! Implements DFS search in a graph.
class dfs_search : private graph_search_base {
public :
    //! \name Types
    //! @{

    using graph_search_base::distance_table_t;

    using graph_search_base::predecessor_table_t;

    //! @}

public :
    //! \name Constructors
    //! @{

    dfs_search(graph& g) : graph_search_base(g) {}

    //! @}
    
public :
    //! \name Operations
    //! @{

    //! \brief Executes a DFS in the graph.
    //! \param start_vertex Vertex where the search is initiated.
    //! \param reversed If true, search is executed on the transposed graph.
    //! \param edge_discovery_fn Functor called when examining an edge.
    //! \param vertex_discover_fn Functor called when a vertex is explored for the
    //!  first time.
    //! \param vertex_end_fn Functor called when all nodes in a vertex's adjacency
    //! list have been explored.
    template<
        typename edge_discovery_fn,
        typename vertex_discover_fn,
        typename vertex_end_fn
    >
    void execute(
        int32_t start_vertex, 
        bool reversed = false,
        edge_discovery_fn edge_fn = null_edge_discovery_fn(),
        vertex_discover_fn vx_discover_fn = null_vertex_fn(),
        vertex_end_fn vx_end_fn = null_vertex_fn()
        );

    using graph_search_base::initialize;

    //! @}

public :
    //! \name Attributes
    //! @{
    
    using graph_search_base::get_distance;

    using graph_search_base::get_predecessor;

    using graph_search_base::get_color;

    //! @}
};

template<
    typename edge_fn, 
    typename vertex_discover_fn, 
    typename vertex_end_fn
> 
void dfs_search::execute(
    int32_t start_vertex, 
    bool reversed,
    edge_fn edge_explore_fn,
    vertex_discover_fn vx_explore_fn,
    vertex_end_fn vx_end_fn
    ) {

    //
    // Saves search state for a given vertex.
    struct dfs_frame_t {
        const graph_vertex*                     vert;
        graph_vertex::arc_list_const_iterator   itr_adj;
        graph_vertex::arc_list_const_iterator   itr_adj_end;
    };
    
    //
    // Pointer to functions that return iterators, suitable for
    // traversing the graph in forward/reversed direction.
    typedef graph_vertex::arc_list_const_iterator 
        (graph_vertex::*const_itr_fn)() const;

    const_itr_fn itr_fn[4] = { 
        &graph_vertex::cbegin, &graph_vertex::cend,
        &graph_vertex::crbegin, &graph_vertex::crend,
    };
    
    //
    // Initialize search data.
    std::stack<dfs_frame_t> dfs_stack;
    dfs_frame_t start_vertex_frame;
    start_vertex_frame.vert = &graph_.get_vertex(start_vertex);
    start_vertex_frame.itr_adj = 
        (start_vertex_frame.vert->*itr_fn[reversed * 2])();
    start_vertex_frame.itr_adj_end =
        (start_vertex_frame.vert->*itr_fn[reversed * 2 + 1])();
    dfs_stack.push(start_vertex_frame);

    colors_[start_vertex] = node_color::gray;
    distances_[start_vertex] = 0;

    while (!dfs_stack.empty()) {
        auto& current_frame = dfs_stack.top();
        const int32_t src_vert = current_frame.vert->id;
        bool explore_further = false;

        while (current_frame.itr_adj != current_frame.itr_adj_end) {
            const int32_t dst_vert = current_frame.itr_adj->destination;

            //
            // Explore edge(src_vert, dst_vert)
            edge_explore_fn(src_vert, dst_vert);

            if (colors_[dst_vert] == node_color::white) {
                //
                // New node discovered
                colors_[dst_vert] = node_color::gray;
                distances_[dst_vert] = distances_[src_vert] + 1;
                predecessors_[dst_vert] = src_vert;

                //
                // Explore node (dst_vert)
                vx_explore_fn(src_vert, dst_vert);

                //
                // Save state, since we need to explore deeper.
                dfs_frame_t dst_frame;
                dst_frame.vert = &graph_.get_vertex(dst_vert);
                dst_frame.itr_adj = (dst_frame.vert->*itr_fn[reversed * 2])();
                dst_frame.itr_adj_end = 
                    (dst_frame.vert->*itr_fn[reversed * 2 + 1])();

                dfs_stack.push(dst_frame);
                ++current_frame.itr_adj;

                explore_further = true;
                break;
            }

            ++current_frame.itr_adj;
        }

        //
        // Start exploring from the last discovered vertex.
        if (explore_further)
            continue;

        //
        // All nodes from src_vert's adjacency list have been explored,
        // we mark src_vert as fully explored.
        colors_[src_vert] = node_color::black;
        //
        // Call finish function.
        vx_end_fn(src_vert);
        dfs_stack.pop();
    }
}

struct my_edge_fn {
    void operator()(int32_t src, int32_t dst) {
        printf("\n(%d ~> %d)", src, dst);
    }
};

struct my_vertex_done_fn {
    void operator()(int32_t src) {
        printf("\n(%d) explored.", src);
    }
};

struct topo_sort {
    std::vector<int32_t>    entry_times_;
    std::vector<int32_t>    exit_times_;
    std::list<int32_t>      sorted_vertices_;
    int32_t                 glob_time_;

    struct topo_vertex_discovery_fn {
        topo_vertex_discovery_fn(topo_sort* sstate)
            : sort_state_data(sstate) {}

        void operator()(int32_t /* src_vert */, int32_t dst_vert) const {
            sort_state_data->entry_times_[dst_vert] = 
                ++sort_state_data->glob_time_;
        }

        topo_sort* sort_state_data;
    };

    struct topo_vertex_end_fn {
        topo_vertex_end_fn(topo_sort* sstate)
            : sort_state_data(sstate) {}

        void operator()(int32_t vertex_id) const {
            sort_state_data->exit_times_[vertex_id] =
                ++sort_state_data->glob_time_;
            sort_state_data->sorted_vertices_.push_front(vertex_id);
        }

        topo_sort* sort_state_data;
    };

    topo_sort() : glob_time_(0) {}

    void initialize(graph& gr) {
        const auto num_vertices = gr.get_vertex_count();
        entry_times_.resize(num_vertices, -1);
        exit_times_.resize(num_vertices, -1);
        sorted_vertices_.clear();
        glob_time_ = 0;
    }

    void execute(graph& g) {
        dfs_search srch(g);

        srch.initialize();
        initialize(g);

        const int32_t num_vertices = (int32_t) g.get_vertex_count();
        for (int32_t i = 1; i <= num_vertices; ++i) {
            if (srch.get_color(i) == node_color::white) {
                entry_times_[i] = ++glob_time_;
                srch.execute(i, false, null_edge_discovery_fn(),
                             topo_vertex_discovery_fn(this),
                             topo_vertex_end_fn(this));
            }
        }
    }
};

void test_topo_sort() {
    const char* const kDataFilePath = "/home/adi.hodos/temp/activities.dat";

    graph g(true);
    if (!g.read_from_file(kDataFilePath)) {
        LAST_ERRNO(fopen);
        return;
    }

    topo_sort ts;
    ts.execute(g);

    printf("\n##### Topological sort : #####");
    std::for_each(ts.sorted_vertices_.begin(), ts.sorted_vertices_.end(),
                  [&ts](int32_t vertex_id) {
        printf("[(%d), (%d, %d)] -> ", vertex_id, ts.entry_times_[vertex_id],
               ts.exit_times_[vertex_id]);
    });
}

class scc_unidirected_graph {
public :
    //! Provides info about vertices in an SCC.
    struct scc_vertex {
        scc_vertex() 
            : scc_id(-1), vertex_id(-1) {}

        scc_vertex(
            int32_t scc, 
            int32_t vertex
            )
            :   scc_id(scc), 
                vertex_id(vertex) {}

        int32_t scc_id; //! SCC that the vertex belongs to
        int32_t vertex_id; //! Vertex id
    };


    scc_unidirected_graph() {}

    int32_t compute_sccs(
        graph& input_graph
        );

    const std::vector<scc_vertex>& get_scc_info() const {
        return sccs_;
    }

private :
    std::vector<scc_vertex>     sccs_;
};

int32_t scc_unidirected_graph::compute_sccs(
    graph& input_graph
    ) {
    assert(input_graph.is_directed());
    const auto num_verts = input_graph.get_vertex_count();

    struct dfs_data {

        struct postvisit_fn {
            postvisit_fn(dfs_data* state_p) : state_ptr(state_p) {}

            void operator()(int32_t vertex_id) const {
                state_ptr->vertex_list.push_back(vertex_id);
            }

            dfs_data*   state_ptr;
        };

        dfs_data(size_t num_verts) {
            vertex_list.reserve(num_verts);
        }
        
        std::vector<int32_t>      vertex_list;
    };

    dfs_search dfs(input_graph);
    dfs_data search_data(num_verts);

    auto null_pre_fn = [](int32_t, int32_t) {};

    dfs.initialize();
    for (size_t i = 1; i < num_verts; ++i) {
        if (dfs.get_color((int32_t) i) == node_color::white) {
            dfs.execute((int32_t) i, true, null_edge_discovery_fn(),
                        null_pre_fn, dfs_data::postvisit_fn(&search_data));
        }
    }

    int32_t scc_num = 0;
    auto itr_vert = search_data.vertex_list.rbegin();
    auto itr_end = search_data.vertex_list.rend();
    dfs.initialize();
    sccs_.reserve(num_verts);

    while (itr_vert != itr_end) {
        if (dfs.get_color(*itr_vert) == node_color::white) {
            ++scc_num;
            dfs.execute(*itr_vert, false, null_edge_discovery_fn(),
                        null_pre_fn, [this, scc_num](int32_t vertex_id) {
                sccs_.push_back({ scc_num, vertex_id });
            });
        }
        ++itr_vert;
    }

    return scc_num;
}

void test_scc() {
    graph g(true);

    const char* const kGraphDataFilePath = "/home/adi/temp/gd2.dat";
    if (!g.read_from_file(kGraphDataFilePath)) {
        LAST_ERRNO(fopen);
        return;
    }

    scc_unidirected_graph scc;
    int num_sccs = scc.compute_sccs(g);

    fprintf(stdout, "\nSCCS : %d", num_sccs);
    const auto& scc_lst = scc.get_scc_info();
    std::for_each(std::begin(scc_lst), std::end(scc_lst), 
                  [](const scc_unidirected_graph::scc_vertex& scc_component) {
        fprintf(stdout, "\n [%d - %d]", scc_component.scc_id, 
                scc_component.vertex_id);
    });
}

int main(int, char**) {
    //test_dfs();
    //test_topo_sort();
    test_scc();
    return 0;
}
