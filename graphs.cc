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
    std::vector<node_color> color_table_;
    distance_table_t        distances_;
    predecessor_table_t     predecessors_;
    
    void initialize() {
        const size_t vertex_count = graph_.get_vertex_count() + 1;
        color_table_.assign(vertex_count, node_color::white);
        distances_.assign(vertex_count, -1);
        predecessors_.assign(vertex_count, -1);
    }

    graph_search_base(graph& g) : graph_(g) {}

    ~graph_search_base() {}
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

class dfs_search : private graph_search_base {
public :
    using graph_search_base::distance_table_t;
    using graph_search_base::predecessor_table_t;
    using graph_search_base::get_distance_table;
    using graph_search_base::get_predecessors_table;

    dfs_search(graph& g) : graph_search_base(g) {}

    void execute(int32_t start_vertex, bool reversed = false);
};

void dfs_search::execute(int32_t start_vertex, bool reversed) {
    graph_search_base::initialize();

    struct stack_frame_t {
        const graph_vertex*                     vertex;
        graph_vertex::arc_list_const_iterator   adj_itr;
        graph_vertex::arc_list_const_iterator   adj_end;
    };
    
    typedef graph_vertex::arc_list_const_iterator 
        (graph_vertex::*const_itr_fn)() const;

    const_itr_fn itr_fn[4] = { 
        &graph_vertex::cbegin, &graph_vertex::cend,
        &graph_vertex::crbegin, &graph_vertex::crend,
    };

    std::stack<stack_frame_t> stk;
    stack_frame_t sf;
    sf.vertex = &graph_.get_vertex(start_vertex);
    sf.adj_itr = (sf.vertex->*itr_fn[reversed * 2])();
    sf.adj_end = (sf.vertex->*itr_fn[reversed * 2 + 1])();
    stk.push(sf);

    color_table_[start_vertex] = node_color::gray;
    distances_[start_vertex] = 0;

    while (!stk.empty()) {
        auto& current_frame = stk.top();
        const int32_t src_vert = current_frame.vertex->id;

        printf("\nVertex %d", src_vert);
        
        bool explore_further = false;
        while (current_frame.adj_itr != current_frame.adj_end) {
            printf("\n(%d ~> %d)", src_vert, current_frame.adj_itr->destination);

            if (color_table_[current_frame.adj_itr->destination] 
                == node_color::white) {
                const int32_t dst_vert = current_frame.adj_itr->destination;
                color_table_[dst_vert] = node_color::gray;
                distances_[dst_vert] = distances_[src_vert] + 1;
                predecessors_[dst_vert] = src_vert;

                stack_frame_t fv;
                fv.vertex = &graph_.get_vertex(dst_vert);
                fv.adj_itr = (fv.vertex->*itr_fn[reversed * 2])();
                fv.adj_end = (fv.vertex->*itr_fn[reversed * 2 + 1])();


                ++current_frame.adj_itr;

                stk.push(fv);
                explore_further = true;
                break;
            }

            ++current_frame.adj_itr;
        }

        if (explore_further)
            continue;

        color_table_[src_vert] = node_color::black;
        stk.pop();
    }
}

class topological_sort : private graph_search_base {
public :
    typedef std::vector<int32_t>    time_table_t;

private :
    time_table_t    entry_times_;
    time_table_t    exit_times_;

    void initialize() {
        graph_search_base::initialize();
        const auto vertex_cnt = graph_.get_vertex_count();
        entry_times_.assign(vertex_cnt, -1);
        exit_times_.assign(vertex_cnt, -1);
    }

public :
    using graph_search_base::distance_table_t;
    using graph_search_base::predecessor_table_t;
    using graph_search_base::get_distance_table;
    using graph_search_base::get_predecessors_table;

    topological_sort(graph& g) : graph_search_base(g) {}

    void execute(bool reversed = false);
};

void topological_sort::execute(bool reversed) {
    initialize();

    struct stack_frame_t {
        graph::vertex_const_iterator            itr_vertex;
        graph_vertex::arc_list_const_iterator   itr_adj_vertex;
        graph_vertex::arc_list_const_iterator   itr_adj_end;
    };
}

void test_dfs() {
    graph g(true);
    if (!g.read_from_file("/home/adihodos/temp/graph.dat"))
        return;

    dfs_search dfs(g);
    printf("\nExecuting dfs, start vertex = 1, direction forward");
    dfs.execute(1);
    printf("\nExecuting dfs, start vertex = 1, direction reversed");
    dfs.execute(1, true);
}

void test_bfs() {
    graph g(true);
    if (!g.read_from_file("/home/adihodos/temp/graph.dat"))
        return;

    bfs_search bfs(g);
    printf("\nExecuting bfs, start vertex = 1, direction forward");
    bfs.execute(1);
    printf("\nExecuting bfs, start vertex = 1, direction reversed");
    bfs.execute(1, true);
}

int main(int, char**) {
    test_bfs();
    test_dfs();
    return 0;
}
